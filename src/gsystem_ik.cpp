// -------------------------------------------------------------------------------
// Copyright (c) 2012, Junggon Kim
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------------------------------------------------------------

#include <algorithm>
#include <fstream>
#include "gsystem_ik.h"
#include "rmatrix3j.h"


//=============================================================
//                 GSystemIK
//=============================================================

bool GSystemIK::buildConstrIK_dq(RMatrix &J, RMatrix &V, std::vector<GBody*> pbodies, std::vector<Vec3> pos, std::vector<se3> V_target, std::vector< std::vector<int> > idxC)
{
	int i, j, k;
	int cnt;			// a counter
	int nb;				// number of bodies
	int ncik;			// number of IK constraints
	std::list<GCoordinate*>::iterator iter_pcoord, iter_pcoord2;

	nb = int(pbodies.size());

	if ( pos.size() != nb || V_target.size() != nb || idxC.size() != nb ) return false;

	// counts the number of IK constraints
	ncik = 0;
	for (i=0; i<nb; i++) {
		for ( j=0; j<int(idxC[i].size()); j++) { 
			if ( idxC[i][j] < 0 || idxC[i][j] > 5 ) return false;
		}
		ncik += int(idxC[i].size());
	}
	J.SetZero(ncik, getNumCoordinates());
	V.SetZero(ncik, 1);

	// build J, V
	cnt = 0;
	for (i=0; i<nb; i++) {
		// update Jacobian
		pbodies[i]->fwdJointChain.setM(SE3(pos[i]));
		pbodies[i]->fwdJointChain.setJointLoopConstraintType(GConstraintJointLoop::JOINTLOOP_ORIENTATION_POSITION);
		pbodies[i]->fwdJointChain.update_J();

		// transformed Jacobian, Jg = [R 0; 0 R] * J
		RMatrix Jg(pbodies[i]->fwdJointChain.J.RowSize(), pbodies[i]->fwdJointChain.J.ColSize());
		//RMatrix R = convert_to_RMatrix(pbodies[i]->T_global.GetRotation());
		RMatrix R(3,3); matSet(R.GetPtr(), pbodies[i]->T_global.GetRotation().GetArray(), 9);
		Jg.Push(0, 0, R * pbodies[i]->fwdJointChain.J.Sub(0, 2, 0, pbodies[i]->fwdJointChain.J.ColSize()-1));
		Jg.Push(3, 0, R * pbodies[i]->fwdJointChain.J.Sub(3, 5, 0, pbodies[i]->fwdJointChain.J.ColSize()-1));

		// build J
		for (j=0, iter_pcoord = pbodies[i]->fwdJointChain.pCoordinates.begin(); iter_pcoord != pbodies[i]->fwdJointChain.pCoordinates.end(); j++, iter_pcoord++) {

			// find index of pbodies[i]->fwdJointChain.pCoordinates[j] in pCoordinates
			int idx = -1;
			for (k=0, iter_pcoord2 = pCoordinates.begin(); iter_pcoord2 != pCoordinates.end(); k++, iter_pcoord2++) {
				if ( *iter_pcoord2 == *iter_pcoord ) { idx = k; break; }
			}
			if ( idx < 0 ) return false;

			// insert j-th column of the body Jacobian to the right place
			for (k=0; k<int(idxC[i].size()); k++) {
				J(cnt+k, idx) = Jg(idxC[i][k], j);
			}
		}

		// build V
		for (j=0; j<int(idxC[i].size()); j++) {
			V(cnt+j, 0) = V_target[i][idxC[i][j]];
		}

		cnt += int(idxC[i].size());
	}

	return true;
}

bool GSystemIK::solveIK_dq(RMatrix &dq, std::vector<GBody*> pbodies_primary, std::vector<GBody*> pbodies_secondary, std::vector<Vec3> p_primary, std::vector<Vec3> p_secondary, std::vector<se3> V_primary, std::vector<se3> V_secondary, std::vector< std::vector<int> > idxC_primary, std::vector< std::vector<int> > idxC_secondary, gReal alpha_primary, gReal alpha_secondary)
{
	RMatrix Jp, Js;		// Jacobian matrices for primary/secondary constraints
	RMatrix Vp, Vs;		// the righthand side of the constraints

	if ( !buildConstrIK_dq(Jp, Vp, pbodies_primary, p_primary, V_primary, idxC_primary) ) return false;
	if ( !buildConstrIK_dq(Js, Vs, pbodies_secondary, p_secondary, V_secondary, idxC_secondary) ) return false;

	dq.SetZero(getNumCoordinates(), 1);
	if ( Jp.RowSize() > 0 ) {
		RMatrix dq0, N;
		dq0 = srInv(Jp, N, alpha_primary) * Vp;
		if ( Js.RowSize() > 0 ) {
			dq = dq0 + N * ( srInv(Js * N, alpha_secondary) * (Vs - Js * dq0) );
		} else {
			dq = dq0;
		}
	} else {
		if ( Js.RowSize() > 0 ) {
			dq = srInv(Js, alpha_secondary) * Vs;
		}
	}

	return true;
}

/*
bool GSystemIK::solveIK_dq(RMatrix &dq, vector<GBody*> pbodies_primary, vector<GBody*> pbodies_secondary, vector<Vec3> p_primary, vector<Vec3> p_secondary, vector<se3> V_primary, vector<se3> V_secondary, vector< vector<int> > idxC_primary, vector< vector<int> > idxC_secondary, vector<GCoordinate*> pcoords_disabled, gReal alpha)
{
	int idx;
	vector<GCoordinate*>::iterator viter_pcoord;
	std::list<GCoordinate*>::iterator iter_pcoord;
	RMatrix Jp, Js;		// Jacobian matrices for primary/secondary constraints
	RMatrix Vp, Vs;		// the righthand isde of the constraints

	if ( !buildConstrIK_dq(Jp, Vp, pbodies_primary, p_primary, V_primary, idxC_primary) ) return false;
	if ( !buildConstrIK_dq(Js, Vs, pbodies_secondary, p_secondary, V_secondary, idxC_secondary) ) return false;

	for (viter_pcoord = pcoords_disabled.begin(); viter_pcoord != pcoords_disabled.end(); viter_pcoord++) {
		for (idx = 0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); idx++, iter_pcoord++) {
			if ( *iter_pcoord == *viter_pcoord ) {
				for (int i=0; i<Jp.RowSize(); i++) { Jp(i,idx) = 0.0; }
				for (int i=0; i<Js.RowSize(); i++) { Js(i,idx) = 0.0; }
			}
		}
	}

	if ( Jp.RowSize() > 0 ) {
		RMatrix dq0, N;
		//dq0 = srInv(Jp, N, alpha) * Vp;
		dq0 = pInv(Jp, N) * Vp;
		dq = dq0 + N * ( srInv(Js * N, alpha) * (Vs - Js * dq0) );
	} else {
		dq = pInv(Js) * Vs;
	}

	return true;
}
*/

bool GSystemIK::solveIK_dq(RMatrix &dq, std::vector<GBody*> pbodies_primary, std::vector<GBody*> pbodies_secondary, std::vector<Vec3> p_primary, std::vector<Vec3> p_secondary, std::vector<se3> V_primary, std::vector<se3> V_secondary, std::vector< std::vector<int> > idxC_primary, std::vector< std::vector<int> > idxC_secondary, std::vector<GCoordinate*> pcoords_prescribed, gReal alpha_primary, gReal alpha_secondary)
{
	if ( pcoords_prescribed.size() == 0 ) {
		return solveIK_dq(dq, pbodies_primary, pbodies_secondary, p_primary, p_secondary, V_primary, V_secondary, idxC_primary, idxC_secondary, alpha_primary, alpha_secondary);
	}

	std::vector<GCoordinate*>::iterator viter_pcoord;
	std::list<GCoordinate*>::iterator iter_pcoord;
	RMatrix Jp, Js;		// Jacobian matrices for primary/secondary constraints
	RMatrix Vp, Vs;		// the righthand isde of the constraints

	if ( !buildConstrIK_dq(Jp, Vp, pbodies_primary, p_primary, V_primary, idxC_primary) ) return false;
	if ( !buildConstrIK_dq(Js, Vs, pbodies_secondary, p_secondary, V_secondary, idxC_secondary) ) return false;

	int num_prescribed = 0;
	std::vector<bool> b_prescribed;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		if ( find(pcoords_prescribed.begin(), pcoords_prescribed.end(), *iter_pcoord) == pcoords_prescribed.end() ) {
			b_prescribed.push_back(false);
		} else {
			b_prescribed.push_back(true);
			num_prescribed++;
		}
	}

	std::vector<GCoordinate*> pcoords_all(pCoordinates.begin(), pCoordinates.end());
	int nx = getNumCoordinates() - num_prescribed;
	RMatrix Jp_(Jp.RowSize(), nx), Js_(Js.RowSize(), nx);
	RMatrix Vp_(Vp), Vs_(Vs);

	int cnt = 0;
	for (int i=0; i<getNumCoordinates(); i++) {
		if ( b_prescribed[i] ) {
			Vp_ -= pcoords_all[i]->dq * Jp.Sub(0, Jp.RowSize()-1, i, i);
			Vs_ -= pcoords_all[i]->dq * Js.Sub(0, Js.RowSize()-1, i, i);
		} else {
			Jp_.Push(0, cnt, Jp.Sub(0, Jp.RowSize()-1, i, i));
			Js_.Push(0, cnt, Js.Sub(0, Js.RowSize()-1, i, i));
			cnt++;
		}
	}

	RMatrix x;
	x.SetZero(nx,1);
	if ( Jp_.RowSize() > 0 ) {
		RMatrix x0, N;
		x0 = srInv(Jp_, N, alpha_primary) * Vp_;
		if ( Js_.RowSize() > 0 ) {
			x = x0 + N * ( srInv(Js_ * N, alpha_secondary) * (Vs_ - Js_ * x0) );
		} else {
			x = x0;
		}
	} else {
		if ( Js_.RowSize() > 0 ) {
			x = srInv(Js_, alpha_secondary) * Vs_;
		}
	}

	dq.ReSize(getNumCoordinates(), 1);
	cnt = 0;
	for (int i=0; i<getNumCoordinates(); i++) {
		if ( b_prescribed[i] ) {
			dq[i] = pcoords_all[i]->dq;
		} else {
			dq[i] = x[cnt++];
		}
	}

	return true;
}

bool GSystemIK::solveIK_dq(RMatrix &dq, std::vector<GBody*> pbodies_primary, std::vector<GBody*> pbodies_secondary, std::vector<Vec3> p_primary, std::vector<Vec3> p_secondary, std::vector<se3> V_primary, std::vector<se3> V_secondary, std::vector< std::vector<int> > idxC_primary, std::vector< std::vector<int> > idxC_secondary, std::vector<GCoordinate*> pcoords_prescribed, std::ofstream *pfout, gReal alpha_primary, gReal alpha_secondary)
{
	if ( pcoords_prescribed.size() == 0 ) {
		return solveIK_dq(dq, pbodies_primary, pbodies_secondary, p_primary, p_secondary, V_primary, V_secondary, idxC_primary, idxC_secondary, alpha_primary, alpha_secondary);
	}

	std::vector<GCoordinate*>::iterator viter_pcoord;
	std::list<GCoordinate*>::iterator iter_pcoord;
	RMatrix Jp, Js;		// Jacobian matrices for primary/secondary constraints
	RMatrix Vp, Vs;		// the righthand isde of the constraints

	if ( !buildConstrIK_dq(Jp, Vp, pbodies_primary, p_primary, V_primary, idxC_primary) ) return false;
	if ( !buildConstrIK_dq(Js, Vs, pbodies_secondary, p_secondary, V_secondary, idxC_secondary) ) return false;

	int num_prescribed = 0;
	std::vector<bool> b_prescribed;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		if ( find(pcoords_prescribed.begin(), pcoords_prescribed.end(), *iter_pcoord) == pcoords_prescribed.end() ) {
			b_prescribed.push_back(false);
		} else {
			b_prescribed.push_back(true);
			num_prescribed++;
		}
	}

	std::vector<GCoordinate*> pcoords_all(pCoordinates.begin(), pCoordinates.end());
	int nx = getNumCoordinates() - num_prescribed;
	RMatrix Jp_(Jp.RowSize(), nx), Js_(Js.RowSize(), nx);
	RMatrix Vp_(Vp), Vs_(Vs);

	int cnt = 0;
	for (int i=0; i<getNumCoordinates(); i++) {
		if ( b_prescribed[i] ) {
			Vp_ -= pcoords_all[i]->dq * Jp.Sub(0, Jp.RowSize()-1, i, i);
			Vs_ -= pcoords_all[i]->dq * Js.Sub(0, Js.RowSize()-1, i, i);
		} else {
			Jp_.Push(0, cnt, Jp.Sub(0, Jp.RowSize()-1, i, i));
			Js_.Push(0, cnt, Js.Sub(0, Js.RowSize()-1, i, i));
			cnt++;
		}
	}

	RMatrix x;
	x.SetZero(nx,1);
	if ( Jp_.RowSize() > 0 ) {
		RMatrix x0, N;
		x0 = srInv(Jp_, N, alpha_primary) * Vp_;
		if ( Js_.RowSize() > 0 ) {
			x = x0 + N * ( srInv(Js_ * N, alpha_secondary) * (Vs_ - Js_ * x0) );
		} else {
			x = x0;
		}
	} else {
		if ( Js_.RowSize() > 0 ) {
			x = srInv(Js_, alpha_secondary) * Vs_;
		}
	}

	dq.ReSize(getNumCoordinates(), 1);
	cnt = 0;
	for (int i=0; i<getNumCoordinates(); i++) {
		if ( b_prescribed[i] ) {
			dq[i] = pcoords_all[i]->dq;
		} else {
			dq[i] = x[cnt++];
		}
	}

	if ( pfout != NULL ) {
		*pfout << "x = " << x << "Jp_ = " << Jp_ << "Js_ = " << Js_ << "Vp_ = " << Vp_ << "Vs_ = " << Vs_;
		*pfout << "Jp = " << Jp << "Js = " << Js << "Vp = " << Vp << "Vs = " << Vs;
		*pfout << "dq_prescribed = (";
		for (int i=0; i<(int)pcoords_prescribed.size(); i++) {
			*pfout << pcoords_prescribed[i]->dq << ", ";
		}
		*pfout << std::endl;
	}

	return true;
}
