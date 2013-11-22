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

#include <vector>
#include <algorithm>
#include "gsystem_constrained.h"
#include "gbody.h"
#include "gjoint.h"
#include "gcoordinate.h"
#include "gconstraint.h"
#include "gconstraint_jointloop.h"

using namespace std;

bool isZero(RMatrix A, gReal eps)
{
	bool iszero = true; 
	for (int i=0; i<A.Size(); i++) { if ( fabs(A[i]) > eps ) { iszero = false; break; } }
	return iszero;
}

//=============================================================
//                 GSystemConstrained
//=============================================================

bool GSystemConstrained::buildSystem(GBody *pGround_, bool b_scan_fwd_joint_chain_)
{
	pConstraints.clear();
	if ( !GSystem::buildSystem(pGround_, b_scan_fwd_joint_chain_) ) return false;
	if ( !_findClosedJointLoopConstraints() ) return false; 
	if ( !_getReady() ) return false;
	return true;
}

void GSystemConstrained::calcDynamics()
{
	vector<GCoordinate*> pcoordsU = getPrescribedCoordinates(), pcoordsV = getUnprescribedCoordinates(); // prescribed and unprescribed coordinates
	int n = getNumCoordinates(), nu = pcoordsU.size(), nv = pcoordsV.size(), nc = getNumConstraints();

	// if nc == 0, call unconstrained dynamics and return
	if ( nc == 0 ) { 
		GSystem::calcDynamics();
		_lambda = Zeros(0,0); 
		return; 
	}

	// memories for calculation
	RMatrix J(nc,n), dJdt(nc,n), Jv(nc,nv), Ju(nc,nu), invMvv_Jvt(nv,nc);
	RMatrix ddqu(nu,1), tauv(nv,1), ddqr(n,1), ddqvr(nv,1), ddqvc(nv,1), taur(n,1), tauur(nu,1), tauuc(nu,1), tauuc2(nu,1), dq = get_dq(), lambda(nc,1), Jvt_lambda;

	// saves the input of the hybrid dynamics problem (ddqu, tauv)
	for (int i=0; i<nu; i++) { ddqu[i] = pcoordsU[i]->ddq; }
	for (int i=0; i<nv; i++) { tauv[i] = pcoordsV[i]->tau; }

	// solves unconstrained dynamics and save the result: M * [ddqu; ddqvr] + b = [tauur; tauv]
	GSystem::calcDynamics();
	get_ddq(ddqr.GetPtr());
	get_tau(taur.GetPtr());
	for (int i=0; i<nv; i++) { ddqvr[i] = pcoordsV[i]->ddq; }
	for (int i=0; i<nu; i++) { tauur[i] = pcoordsU[i]->tau; }

	// calculates constraint Jacobian and its derivative
	_update_J(); 
	_copy_J(J, pCoordinates);
	_copy_J(Jv, pcoordsV);
	_copy_J(Ju, pcoordsU);
	_update_dJdt();
	_copy_dJdt(dJdt, pCoordinates);

	// calculates constraint forces (lambda)
	initDynamicsWithZeroGravityAndVelocity();
	invMvv_Jvt.SetZero();
	for (int i=0; i<nc; i++) {
		for (int j=0; j<nu; j++) { pcoordsU[j]->ddq = 0; }
		for (int j=0; j<nv; j++) { pcoordsV[j]->tau = Jv(i,j); }
		calcDynamicsWithZeroGravityAndVelocity();
		for (int j=0; j<nv; j++) { invMvv_Jvt(j,i) = pcoordsV[j]->ddq; }
	}
	solve_Ax_b_pInv(lambda, Jv*invMvv_Jvt, -(J*ddqr+dJdt*dq), tolerance_pinv);

	// solves constrained dynamics: M * [0; ddqvc] = [tauuc; ~Jv*lambda]
	Jvt_lambda = ~Jv*lambda;
	for (int i=0; i<nu; i++) { pcoordsU[i]->ddq = 0; }
	for (int i=0; i<nv; i++) { pcoordsV[i]->tau = Jvt_lambda[i]; }
	calcDynamicsWithZeroGravityAndVelocity();
	for (int i=0; i<nu; i++) { tauuc[i] = pcoordsU[i]->tau; }
	for (int i=0; i<nv; i++) { ddqvc[i] = pcoordsV[i]->ddq; }
	tauuc2 = tauuc - ~Ju*lambda;

	// final torque and acceleration
	for (int i=0; i<nu; i++) { pcoordsU[i]->ddq = ddqu[i]; pcoordsU[i]->tau = tauur[i] + tauuc2[i]; }
	for (int i=0; i<nv; i++) { pcoordsV[i]->tau = tauv[i]; pcoordsV[i]->ddq = ddqvr[i] + ddqvc[i]; }

	// save the constraint force
	_lambda = lambda; 
}

void GSystemConstrained::calcInverseDynamics()
{
	int n = getNumCoordinates(), nc = getNumConstraints();

	// if nc == 0, call unconstrained inverse dynamics and return
	if ( nc == 0 ) { 
		GSystem::calcInverseDynamics();
		_lambda = Zeros(0,0); 
		return; 
	}

	// memories for calculation
	RMatrix J(nc,n), Jt(n,nc), taur(n,1), tau(n,1), lambda(nc,1);

	// solve unconstrained inverse dynamics and save the result to taur
	GSystem::calcInverseDynamics();
	get_tau(taur.GetPtr());

	// calculates constraint Jacobian 
	_update_J(); 
	_copy_J(J, pCoordinates);
	Jt = ~J;

	// calculates constraint forces which minimizes the joint torque while satisfying the constraints
	solve_Ax_b_pInv(lambda, Jt, taur, tolerance_pinv);

	// solve constrained inverse dynamics
	tau = taur - Jt*lambda;

	// final torque
	set_tau(tau);

	// save the constraint force
	_lambda = lambda;
}

bool GSystemConstrained::stepSimulation(gReal h_)
{
	vector<GCoordinate*> pcoordsV = getUnprescribedCoordinates();

	// apply force elements to the bodies and joints
	// ** this will add forces/torques to the bodies and joints 
	// ** therefore, initBodyForcesAndJointTorques() must be called before calling stepSimulation() for initialization of the forces/torques
	for (std::list<GForce *>::iterator iter_psd = pForces.begin(); iter_psd != pForces.end(); iter_psd++) {
		(*iter_psd)->applyForce(); 
	}

	// calculate joint accelerations (for unprescribed coordinates) and torques (for prescribed coordinates)
	calcDynamics();

	// mixed Euler integration
	for (std::list<GCoordinate *>::iterator iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dq += h_ * (*iter_pcoord)->ddq; // velocity update
	}
	if ( !_adjust_dqv(pcoordsV) ) { cerr << "error:: failed in enforcing constraints in velocity level!" << endl; return false; } // adjust dqv
	for (std::list<GCoordinate *>::iterator iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->q += h_ * (*iter_pcoord)->dq;   // position update
	}
	if ( !_adjust_qv(pcoordsV) ) { cerr << "error:: failed in enforcing constraints in position level!" << endl; return false; } // adjust qv

	// update kinematics (updates positions and velocities of the bodies)
	update_joint_local_info();
	fsFwdRecursion_a();

	return true;
}

bool GSystemConstrained::addConstraint(GConstraint *pConstraint_)
{
	if ( pConstraint_ == NULL ) return false;
	if ( find(pConstraints.begin(), pConstraints.end(), pConstraint_) != pConstraints.end() ) return false;
	pConstraints.push_back(pConstraint_);
	return true;
}


bool GSystemConstrained::removeConstraint(GConstraint *pConstraint_)
{
	if ( pConstraint_ == NULL ) return false;
	if ( find(pConstraints.begin(), pConstraints.end(), pConstraint_) == pConstraints.end() ) return false;
	pConstraints.remove(pConstraint_);
	return true;
}


void GSystemConstrained::removeAllConstraints()
{
	pConstraints.clear();
}

bool GSystemConstrained::checkConstrainedKinematics(gReal eps)
{
	RMatrix C, J, dJdt, dq = get_dq(), ddq = get_ddq();
	updateKinematics();
	getConstraintFunction(C);
	getConstraintJacobian(J);
	getConstraintJacobianDerivative(dJdt);
	if ( !isZero(C, eps) ) { cerr << "C != 0" << endl; return false; }
	if ( !isZero(J*dq, eps) ) { cerr << "J*dq != 0" << endl; return false; }
	if ( !isZero(J*ddq+dJdt*dq, eps) ) { cerr << "J*ddq + dJdt*dq != 0" << endl; return false; }
	return true;
}

bool GSystemConstrained::checkConstrainedDynamics(gReal eps)
{
	RMatrix J, dJdt, M, b, dq = get_dq(), ddq = get_ddq(), tau = get_tau();	
	updateKinematics();
	getConstraintJacobian(J);
	getConstraintJacobianDerivative(dJdt);
	getEquationsOfMotion(M, b);
	if ( !isZero(M*ddq+b-tau-~J*_lambda, eps) ) { cerr << "M*ddq + b != tau + J'*lambda" << endl; return false; }
	if ( !isZero(J*ddq+dJdt*dq, eps) ) { cerr << "J*ddq + dJdt*dq != 0" << endl; return false; }
	return true;
}

bool GSystemConstrained::enforceConstraints(std::vector<GCoordinate*> pcoords)
{
	if ( !_adjust_qv(pcoords) ) return false;
	if ( !_adjust_dqv(pcoords) ) return false;
	if ( !_adjust_ddqv(pcoords) ) return false;
	return true;
}

//void GSystemConstrained::checkConstraintSolving(std::vector<GCoordinate*> pcoords)
//{
//	RMatrix Jv;
//	_update_J();
//	_copy_J(Jv, pcoords);
//	int rank = Rank(Jv), nc = getNumCoordinates();
//	if ( rank < nc ) {
//		cout << " Rank of Jacobian (" << rank << ") < number of constraints (" << nc << ")!" << endl;
//	}
//}

bool GSystemConstrained::_findClosedJointLoopConstraints()
{
	list<GJoint *> ploopjoints, pcutjoints;
	list<GJoint *>::iterator iter_pjoint;
	list<GJoint *>::iterator iter_pcutjoint;
	list<GCoordinate *>::iterator iter_pcoord;
	list<GConstraintJointLoop>::iterator iter_constr_jointloop;

	// find cut joints
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		if ( (*iter_pjoint)->isCut() ) pcutjoints.push_back(*iter_pjoint);
	}

	// set closedJointLoopConstraints
	closedJointLoopConstraints.clear();
	closedJointLoopConstraints.resize(pcutjoints.size());	// number of joint loops = number of cut joints

	for (iter_constr_jointloop = closedJointLoopConstraints.begin(), iter_pcutjoint = pcutjoints.begin(); iter_constr_jointloop != closedJointLoopConstraints.end(); iter_constr_jointloop++, iter_pcutjoint++) {
		// find joint loop for each cut joint
		if ( !_findClosedJointLoop(*iter_pcutjoint, ploopjoints) ) return false;

		// set joint loop constraints 
		(*iter_constr_jointloop).setJoints(ploopjoints);
		(*iter_constr_jointloop).setT(SE3());	// identity

		if ( !addConstraint(&(*iter_constr_jointloop)) ) return false;
	}

	// create dummy bodies and connect them with the cut joints
	for (iter_pcutjoint = pcutjoints.begin(); iter_pcutjoint != pcutjoints.end(); iter_pcutjoint++) {
		GBody *pcutbody = _createCutBody(*iter_pcutjoint);
		if ( !pcutbody ) return false;
		pCutBodies.push_back(pcutbody);
		pBodies.push_back(pcutbody);
	}

	return true;
}

 
bool GSystemConstrained::_findClosedJointLoop(GJoint *pCutJoint_, list<GJoint *> &loopjoints_)
{
	if ( !_findJointLoop(pCutJoint_, loopjoints_) ) return false;
	if ( (*loopjoints_.begin())->pLeftBody != pCutJoint_->pRightBody ) return false;
	return true;
}

 
bool GSystemConstrained::_findJointLoop(GJoint *pEndJoint_, list<GJoint *> &loopJoints_)
{
	GJoint *pjoint;
	list<GJoint *>::iterator iter_pjoint;

	if ( pEndJoint_ == NULL ) return false;

	loopJoints_.clear();

	pjoint = pEndJoint_;

	while (1) {

		for (iter_pjoint = pjoint->pLeftBody->pJoints.begin(); iter_pjoint != pjoint->pLeftBody->pJoints.end(); iter_pjoint++) {
			if ( *iter_pjoint != NULL && !(*iter_pjoint)->isCut() && (*iter_pjoint)->pRightBody == pjoint->pLeftBody && (*iter_pjoint) != pEndJoint_ ) {
				loopJoints_.push_front(*iter_pjoint);
				break;
			}
		}
		if ( (*iter_pjoint)->pLeftBody == pEndJoint_->pRightBody ) { break; }
		if ( (*iter_pjoint)->pLeftBody == pGround ) { break; }
		pjoint = *iter_pjoint;
	}

	loopJoints_.push_back(pEndJoint_);

	return true;
}

GBody* GSystemConstrained::_createCutBody(GJoint *pcutjoint)
{
	if ( !pcutjoint || !pcutjoint->isCut() || !pcutjoint->pRightBody ) return NULL;

	GBody *pbodyrightoriginal = pcutjoint->pRightBody;
	GBody *pbodyrightnew = new GBody; // create a dummy body
	pcutjoint->pRightBody = pbodyrightnew;
	pbodyrightnew->pJoints.push_back(pcutjoint);
	pbodyrightnew->pBaseJoint = pcutjoint;
	pbodyrightnew->pParentBody = pcutjoint->pLeftBody;
	pbodyrightnew->I = pbodyrightoriginal->I;
	if ( pbodyrightoriginal == pGround ) {
		pbodyrightnew->I = Inertia(1, 1, 1, 1);
	} else {
		// mass and inertia of the original body are shared by the two bodies
		pbodyrightnew->I *= 0.5;
		pbodyrightoriginal->I *= 0.5;
	}
	pbodyrightnew->setName("cut body");
	pcutjoint->pLeftBody->pChildBodies.push_back(pbodyrightnew);

	return pbodyrightnew;
}

bool GSystemConstrained::_adjust_qv(std::vector<GCoordinate*> pcoords)
{
	// -----------
	// adjusts qv by using Newton-Raphson method (while holding qu.)
	// C = C0 + Jv * del_qv = 0 where Jv = dC/dqv and del_qv = qv - qv0
	// -----------

	if ( getNumConstraints() == 0 ) return true;

	int cnt = 0;
	RMatrix C, Jv, del_qv((int)pcoords.size(),1);

	// Newton-Raphson iteration
	while (1) { 
		_update_C();
		_copy_C(C);
		
		// stop iteration if fabs(C[i]) is very small
		if ( isZero(C, tolerance_C) ) break;

		_update_J();
		_copy_J(Jv, pcoords);
		//if ( !SolveAxEqualB(Jv, del_qv, -C) ) return false; // solve C0 + Jv * del_qv = 0
		solve_Ax_b_pInv(del_qv, Jv, -C, tolerance_pinv); // solve C0 + Jv * del_qv = 0

		for (int i=0; i<(int)pcoords.size(); i++) {
			pcoords[i]->q += del_qv[i];
		}

		// perturb qv if del_qv is very small but fabs(C[i]) is not small
		if ( FNorm(del_qv) < tolerance_delq && !isZero(C, tolerance_Cmax) ) {
			for (int i=0; i<(int)pcoords.size(); i++) {
				pcoords[i]->q += drand(perturb_amount);
			}
		}
		
		// in case of too many iteration, stop iteration and return false
		if ( cnt++ > maxIterNum ) {
			return false;
		}
	}

	return true;
}

bool GSystemConstrained::_adjust_dqv(std::vector<GCoordinate*> pcoords, bool bupdate_J)
{
	// -----------
	// adjusts dqv to enforce the constraint: J * dq = 0 (while holding q and dqu.)
	// J * dq = Ju * dqu + Jv * (dqv0 + del_dqv) = J * dq0 + Jv * del_dqv = 0
	// -----------

	if ( getNumConstraints() == 0 ) return true;

	RMatrix dq0 = get_dq(), del_dqv((int)pcoords.size(),1), J, Jv; 
	if ( bupdate_J ) { 
		_update_J(); 
	}
	_copy_J(J, pCoordinates);
	_copy_J(Jv, pcoords);
	//if ( !SolveAxEqualB(Jv, del_dqv, -(J*dq)) ) return false; // solve J * dq0 + Jv * del_dqv = 0
	solve_Ax_b_pInv(del_dqv, Jv, -(J*dq0), tolerance_pinv); // solve J * dq0 + Jv * del_dqv = 0 
	for (int i=0; i<(int)pcoords.size(); i++) {
		pcoords[i]->dq += del_dqv[i];
	}

	return true;
}

bool GSystemConstrained::_adjust_ddqv(std::vector<GCoordinate*> pcoords, bool bupdate_J, bool bupdate_dJdt)
{
	// -----------
	// adjusts ddqv to enforce the constraint: J * ddq + dJdt * dq = 0  (while holding q, dq and ddqu.)
	// J * ddq + dJdt * dq = Ju * ddqu + Jv * (ddqv0 + del_ddqv) + dJdt * dq = J * ddq0 + Jv * del_ddqv + dJdt * dq = 0
	// -----------

	if ( getNumConstraints() == 0 ) return true;

	RMatrix dq = get_dq(), ddq0 = get_ddq(), del_ddqv((int)pcoords.size(),1), J, Jv, dJdt;
	if ( bupdate_J ) { 
		_update_J(); 
	}
	_copy_J(J, pCoordinates);
	_copy_J(Jv, pcoords);
	if ( bupdate_dJdt ) {
		_update_dJdt();
	}
	_copy_dJdt(dJdt, pCoordinates);
	//if ( !SolveAxEqualB(Jv, del_ddqv, -(J*ddq0+dJdt*dq)) ) return false; // solve J * ddq0 + Jv * del_ddqv + dJdt * dq = 0
	solve_Ax_b_pInv(del_ddqv, Jv, -(J*ddq0+dJdt*dq), tolerance_pinv); // solve J * ddq0 + Jv * del_ddqv + dJdt * dq = 0
	for (int i=0; i<(int)pcoords.size(); i++) {
		pcoords[i]->ddq += del_ddqv[i];
	}

	return true;
}

//// calculate equivalent independent coordinates force: (tau_u, tau_v) --> tau_u_equivalent
//bool GSystemConstrained::_calculateEquivalentIndependentCoordinatesForce(list<gReal> &tauu)
//{
//	int i, nu, nv;
//	static RMatrix tau_u, tau_v;
//	list<GCoordinate *>::iterator iter_pcoord;
//	
//	nu = getNumCoordinatesIndependent();
//	nv = getNumCoordinatesDependent();
//
//	tauu.clear();
//	
//	if ( nv == 0 ) {
//		for (iter_pcoord = pCoordinatesU.begin(); iter_pcoord != pCoordinatesU.end(); iter_pcoord++) {
//			tauu.push_back((*iter_pcoord)->tau);
//		}
//		return true;
//	}
//
//	tau_u.ReNew(nu, 1);
//	tau_v.ReNew(nv, 1);
//	
//	for (iter_pcoord = pCoordinatesU.begin(), i=0; iter_pcoord != pCoordinatesU.end(); iter_pcoord++, i++) {
//		tau_u[i] = (*iter_pcoord)->tau;
//	}
//	for (iter_pcoord = pCoordinatesV.begin(), i=0; iter_pcoord != pCoordinatesV.end(); iter_pcoord++, i++) {
//		tau_v[i] = (*iter_pcoord)->tau;
//	}
//
//	// tau_u = tau_u - Ju^T Jv^{-T} tau_v
//	tau_u -= Ju ^ ( Jv & tau_v );
//
//	for (iter_pcoord = pCoordinatesU.begin(), i=0; iter_pcoord != pCoordinatesU.end(); iter_pcoord++, i++) {
//		tauu.push_back(tau_u[i]);
//	}
//
//	return true;
//}

void GSystemConstrained::getConstraintFunction(RMatrix &C)
{
	_update_C();
	_copy_C(C);
}

void GSystemConstrained::getConstraintJacobian(RMatrix &J, std::vector<GCoordinate*> pcoords)
{
	_update_J();
	_copy_J(J, pcoords);
}

void GSystemConstrained::getConstraintJacobianDerivative(RMatrix &dJdt, std::vector<GCoordinate*> pcoords)
{
	_update_dJdt();
	_copy_dJdt(dJdt, pcoords);
}

void GSystemConstrained::_update_C()
{
	for (list<GConstraint *>::iterator iter_pconstr = pConstraints.begin(); iter_pconstr != pConstraints.end(); iter_pconstr++) {
		(*iter_pconstr)->update_C();
	}
}

void GSystemConstrained::_update_J()
{
	for (list<GConstraint *>::iterator iter_pconstr = pConstraints.begin(); iter_pconstr != pConstraints.end(); iter_pconstr++) {
		(*iter_pconstr)->update_J();
	}
}

void GSystemConstrained::_update_dJdt()
{
	for (list<GConstraint *>::iterator iter_pconstr = pConstraints.begin(); iter_pconstr != pConstraints.end(); iter_pconstr++) {
		(*iter_pconstr)->update_dJdt();
	}
}

void GSystemConstrained::_copy_C(RMatrix &C)
{
	int r = 0;

	C.SetZero(getNumConstraints(), 1);

	for (list<GConstraint *>::iterator iter_pconstr = pConstraints.begin(); iter_pconstr != pConstraints.end(); iter_pconstr++) {
		(*iter_pconstr)->get_C(&C[r]); //C.Push(r, 0, (*iter_pconstr)->get_C());
		r += (*iter_pconstr)->constrNum;
	}
}

void GSystemConstrained::_copy_J(RMatrix &J, vector<GCoordinate*> pcoords)
{
	int j, r=0, c, nc = getNumConstraints();
	list<GConstraint *>::iterator iter_pconstr;
	list<GCoordinate *>::iterator iter_pcoord_constr;
	vector<GCoordinate *>::iterator iter_pcoord;

	J.SetZero(nc, pcoords.size());

	for (iter_pconstr = pConstraints.begin(); iter_pconstr != pConstraints.end(); iter_pconstr++) {
		for (iter_pcoord_constr = (*iter_pconstr)->pCoordinates.begin(), j=0; iter_pcoord_constr != (*iter_pconstr)->pCoordinates.end(); iter_pcoord_constr++, j++) {
			for (iter_pcoord = pcoords.begin(), c=0; iter_pcoord != pcoords.end(); iter_pcoord++, c++) {
				if ( (*iter_pcoord) == (*iter_pcoord_constr) ) {
					(*iter_pconstr)->get_J(&J[r+c*nc], j); // J.Push(r, c, (*iter_pconstr)->get_J(j));
					break;
				}
			}
		}
		r += (*iter_pconstr)->constrNum;
	}
}

void GSystemConstrained::_copy_dJdt(RMatrix &dJdt, vector<GCoordinate*> pcoords)
{
	int j, r=0, c, nc = getNumConstraints();
	list<GConstraint *>::iterator iter_pconstr;
	list<GCoordinate *>::iterator iter_pcoord_constr;
	vector<GCoordinate *>::iterator iter_pcoord;

	dJdt.SetZero(nc, pcoords.size());

	for (iter_pconstr = pConstraints.begin(); iter_pconstr != pConstraints.end(); iter_pconstr++) {
		for (iter_pcoord_constr = (*iter_pconstr)->pCoordinates.begin(), j=0; iter_pcoord_constr != (*iter_pconstr)->pCoordinates.end(); iter_pcoord_constr++, j++) {
			for (iter_pcoord = pcoords.begin(), c=0; iter_pcoord != pcoords.end(); iter_pcoord++, c++) {
				if ( (*iter_pcoord) == (*iter_pcoord_constr) ) {
					(*iter_pconstr)->get_dJdt(&dJdt[r+c*nc], j); // dJdt.Push(r, c, (*iter_pconstr)->get_dJdt(j));
					break;
				}
			}
		}
		r += (*iter_pconstr)->constrNum;
	}
}
