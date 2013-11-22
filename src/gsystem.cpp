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
#include <sstream>
#include "gsystem.h"
#include "liegroup_rmatrix3_ext.h"

bool gTraceJointsBackward(GJoint *pEndJoint_, GBody *pLeftBodyOfEndJoint_, std::list<GJoint *> &pTracedJoints_);

//=============================================================
//                 GSystem
//=============================================================


GSystem::GSystem()
{
	pGround = NULL;
}


bool GSystem::setGravity(Vec3 g_)
{
	if ( pGround == NULL ) return false;

	pGround->dV = se3(Vec3(0,0,0),-g_);

	return true;
}

Vec3 GSystem::getGravity()
{
	return -pGround->dV.GetV();
}

GJoint* GSystem::getJoint(std::string name_)
{
	std::list<GJoint*>::iterator iter_pjoint;
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		if ( (*iter_pjoint)->getName() == name_ ) return *iter_pjoint;
	}
	return NULL;
}

GJoint* GSystem::getJoint(int idx_)
{
	if ( idx_ >= int(pJoints.size()) ) return NULL;
	std::list<GJoint*>::iterator iter_pjoint = pJoints.begin();
	std::advance(iter_pjoint, idx_); 
	return *iter_pjoint;
}

std::vector<GBody*> GSystem::getBodies()
{
	return std::vector<GBody*>(pBodies.begin(), pBodies.end());
}

std::vector<GJoint*> GSystem::getJoints()
{
	return std::vector<GJoint*>(pJoints.begin(), pJoints.end());
}

GBody* GSystem::getBody(std::string name_)
{
	std::list<GBody*>::iterator iter_pbody;
	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		if ( (*iter_pbody)->getName() == name_ ) return *iter_pbody;
	}
	return NULL;
}

GBody* GSystem::getBody(int idx_)
{
	if ( idx_ >= int(pBodies.size()) ) return NULL;
	std::list<GBody*>::iterator iter_pbody = pBodies.begin();
	std::advance(iter_pbody, idx_);
	return *iter_pbody;
}

int GSystem::getIndexJoint(std::string name_)
{
	int i;
	std::list<GJoint*>::iterator iter_pjoint;
	for (i=0, iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++, i++) {
		if ( (*iter_pjoint)->getName() == name_ ) return i;
	}
	return -1;
}

int GSystem::getIndexBody(std::string name_)
{
	int i;
	std::list<GBody*>::iterator iter_pbody;
	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++, i++) {
		if ( (*iter_pbody)->getName() == name_ ) return i;
	}
	return -1;
}

gReal GSystem::getMass()
{
	gReal mass = 0;
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		mass += (*iter_pbody)->getMass();
	}
	
	return mass;
}

Vec3 GSystem::getPositionCOMGlobal()
{
	Vec3 p(0,0,0);
	std::list<GBody *>::iterator iter_pbody;

	gReal mass_total = getMass();

	if ( mass_total > 1E-8 ) {
		for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
			p += (*iter_pbody)->getMass() * (*iter_pbody)->getPositionCOMGlobal();
		}
		p *= ((gReal)1./mass_total);
	}

	return p;
}

Vec3 GSystem::getVelocityCOMGlobal()
{
	Vec3 v(0,0,0);
	std::list<GBody *>::iterator iter_pbody;

	gReal mass_total = getMass();

	if ( mass_total > 1E-8 ) {
		for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
			v += (*iter_pbody)->getMass() * (*iter_pbody)->getVelocityCOMGlobal();
		}
		v *= ((gReal)1./mass_total);
	}

	return v;
}

Vec3 GSystem::getAccelerationCOMGlobal()
{
	Vec3 a(0,0,0);
	std::list<GBody *>::iterator iter_pbody;

	gReal mass_total = getMass();

	if ( mass_total > 1E-8 ) {
		for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
			a += (*iter_pbody)->getMass() * (*iter_pbody)->getAccelerationCOMGlobal();
		}
		a *= ((gReal)1./mass_total);
	}

	return a;
}

dse3 GSystem::getMomentumGlobal()
{
	dse3 H(0,0,0,0,0,0);
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		H += (*iter_pbody)->getMomentumGlobal();
	}

	return H;
}

dse3 GSystem::getMomentumCOM()
{
	return dAd(SE3(getPositionCOMGlobal()), getMomentumGlobal());
}

gReal GSystem::getGravitationalPotentialEnergy()
{
	gReal e = 0;
	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		e += Inner((*iter_pbody)->getPositionCOMGlobal(), -getGravity()) * (*iter_pbody)->getMass(); // e += <c,-mg>
	}
	return e;
}

gReal GSystem::getKineticEnergy()
{
	gReal e = 0;
	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		e += 0.5 * (*iter_pbody)->V.InnerProductWith( ((*iter_pbody)->I * (*iter_pbody)->V).GetArray() ); // e += 0.5*<V,I*V>
	}
	return e;
}

Vec3 GSystem::getDerivative_PositionCOMGlobal_Dq_2(GCoordinate *pCoordinate_)
{
	Vec3 DpDq(0,0,0), w, v, c;
	std::list<GBody *>::iterator iter_pbody;

	int idx = getIndexOfCoordinate(pCoordinate_);
	if ( idx < 0 ) return Vec3(0,0,0);

	gReal mass_total = getMass();

	if ( mass_total > 1E-8 ) {

		for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {

			w[0] = (*iter_pbody)->Jacobian(0,idx);  v[0] = (*iter_pbody)->Jacobian(3,idx);
			w[1] = (*iter_pbody)->Jacobian(1,idx);	v[1] = (*iter_pbody)->Jacobian(4,idx);
			w[2] = (*iter_pbody)->Jacobian(2,idx);	v[2] = (*iter_pbody)->Jacobian(5,idx);

			if ( w[0] != 0 || w[1] != 0 || w[2] != 0 || v[0] != 0 || v[1] != 0 || v[2] != 0 ) {
				c = (*iter_pbody)->getPositionCOM();
				DpDq += (*iter_pbody)->getMass() * ( (*iter_pbody)->T_global.GetRotation() * ( Cross(w, c) + v ) );
			} 
		}

		DpDq *= ((gReal)1./mass_total);

	}

	return DpDq;
}

Vec3 GSystem::getDerivative_PositionCOMGlobal_Dq(GCoordinate *pCoordinate_)
{
	Vec3 DpDq(0,0,0);
	std::list<GBody *>::iterator iter_pbody;

	gReal mass_total = getMass();

	if ( mass_total > 1E-8 ) {

		for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
			DpDq += (*iter_pbody)->getMass() * (*iter_pbody)->getDerivative_PositionCOMGlobal_Dq(pCoordinate_);
		}
		DpDq *= ((gReal)1./mass_total);

	}

	return DpDq;
}

dse3 GSystem::getDerivative_MomentumGlobal_Dq(GCoordinate *pCoordinate_)
{
	dse3 DHDq(0,0,0,0,0,0);
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		DHDq += (*iter_pbody)->getDerivative_MomentumGlobal_Dq(pCoordinate_);
	}

	return DHDq;
}

dse3 GSystem::getDerivative_MomentumGlobal_Ddq(GCoordinate *pCoordinate_)
{
	dse3 DHDdq(0,0,0,0,0,0);
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		DHDdq += (*iter_pbody)->getDerivative_MomentumGlobal_Ddq(pCoordinate_);
	}

	return DHDdq;
}

void GSystem::calcDerivative_PositionCOMGlobal_Dq(std::vector<GCoordinate *> pCoordinates_, RMatrix &DpDq_)
{
	int i;
	std::list<GBody *>::iterator iter_pbody;
	std::vector<GCoordinate *>::iterator iter_pcoord;
	std::vector<SE3> M(pBodies.size());
	std::vector<GConstraintJointLoop::JOINTLOOP_CONSTRAINT_TYPE> jlc_type(pBodies.size());
	Vec3 DpDqi;

	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {

		// save previous settings
		M[i] = (*iter_pbody)->fwdJointChain.M1;
		jlc_type[i] = (*iter_pbody)->fwdJointChain.jointLoopConstraintType;
		
		// prerequisites for (*iter_pbody)->getDerivative_PositionCOMGlobal_Dq(...)
		(*iter_pbody)->fwdJointChain.setM(SE3());
		(*iter_pbody)->fwdJointChain.jointLoopConstraintType = GConstraintJointLoop::JOINTLOOP_ORIENTATION_POSITION;
		(*iter_pbody)->fwdJointChain.update_J();
	}

	// calculate the derivative
	DpDq_.SetZero(3, int(pCoordinates_.size()));
	for (i=0, iter_pcoord = pCoordinates_.begin(); iter_pcoord != pCoordinates_.end(); i++, iter_pcoord++) {
		DpDqi = getDerivative_PositionCOMGlobal_Dq(*iter_pcoord);
		//DpDq_.Push(0, i, convert_to_RMatrix(DpDqi));
		matSet(&DpDq_[3*i], DpDqi.GetArray(), 3);
	}

	// restore the previous settings
	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {
		(*iter_pbody)->fwdJointChain.setM(M[i]);
		(*iter_pbody)->fwdJointChain.setJointLoopConstraintType(jlc_type[i]);
	}
}

void GSystem::calcDerivative_PositionCOMGlobal_Dq(RMatrix &DpDq_)
{
	int i;
	std::list<GBody *>::iterator iter_pbody;
	std::list<GCoordinate *>::iterator iter_pcoord;
	std::vector<SE3> M(pBodies.size());
	std::vector<GConstraintJointLoop::JOINTLOOP_CONSTRAINT_TYPE> jlc_type(pBodies.size());
	Vec3 DpDqi;

	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {

		// save previous settings
		M[i] = (*iter_pbody)->fwdJointChain.M1;
		jlc_type[i] = (*iter_pbody)->fwdJointChain.jointLoopConstraintType;
		
		// prerequisites for (*iter_pbody)->getDerivative_PositionCOMGlobal_Dq(...)
		(*iter_pbody)->fwdJointChain.setM(SE3());
		(*iter_pbody)->fwdJointChain.jointLoopConstraintType = GConstraintJointLoop::JOINTLOOP_ORIENTATION_POSITION;
		(*iter_pbody)->fwdJointChain.update_J();
	}

	// calculate the derivative
	DpDq_.SetZero(3, getNumCoordinates());
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); i++, iter_pcoord++) {
		DpDqi = getDerivative_PositionCOMGlobal_Dq(*iter_pcoord);
		//DpDq_.Push(0, i, convert_to_RMatrix(DpDqi));
		matSet(&DpDq_[3*i], DpDqi.GetArray(), 3);
	}

	// restore the previous settings
	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {
		(*iter_pbody)->fwdJointChain.setM(M[i]);
		(*iter_pbody)->fwdJointChain.setJointLoopConstraintType(jlc_type[i]);
	}
}

void GSystem::calcDerivative_PositionCOMGlobal_Dq_2(RMatrix &DpDq_)
{
	int i;
	std::list<GCoordinate *>::iterator iter_pcoord;
	Vec3 DpDqi;

	DpDq_.SetZero(3, getNumCoordinates());
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); i++, iter_pcoord++) {
		DpDqi = getDerivative_PositionCOMGlobal_Dq_2(*iter_pcoord);
		DpDq_(0,i) = DpDqi[0];
		DpDq_(1,i) = DpDqi[1];
		DpDq_(2,i) = DpDqi[2];
	}
}

void GSystem::calcDerivative_MomentumGlobal_Dq_Ddq(std::vector<GCoordinate*> pCoordinates_, RMatrix &DHgDq_, RMatrix &DHgDdq_)
{
	int i;
	std::list<GBody *>::iterator iter_pbody;
	std::vector<GCoordinate *>::iterator iter_pcoord;
	std::vector<SE3> M(pBodies.size());
	std::vector<GConstraintJointLoop::JOINTLOOP_CONSTRAINT_TYPE> jlc_type(pBodies.size());
	dse3 DHDqi, DHDdqi;

	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {

		// save previous settings
		M[i] = (*iter_pbody)->fwdJointChain.M1;
		jlc_type[i] = (*iter_pbody)->fwdJointChain.jointLoopConstraintType;
		
		// prerequisites for (*iter_pbody)->getDerivative_MomentumGlobal_Dq(...) and (*iter_pbody)->getDerivative_MomentumGlobal_Ddq(...)
		(*iter_pbody)->fwdJointChain.setM(SE3());
		(*iter_pbody)->fwdJointChain.jointLoopConstraintType = GConstraintJointLoop::JOINTLOOP_ORIENTATION_POSITION;
		(*iter_pbody)->fwdJointChain.update_J();
	}

	// calculate the derivatives
	DHgDq_.SetZero(6, int(pCoordinates_.size()));
	DHgDdq_.SetZero(6, int(pCoordinates_.size()));
	for (i=0, iter_pcoord = pCoordinates_.begin(); iter_pcoord != pCoordinates_.end(); i++, iter_pcoord++) {
		DHDqi = getDerivative_MomentumGlobal_Dq(*iter_pcoord);
		DHDdqi = getDerivative_MomentumGlobal_Ddq(*iter_pcoord);
		//DHgDq_.Push(0, i, convert_to_RMatrix(DHDqi));
		//DHgDdq_.Push(0, i, convert_to_RMatrix(DHDdqi));
		matSet(&DHgDq_[6*i], DHDqi.GetArray(), 6);
		matSet(&DHgDdq_[6*i], DHDdqi.GetArray(), 6);
	}

	// restore the previous settings
	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {
		(*iter_pbody)->fwdJointChain.setM(M[i]);
		(*iter_pbody)->fwdJointChain.setJointLoopConstraintType(jlc_type[i]);
	}
}

void GSystem::calcDerivative_MomentumGlobal_Dq_Ddq(RMatrix &DHgDq_, RMatrix &DHgDdq_)
{
	int i;
	std::list<GBody *>::iterator iter_pbody;
	std::list<GCoordinate *>::iterator iter_pcoord;
	std::vector<SE3> M(pBodies.size());
	std::vector<GConstraintJointLoop::JOINTLOOP_CONSTRAINT_TYPE> jlc_type(pBodies.size());
	dse3 DHDqi, DHDdqi;

	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {

		// save previous settings
		M[i] = (*iter_pbody)->fwdJointChain.M1;
		jlc_type[i] = (*iter_pbody)->fwdJointChain.jointLoopConstraintType;
		
		// prerequisites for (*iter_pbody)->getDerivative_MomentumGlobal_Dq(...) and (*iter_pbody)->getDerivative_MomentumGlobal_Ddq(...)
		(*iter_pbody)->fwdJointChain.setM(SE3());
		(*iter_pbody)->fwdJointChain.jointLoopConstraintType = GConstraintJointLoop::JOINTLOOP_ORIENTATION_POSITION;
		(*iter_pbody)->fwdJointChain.update_J();
	}

	// calculate the derivatives
	DHgDq_.SetZero(6, getNumCoordinates());
	DHgDdq_.SetZero(6, getNumCoordinates());
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); i++, iter_pcoord++) {
		DHDqi = getDerivative_MomentumGlobal_Dq(*iter_pcoord);
		DHDdqi = getDerivative_MomentumGlobal_Ddq(*iter_pcoord);
		//DHgDq_.Push(0, i, convert_to_RMatrix(DHDqi));
		//DHgDdq_.Push(0, i, convert_to_RMatrix(DHDdqi));
		matSet(&DHgDq_[6*i], DHDqi.GetArray(), 6);
		matSet(&DHgDdq_[6*i], DHDdqi.GetArray(), 6);
	}

	// restore the previous settings
	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {
		(*iter_pbody)->fwdJointChain.setM(M[i]);
		(*iter_pbody)->fwdJointChain.setJointLoopConstraintType(jlc_type[i]);
	}
}

void GSystem::calcDerivative_MomentumCOM_Dq_Ddq(std::vector<GCoordinate*> pCoordinates_, RMatrix &DHcDq_, RMatrix &DHcDdq_)
{
	int i;
	std::list<GBody *>::iterator iter_pbody;
	std::vector<GCoordinate *>::iterator iter_pcoord;
	std::vector<SE3> M(pBodies.size());
	std::vector<GConstraintJointLoop::JOINTLOOP_CONSTRAINT_TYPE> jlc_type(pBodies.size());

	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {

		// save previous settings
		M[i] = (*iter_pbody)->fwdJointChain.M1;
		jlc_type[i] = (*iter_pbody)->fwdJointChain.jointLoopConstraintType;
		
		// prerequisites for (*iter_pbody)->getDerivative_MomentumGlobal_Dq(...) and (*iter_pbody)->getDerivative_MomentumGlobal_Ddq(...)
		(*iter_pbody)->fwdJointChain.setM(SE3());
		(*iter_pbody)->fwdJointChain.jointLoopConstraintType = GConstraintJointLoop::JOINTLOOP_ORIENTATION_POSITION;
		(*iter_pbody)->fwdJointChain.update_J();
	}

	// calculate the derivatives
	Vec3 p = getPositionCOMGlobal();
	dse3 Hg = getMomentumGlobal();
	RMatrix DHgDq(6, int(pCoordinates_.size())), DHgDdq(6, int(pCoordinates_.size()));
	dse3 DHgDqi, DHgDdqi;
	for (i=0, iter_pcoord = pCoordinates_.begin(); iter_pcoord != pCoordinates_.end(); i++, iter_pcoord++) {
		DHgDqi = getDerivative_MomentumGlobal_Dq(*iter_pcoord);
		DHgDdqi = getDerivative_MomentumGlobal_Ddq(*iter_pcoord);
		//DHgDq.Push(0, i, convert_to_RMatrix(DHgDqi));
		//DHgDdq.Push(0, i, convert_to_RMatrix(DHgDdqi));
		matSet(&DHgDq[6*i], DHgDqi.GetArray(), 6);
		matSet(&DHgDdq[6*i], DHgDdqi.GetArray(), 6);
	}
	RMatrix DdAdDq_Hg(6, int(pCoordinates_.size()));
	Vec3 DpDqi;
	for (i=0, iter_pcoord = pCoordinates_.begin(); iter_pcoord != pCoordinates_.end(); i++, iter_pcoord++) {
		DpDqi = getDerivative_PositionCOMGlobal_Dq(*iter_pcoord);
		//DdAdDq_Hg.Push(0, i, convert_to_RMatrix(dAd(SE3(p), dad(se3(Vec3(0,0,0),DpDqi), Hg))));
		matSet(&DdAdDq_Hg[6*i], dAd(SE3(p), dad(se3(Vec3(0,0,0),DpDqi), Hg)).GetArray(), 6);
	}
	DHcDq_ = dAd(SE3(p), DHgDq) + DdAdDq_Hg;
	DHcDdq_ = dAd(SE3(p), DHgDdq);

	// restore the previous settings
	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {
		(*iter_pbody)->fwdJointChain.setM(M[i]);
		(*iter_pbody)->fwdJointChain.setJointLoopConstraintType(jlc_type[i]);
	}
}

void GSystem::calcDerivative_MomentumCOM_Dq_Ddq(RMatrix &DHcDq_, RMatrix &DHcDdq_)
{
	int i;
	std::list<GBody *>::iterator iter_pbody;
	std::list<GCoordinate *>::iterator iter_pcoord;
	std::vector<SE3> M(pBodies.size());
	std::vector<GConstraintJointLoop::JOINTLOOP_CONSTRAINT_TYPE> jlc_type(pBodies.size());

	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {

		// save previous settings
		M[i] = (*iter_pbody)->fwdJointChain.M1;
		jlc_type[i] = (*iter_pbody)->fwdJointChain.jointLoopConstraintType;
		
		// prerequisites for (*iter_pbody)->getDerivative_MomentumGlobal_Dq(...) and (*iter_pbody)->getDerivative_MomentumGlobal_Ddq(...)
		(*iter_pbody)->fwdJointChain.setM(SE3());
		(*iter_pbody)->fwdJointChain.jointLoopConstraintType = GConstraintJointLoop::JOINTLOOP_ORIENTATION_POSITION;
		(*iter_pbody)->fwdJointChain.update_J();
	}

	// calculate the derivatives
	Vec3 p = getPositionCOMGlobal();
	dse3 Hg = getMomentumGlobal();
	RMatrix DHgDq(6, getNumCoordinates()), DHgDdq(6, getNumCoordinates());
	dse3 DHgDqi, DHgDdqi;
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); i++, iter_pcoord++) {
		DHgDqi = getDerivative_MomentumGlobal_Dq(*iter_pcoord);
		DHgDdqi = getDerivative_MomentumGlobal_Ddq(*iter_pcoord);
		//DHgDq.Push(0, i, convert_to_RMatrix(DHgDqi));
		//DHgDdq.Push(0, i, convert_to_RMatrix(DHgDdqi));
		matSet(&DHgDq[6*i], DHgDqi.GetArray(), 6);
		matSet(&DHgDdq[6*i], DHgDdqi.GetArray(), 6);
	}
	RMatrix DdAdDq_Hg(6, getNumCoordinates());
	Vec3 DpDqi;
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); i++, iter_pcoord++) {
		DpDqi = getDerivative_PositionCOMGlobal_Dq(*iter_pcoord);
		//DdAdDq_Hg.Push(0, i, convert_to_RMatrix(dAd(SE3(p), dad(se3(Vec3(0,0,0),DpDqi), Hg))));
		matSet(&DdAdDq_Hg[6*i], dAd(SE3(p), dad(se3(Vec3(0,0,0),DpDqi), Hg)).GetArray(), 6);
	}
	DHcDq_ = dAd(SE3(p), DHgDq) + DdAdDq_Hg;
	DHcDdq_ = dAd(SE3(p), DHgDdq);

	// restore the previous settings
	for (i=0, iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); i++, iter_pbody++) {
		(*iter_pbody)->fwdJointChain.setM(M[i]);
		(*iter_pbody)->fwdJointChain.setJointLoopConstraintType(jlc_type[i]);
	}
}

void GSystem::setDeriv_Dq(GCoordinate *pCoordinate_)
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->set_bDpAlien(false);
		(*iter_pbody)->setDifferentiatingVariable_Dq(pCoordinate_);
	}
}

void GSystem::setDeriv_Ddq(GCoordinate *pCoordinate_)
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->set_bDpAlien(false);
		(*iter_pbody)->setDifferentiatingVariable_Ddq(pCoordinate_);
	}
}

void GSystem::setDeriv_Dddq(GCoordinate *pCoordinate_)
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->set_bDpAlien(false);
		(*iter_pbody)->setDifferentiatingVariable_Dddq(pCoordinate_);
	}
}

void GSystem::setDeriv_Dtau(GCoordinate *pCoordinate_)
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->set_bDpAlien(false);
		(*iter_pbody)->setDifferentiatingVariable_Dtau(pCoordinate_);
	}
}

void GSystem::setDeriv_DFe(GBody *pBody_, int idx_)
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->set_bDpAlien(false);
		(*iter_pbody)->setDifferentiatingVariable_DFe(pBody_, idx_);
	}
}

bool GSystem::buildSystem(GBody *pGround_, bool b_scan_fwd_joint_chain_)
{
	pCoordinates.clear();
	pBodies.clear();
	pJoints.clear();
	
	if ( pGround_ == NULL ) return false;

	pGround = pGround_;
	
	if ( !_scanBodiesAndJoints(pGround_) ) return false;

	pBodies.pop_front();	// pGround is popped out.

	if ( !_scanCoordinates() ) return false;

	if ( !_findParentChildRelation() ) return false;

	if ( b_scan_fwd_joint_chain_ ) {
		if ( !_findFwdJointChainOfBodies() ) return false;
	}

	if ( !_scanForceElements() ) return false;

	if ( !_getReady() ) return false;

	return true;
}

bool GSystem::buildSystemWith(GBody *pGround_, std::vector<GBody*> pFirstBodiesIn_, bool b_scan_fwd_joint_chain_)
{
	std::list<GJoint *>::iterator iter_pjoint;
	std::vector<GBody *>::iterator iter_pbody;
	std::vector<GBody *> childbodies;

	pCoordinates.clear();
	pBodies.clear();
	pJoints.clear();

	if ( pGround_ == NULL ) return false;

	pGround = pGround_;

	// find child bodies of the ground
	for (iter_pjoint = pGround_->pJoints.begin(); iter_pjoint != pGround_->pJoints.end(); iter_pjoint++) {
		if ( (*iter_pjoint)->pLeftBody == pGround_ && (*iter_pjoint)->pRightBody != pGround_ ) {
			childbodies.push_back((*iter_pjoint)->pRightBody);
		} else if ( (*iter_pjoint)->pLeftBody != pGround_ && (*iter_pjoint)->pRightBody == pGround_ ) {
			childbodies.push_back((*iter_pjoint)->pLeftBody);
		} else {
			;
		}
	}

	pBodies.push_back(pGround_); 

	// start scanning only from the child bodies which are listed in pFirstBodiesIn_
	for (iter_pbody = childbodies.begin(); iter_pbody != childbodies.end(); iter_pbody++) {
		if ( find(pFirstBodiesIn_.begin(), pFirstBodiesIn_.end(), *iter_pbody) != pFirstBodiesIn_.end() ) {
			bool bjointadded = false;
			for (iter_pjoint = (*iter_pbody)->pJoints.begin(); iter_pjoint != (*iter_pbody)->pJoints.end(); iter_pjoint++) {
				if ( (*iter_pjoint)->pLeftBody == pGround_ ) {
					pJoints.push_back(*iter_pjoint);
					bjointadded = true;
					break;
				} else if ( (*iter_pjoint)->pRightBody == pGround_ ) {
					if ( !(*iter_pjoint)->reverse() ) return false;
					pJoints.push_back(*iter_pjoint);
					bjointadded = true;
					break;
				} else {
					;
				}
			}
			if ( !bjointadded ) return false;
			if ( !_scanBodiesAndJoints(*iter_pbody) ) return false;
		}
	}

	pBodies.pop_front(); 

	if ( !_scanCoordinates() ) return false;

	if ( !_findParentChildRelation() ) return false;

	if ( b_scan_fwd_joint_chain_ ) {
		if ( !_findFwdJointChainOfBodies() ) return false;
	}

	if ( !_scanForceElements() ) return false;

	if ( !_getReady() ) return false;

	return true;
}

bool GSystem::buildSystemWithout(GBody *pGround_, std::vector<GBody*> pFirstBodiesOut_, bool b_scan_fwd_joint_chain_)
{
	std::list<GJoint *>::iterator iter_pjoint;
	std::vector<GBody *>::iterator iter_pbody;
	std::vector<GBody *> childbodies;

	pCoordinates.clear();
	pBodies.clear();
	pJoints.clear();

	if ( pGround_ == NULL ) return false;

	pGround = pGround_;

	// find child bodies of the ground
	for (iter_pjoint = pGround_->pJoints.begin(); iter_pjoint != pGround_->pJoints.end(); iter_pjoint++) {
		if ( (*iter_pjoint)->pLeftBody == pGround_ && (*iter_pjoint)->pRightBody != pGround_ ) {
			childbodies.push_back((*iter_pjoint)->pRightBody);
		} else if ( (*iter_pjoint)->pLeftBody != pGround_ && (*iter_pjoint)->pRightBody == pGround_ ) {
			childbodies.push_back((*iter_pjoint)->pLeftBody);
		} else {
			;
		}
	}

	pBodies.push_back(pGround_); 

	// start scanning only from the child bodies which are NOT listed in pFirstBodiesOut_
	for (iter_pbody = childbodies.begin(); iter_pbody != childbodies.end(); iter_pbody++) {
		if ( find(pFirstBodiesOut_.begin(), pFirstBodiesOut_.end(), *iter_pbody) == pFirstBodiesOut_.end() ) {
			bool bjointadded = false;
			for (iter_pjoint = (*iter_pbody)->pJoints.begin(); iter_pjoint != (*iter_pbody)->pJoints.end(); iter_pjoint++) {
				if ( (*iter_pjoint)->pLeftBody == pGround_ ) {
					pJoints.push_back(*iter_pjoint);
					bjointadded = true;
					break;
				} else if ( (*iter_pjoint)->pRightBody == pGround_ ) {
					if ( !(*iter_pjoint)->reverse() ) return false;
					pJoints.push_back(*iter_pjoint);
					bjointadded = true;
					break;
				} else {
					;
				}
			}
			if ( !bjointadded ) return false;
			if ( !_scanBodiesAndJoints(*iter_pbody) ) return false;
		}
	}

	pBodies.pop_front(); 

	if ( !_scanCoordinates() ) return false;

	if ( !_findParentChildRelation() ) return false;

	if ( b_scan_fwd_joint_chain_ ) {
		if ( !_findFwdJointChainOfBodies() ) return false;
	}

	if ( !_scanForceElements() ) return false;

	if ( !_getReady() ) return false;

	return true;
}

void GSystem::updateKinematics()
{
	update_joint_local_info();

	neFwdRecursion_a();

	update_joint_global_info();
}

void GSystem::updateSystemJacobianOfAllBodies()
{
	int idx;
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->Jacobian.SetZero(6, getNumCoordinates());
	}

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		if ( (*iter_pbody)->pBaseJoint->getNumCoordinates() > 0 ) {
			idx = getIndexOfCoordinate(*((*iter_pbody)->pBaseJoint->pCoordinates.begin()));
			if ( idx >= 0 ) {
				(*iter_pbody)->Jacobian.Push(0, idx, (*iter_pbody)->S);
				update_Jacobian_child_bodies(*iter_pbody, idx, (*iter_pbody)->S);
			}
		}
	}
}

void GSystem::update_Jacobian_child_bodies(GBody *pbody_, int idx_, const RMatrix &S_)
{
	std::list<GBody *>::iterator iter_pbody;
	RMatrix S2;
	for (iter_pbody = pbody_->pChildBodies.begin(); iter_pbody != pbody_->pChildBodies.end(); iter_pbody++) {
		S2 = Ad((*iter_pbody)->invT, S_);
		(*iter_pbody)->Jacobian.Push(0, idx_, S2);
		update_Jacobian_child_bodies(*iter_pbody, idx_, S2);	// recursive call until reaching end of branch
	}
}

void GSystem::calcDynamics()
{
	update_joint_local_info();

	fsFwdRecursion_a();
	fsBwdRecursion_b();
	fsFwdRecursion_c();
}

void GSystem::calcInverseDynamics()
{
	update_joint_local_info();

	neFwdRecursion_a();
	neBwdRecursion_b();
}

void GSystem::calcForwardDynamics()
{
	int i;
	std::list<GJoint *>::iterator iter_pjoint;
	std::vector<bool> isprescribed(pJoints.size());

	// save, and set all joints unprescribed
	for (i=0, iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); i++, iter_pjoint++) {
		isprescribed[i] = (*iter_pjoint)->isPrescribed();
	}
	setAllJointsPrescribed(false);

	calcDynamics();

	// restore
	for (i=0, iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); i++, iter_pjoint++) {
		(*iter_pjoint)->setPrescribed(isprescribed[i]);
	}
}

void GSystem::getMassMatrix(RMatrix &M)
{
	RMatrix ddq = get_ddq(), tau = get_tau(); // save current ddq and tau
	int n = getNumCoordinates();
	M.ReNew(n,n);
	set_ddq(Zeros(n,1));
	GSystem::calcInverseDynamics();
	RMatrix b0 = get_tau();
	for (int i=0; i<n; i++) {
		RMatrix unit = Zeros(n,1);
		unit[i] = 1;
		set_ddq(unit);
		GSystem::calcInverseDynamics();
		get_tau(&M[i*n]);
		for (int j=0; j<n; j++) {
			M[i*n+j] -= b0[j];
		}
	}
	set_ddq(ddq); set_tau(tau); // restore ddq and tau
}

void GSystem::getEquationsOfMotion(RMatrix &M, RMatrix &b)
{
	RMatrix ddq = get_ddq(), tau = get_tau(); // save current ddq and tau
	int n = getNumCoordinates();
	M.ReNew(n,n);
	set_ddq(Zeros(n,1));
	GSystem::calcInverseDynamics();
	b = get_tau();
	for (int i=0; i<n; i++) {
		RMatrix unit = Zeros(n,1);
		unit[i] = 1;
		set_ddq(unit);
		GSystem::calcInverseDynamics();
		get_tau(&M[i*n]);
		for (int j=0; j<n; j++) {
			M[i*n+j] -= b[j];
		}
	}
	set_ddq(ddq); set_tau(tau); // restore ddq and tau
}

bool GSystem::calcProductOfInvMassAndMatrix(RMatrix &invM_A, const RMatrix &A)
{
	if ( A.RowSize() != pCoordinates.size() ) return false;
	invM_A.SetZero(A.RowSize(), A.ColSize());

	int i;
	std::list<GJoint *>::iterator iter_pjoint;
	std::vector<bool> isprescribed(pJoints.size());

	// save current info
	for (i=0, iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); i++, iter_pjoint++) {
		isprescribed[i] = (*iter_pjoint)->isPrescribed();
	}
	Vec3 g = getGravity();

	// set all joint unprescribed and set zero gravity
	setAllJointsPrescribed(false);
	setGravity(Vec3(0,0,0));

	update_joint_local_info_short();

	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->update_base_joint_info();
		(*iter_pbody)->update_T();
		(*iter_pbody)->set_eta_zero();
	}

	for (std::list<GBody *>::reverse_iterator riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
		(*riter_pbody)->update_aI();
		(*riter_pbody)->update_Psi();
		(*riter_pbody)->update_Pi();
	}

	for (i=0; i<A.ColSize(); i++) {
		set_ddq(Zeros(pCoordinates.size(),1)); // this isn't necessary for real tree structure systems, but works for the cut joints in closed-loop
		set_tau(&(A[i*A.RowSize()]));
		
		for (std::list<GBody *>::reverse_iterator riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
			(*riter_pbody)->update_aB_zeroV_zeroeta();
			(*riter_pbody)->update_beta_zeroeta();
		}
		
		for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
			(*iter_pbody)->update_ddq();
			(*iter_pbody)->update_dV(true);
		}
		
		get_ddq(&(invM_A[i*invM_A.RowSize()]));
	}

	// restore
	for (i=0, iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); i++, iter_pjoint++) {
		(*iter_pjoint)->setPrescribed(isprescribed[i]);
	}
	setGravity(g);

	return true;
}

bool GSystem::calcProductOfInvMassAndMatrix2(RMatrix &invM_A, const RMatrix &A)
{
	if ( A.RowSize() != pCoordinates.size() ) return false;
	invM_A.SetZero(A.RowSize(), A.ColSize());

	int i;
	std::list<GJoint *>::iterator iter_pjoint;
	std::vector<bool> isprescribed(pJoints.size());

	// save current info
	for (i=0, iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); i++, iter_pjoint++) {
		isprescribed[i] = (*iter_pjoint)->isPrescribed();
	}

	setAllJointsPrescribed(false);

	initDynamicsWithZeroGravityAndVelocity();
	for (i=0; i<A.ColSize(); i++) {
		set_tau(&(A[i*A.RowSize()]));
		calcDynamicsWithZeroGravityAndVelocity();
		get_ddq(&(invM_A[i*invM_A.RowSize()]));
	}

	// restore
	for (i=0, iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); i++, iter_pjoint++) {
		(*iter_pjoint)->setPrescribed(isprescribed[i]);
	}

	return true;
}

void GSystem::initDynamicsWithZeroGravityAndVelocity()
{
	update_joint_local_info_short();

	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->update_base_joint_info();
		(*iter_pbody)->update_T();
		(*iter_pbody)->set_eta_zero();
	}
	
	for (std::list<GBody *>::reverse_iterator riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
		(*riter_pbody)->update_aI();
		(*riter_pbody)->update_Psi();
		(*riter_pbody)->update_Pi();
	}
}

void GSystem::calcDynamicsWithZeroGravityAndVelocity()
{
	Vec3 g = getGravity();
	setGravity(Vec3(0,0,0));

	double ddq[6];
	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->pBaseJoint->get_ddq(ddq);
		matSet_multAB((*iter_pbody)->Sddq.GetArray(), (*iter_pbody)->S.GetPtr(), ddq, 6, (*iter_pbody)->bjDOF, (*iter_pbody)->bjDOF, 1);
	}

	for (std::list<GBody *>::reverse_iterator riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
		(*riter_pbody)->update_aB_zeroV_zeroeta();
		(*riter_pbody)->update_beta_zeroeta();
	}
	
	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		if ( (*iter_pbody)->pBaseJoint->isPrescribed() ) {
			(*iter_pbody)->update_dV(false);
			(*iter_pbody)->update_F_fs();
			(*iter_pbody)->update_tau();
		} else {
			(*iter_pbody)->update_ddq();
			(*iter_pbody)->update_dV(true);
			//(*iter_pbody)->update_F_fs(); // we do not need this.
		}
	}

	setGravity(g);
}

//void GSystem::initInverseDynamicsWithZeroGravityAndVelocity()
//{
//	update_joint_local_info_short();
//
//	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
//		(*iter_pbody)->update_base_joint_info();
//		(*iter_pbody)->update_T();
//		(*iter_pbody)->set_eta_zero();
//	}
//}
//
//void GSystem::calcInverseDynamicsWithZeroGravityAndVelocity()
//{
//	Vec3 g = getGravity();
//	setGravity(Vec3(0,0,0));
//
//	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
//		(*iter_pbody)->update_dV(false);
//	}
//	for (std::list<GBody *>::reverse_iterator riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
//		(*riter_pbody)->update_F();
//		(*riter_pbody)->update_tau();
//	}
//
//	setGravity(g);
//}

void GSystem::diffDynamics()
{
	fsFwdRecursion_DaDp();
	fsBwdRecursion_DbDp();
	fsFwdRecursion_DcDp();
}

void GSystem::initBodyForces()
{
	for (std::list<GBody *>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->initExternalForce();
	}
}

void GSystem::initJointTorques()
{
	for (std::list<GCoordinate *>::iterator iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->tau = (gReal)0.0;
	}
}

void GSystem::initBodyForcesAndJointTorques()
{
	initBodyForces();
	initJointTorques();
}

bool GSystem::stepSimulation(gReal h_)
{
	// apply force elements to the bodies and joints
	// ** this will add forces/torques to the bodies and joints 
	// ** therefore, initBodyForcesAndJointTorques() must be called before calling stepSimulation() for initialization of the forces/torques
	for (std::list<GForce *>::iterator iter_psd = pForces.begin(); iter_psd != pForces.end(); iter_psd++) {
		(*iter_psd)->applyForce(); 
	}

	// calculate joint accelerations (for unprescribed coordinates) and torques (for prescribed coordinates)
	fsBwdRecursion_b();
	fsFwdRecursion_c();

	// mixed Euler integration
	for (std::list<GCoordinate *>::iterator iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dq += h_ * (*iter_pcoord)->ddq; // velocity update
		(*iter_pcoord)->q += h_ * (*iter_pcoord)->dq;   // position update
	}

	// update kinematics (updates positions and velocities of the bodies)
	update_joint_local_info();
	fsFwdRecursion_a();

	return true;
}

void GSystem::update_joint_local_info_short()
{
	std::list<GJoint *>::iterator iter_pjoint;

	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		(*iter_pjoint)->update_short();
	}
}

void GSystem::update_joint_local_info()
{
	std::list<GJoint *>::iterator iter_pjoint;

	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		(*iter_pjoint)->update();
	}
}

void GSystem::update_joint_global_info()
{
	// update GJoint::T_global using the connected body transform
	std::list<GJoint *>::iterator iter_pjoint;
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		(*iter_pjoint)->T_global = (*iter_pjoint)->pLeftBody->T_global;
		(*iter_pjoint)->T_global *= (*iter_pjoint)->T_left;
	}
}

void GSystem::neFwdRecursion_a()
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->neDynaRecursion_a();
	}
}

void GSystem::neBwdRecursion_b()
{
	std::list<GBody *>::reverse_iterator riter_pbody;

	for (riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
		(*riter_pbody)->neDynaRecursion_b();
	}
}

void GSystem::fsFwdRecursion_a()
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->fsDynaRecursion_a();
	}
}

void GSystem::fsBwdRecursion_b()
{
	std::list<GBody *>::reverse_iterator riter_pbody;

	for (riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
		(*riter_pbody)->fsDynaRecursion_b();
	}
}

void GSystem::fsFwdRecursion_c()
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->fsDynaRecursion_c();
	}
}

void GSystem::fsFwdRecursion_DaDp()
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->fsDynaRecursion_DaDp();
	}
}

void GSystem::fsBwdRecursion_DbDp()
{
	std::list<GBody *>::reverse_iterator riter_pbody;

	for (riter_pbody = pBodies.rbegin(); riter_pbody != pBodies.rend(); riter_pbody++) {
		(*riter_pbody)->fsDynaRecursion_DbDp();
	}
}

void GSystem::fsFwdRecursion_DcDp()
{
	std::list<GBody *>::iterator iter_pbody;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->fsDynaRecursion_DcDp();
	}

}

bool GSystem::_scanBodiesAndJoints(GBody *pbody_)
{
	GJoint *pjoint_;
	GBody *pbody_next_;
	std::list<GJoint*>::iterator iter_pjoint;
	
	if ( pbody_ == NULL ) return false;
	if ( find(pBodies.begin(), pBodies.end(), pbody_) != pBodies.end() ) return false;

	// add pbody_ to system
	pBodies.push_back(pbody_);
	
	for (iter_pjoint = pbody_->pJoints.begin(); iter_pjoint != pbody_->pJoints.end(); iter_pjoint++) {

		if ( (*iter_pjoint) == NULL ) return false;

		pjoint_ = *iter_pjoint;

		if ( find(pJoints.begin(), pJoints.end(), pjoint_) == pJoints.end() ) {
		
			// if needed, swap the bodies of pjoint_
			if ( pjoint_->pRightBody == pbody_ ) { 
				if ( !pjoint_->reverse() ) return false; 
			}

			// add pjoint_ to system
			pJoints.push_back(pjoint_);

			pbody_next_ = pjoint_->pRightBody;

			if ( find(pBodies.begin(), pBodies.end(), pbody_next_) == pBodies.end() ) {
				pjoint_->bCut = false;
				if ( !_scanBodiesAndJoints(pbody_next_) ) return false;
			} else {
				pjoint_->bCut = true;
			}
		}

	}

	return true;
}


bool GSystem::_scanCoordinates()
{
	std::list<GBody *>::iterator iter_pbody;
	std::list<GJoint *>::iterator iter_pjoint;
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		for (iter_pcoord = (*iter_pjoint)->pCoordinates.begin(); iter_pcoord != (*iter_pjoint)->pCoordinates.end(); iter_pcoord++) {
			if ( *iter_pcoord == NULL ) return false;
			if ( find(pCoordinates.begin(), pCoordinates.end(), *iter_pcoord) != pCoordinates.end() ) return false;
			pCoordinates.push_back(*iter_pcoord);
		}
	}

	return true;
}

bool GSystem::_findParentChildRelation()
{
	std::list<GBody *>::iterator iter_pbody;
	std::list<GJoint *>::iterator iter_pjoint;

	// find the base joint, parent body, and child bodies of pBodies
	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {

		(*iter_pbody)->pChildBodies.clear();

		for (iter_pjoint = (*iter_pbody)->pJoints.begin(); iter_pjoint != (*iter_pbody)->pJoints.end(); iter_pjoint++) {

			if ( !(*iter_pjoint)->isCut() ) {

				if ( (*iter_pjoint)->pRightBody == (*iter_pbody) ) {

					(*iter_pbody)->pBaseJoint = *iter_pjoint;

					(*iter_pbody)->pParentBody = (*iter_pjoint)->pLeftBody;

				} else if ( (*iter_pjoint)->pLeftBody == (*iter_pbody) ) {

					(*iter_pbody)->pChildBodies.push_back((*iter_pjoint)->pRightBody);

				} else {

					return false;
				}
			}
		}
	}

	return true;
}

bool GSystem::_findFwdJointChainOfBodies()
{
	std::list<GBody *>::iterator iter_pbody;
	std::list<GJoint *> pjoints;

	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		pjoints.clear();
		if ( !gTraceJointsBackward((*iter_pbody)->pBaseJoint, pGround, pjoints) ) return false;
		(*iter_pbody)->fwdJointChain.setJoints(pjoints);
	}

	return true;
}

bool GSystem::_scanForceElements()
{
	for (std::list<GBody*>::iterator iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		for (std::list<GForce*>::iterator iter_psd = (*iter_pbody)->pForces.begin(); iter_psd != (*iter_pbody)->pForces.end(); iter_psd++) {
			pForces.push_back(*iter_psd);
		}
	}
	for (std::list<GJoint*>::iterator iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		for (std::list<GForce*>::iterator iter_psd = (*iter_pjoint)->pForces.begin(); iter_psd != (*iter_pjoint)->pForces.end(); iter_psd++) {
			pForces.push_back(*iter_psd);
		}
	}
	pForces.unique();
	return true;
}

bool GSystem::_getReady()
{
	std::list<GJoint *>::iterator iter_pjoint;
	std::list<GBody *>::iterator iter_pbody;

	// check joints and bodies (This should be placed after _findParentChildRelation().)
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		if ( !(*iter_pjoint)->getReady() ) return false;
	}
	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		if ( !(*iter_pbody)->getReady() ) return false;
	}

	// update kinematics
	updateKinematics();

	// load joint information on each body (one time call)
	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		(*iter_pbody)->load_base_joint_info();	
	}

	return true;
}

std::string GSystem::getInfoStr()
{
	std::stringstream sstr;
	std::list<GBody *>::iterator iter_pbody;
	std::list<GJoint *>::iterator iter_pjoint;

	sstr << "GSystem:: " << GElement::getInfoStr();
	sstr << "    number of bodies included = " << int(pBodies.size()) << std::endl;
	sstr << "    number of joints included = " << int(pJoints.size()) << std::endl;
	sstr << "    number of coordinates = " << getNumCoordinates() << std::endl;
	sstr << "    ground body = " << pGround->getName() << std::endl;
	sstr << "--------------------------------------------------------------------" << std::endl;
	sstr << "         Bodies" << std::endl;
	sstr << "--------------------------------------------------------------------" << std::endl;
	for (iter_pbody = pBodies.begin(); iter_pbody != pBodies.end(); iter_pbody++) {
		sstr << (*iter_pbody)->getInfoStr();
		sstr << "----------" << std::endl;
	}
	sstr << "--------------------------------------------------------------------" << std::endl;
	sstr << "         Joints" << std::endl;
	sstr << "--------------------------------------------------------------------" << std::endl;
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		sstr << (*iter_pjoint)->getInfoStr();
		sstr << "----------" << std::endl;
	}

	return sstr.str();
}

void GSystem::setAllJointsPrescribed(bool b)
{
	for (std::list<GJoint*>::iterator iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		(*iter_pjoint)->setPrescribed(b);
	}
}

bool gTraceJointsBackward(GJoint *pEndJoint_, GBody *pLeftBodyOfEndJoint_, std::list<GJoint *> &pTracedJoints_)
{
	GJoint *pjoint;
	std::list<GJoint *>::iterator iter_pjoint;

	pTracedJoints_.clear();

	if ( pEndJoint_ == NULL ) return false;

	pjoint = pEndJoint_;

	while (1) {

		if ( pjoint->pLeftBody == pLeftBodyOfEndJoint_ ) { break; }

		for (iter_pjoint = pjoint->pLeftBody->pJoints.begin(); iter_pjoint != pjoint->pLeftBody->pJoints.end(); iter_pjoint++) {
			if ( *iter_pjoint != NULL && !(*iter_pjoint)->isCut() && (*iter_pjoint)->pRightBody == pjoint->pLeftBody ) {
				pTracedJoints_.push_front(*iter_pjoint);
				break;
			}
		}

		if ( iter_pjoint == pjoint->pLeftBody->pJoints.end() ) { return false; }

		pjoint = *iter_pjoint;
	}

	pTracedJoints_.push_back(pEndJoint_);

	return true;
}
