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

#include <list>
#include <algorithm>
#include <string>
#include <sstream>

#include "gbody.h"
#include "gjoint.h"
#include "gcoordinate.h"

#include "liegroup.h"
#include "rmatrix3j.h"
#include "liegroup_rmatrix3_ext.h"

#define EPS_INV_SaIS 1E-12

//=============================================================
//                 GBody
//=============================================================
GBody::GBody()
{
	T_global.SetIdentity();
	pBaseJoint = NULL;
	pParentBody = NULL;

	T.SetIdentity();
	invT.SetIdentity();

	bjDOF = 0;
	S.SetZero(0,0);
	dS.SetZero(0,0);
	Sdq.SetZero();
	dSdq.SetZero();
	Sddq.SetZero();
	DSdqDt.SetZero();

	Jacobian.SetZero(0,0);

	I.SetZero();
	Fe.SetZero();

	V.SetZero();
	dV.SetZero();
	F.SetZero();
	aI.SetZero();
	aB.SetZero();
	eta.SetZero();
	aI_S.SetZero(0,0);
	Psi.SetZero(0,0);
	Pi.SetZero();
	beta.SetZero();

	DSDp.SetZero(0,0);
	DdSDp.SetZero(0,0);
	DhDp.SetZero();
	DFeDp.SetZero();
	DIDp.SetZero();

	bDpAlien = false;
	bnzDqDp = bnzDdqDp = false;
	bnzDSDp = bnzDdSDp = bnzDhDp = bnzDFeDp = bnzDIDp = false;

	DVDp.SetZero();
	DdVDp.SetZero();
	DFDp.SetZero();
	DaIDp.SetZero();
	DaBDp.SetZero();
	DetaDp.SetZero();
	DPsiDp.SetZero(0,0);
	DPiDp.SetZero();
	DbetaDp.SetZero();
}

void GBody::clear()
{
	GElement::clear();

	pJoints.clear();
	pForces.clear();
	pChildBodies.clear();

	T_global.SetIdentity();
	pBaseJoint = NULL;
	pParentBody = NULL;

	T.SetIdentity();
	invT.SetIdentity();

	bjDOF = 0;
	S.SetZero(0,0);
	dS.SetZero(0,0);
	Sdq.SetZero();
	dSdq.SetZero();
	Sddq.SetZero();
	DSdqDt.SetZero();

	Jacobian.SetZero(0,0);

	I.SetZero();
	Fe.SetZero();

	V.SetZero();
	dV.SetZero();
	F.SetZero();
	aI.SetZero();
	aB.SetZero();
	eta.SetZero();
	aI_S.SetZero(0,0);
	Psi.SetZero(0,0);
	Pi.SetZero();
	beta.SetZero();

	DSDp.SetZero(0,0);
	DdSDp.SetZero(0,0);
	DhDp.SetZero();
	DFeDp.SetZero();
	DIDp.SetZero();

	bDpAlien = false;
	bnzDqDp = bnzDdqDp = false;
	bnzDSDp = bnzDdSDp = bnzDhDp = bnzDFeDp = bnzDIDp = false;

	DVDp.SetZero();
	DdVDp.SetZero();
	DFDp.SetZero();
	DaIDp.SetZero();
	DaBDp.SetZero();
	DetaDp.SetZero();
	DPsiDp.SetZero(0,0);
	DPiDp.SetZero();
	DbetaDp.SetZero();
}

bool GBody::getReady()
{
	// check if the base joint of the body exists
	if ( pBaseJoint == NULL ) return false;
	bjDOF = pBaseJoint->getDOF();
	S.SetZero(6, bjDOF);
	dS.SetZero(6, bjDOF);
	aI_S.SetZero(6, bjDOF);
	Psi.SetZero(bjDOF, bjDOF);
	DSDp.SetZero(6, bjDOF);
	DdSDp.SetZero(6, bjDOF);
	DPsiDp.SetZero(bjDOF, bjDOF);
	return true;
}

void GBody::setMass(const gReal &mass_, const Vec3 &p_)
{
	setMass(mass_, 0, 0, 0, 0, 0, 0, SE3(p_));
}

void GBody::setMass(const gReal &mass_, const gReal &ixx_, const gReal &iyy_, const gReal &izz_, const gReal &ixy_, const gReal &ixz_, const gReal &iyz_, const SE3 &T_ref_)
{
	Inertia I_;						// I_ = the generalized inertia w.r.t. {ref}
	I_.SetMass(mass_);
	I_.SetInertia(ixx_,iyy_,izz_,ixy_,ixz_,iyz_);

	I = I_.Transform(Inv(T_ref_));	// I = the generalized inertia w.r.t. {body}
}

void GBody::addMass(const gReal &mass_, const gReal &ixx_, const gReal &iyy_, const gReal &izz_, const gReal &ixy_, const gReal &ixz_, const gReal &iyz_, const SE3 &T_ref_)
{
	Inertia I_;									// I_ = the generalized inertia w.r.t. {ref}
	I_.SetMass(mass_);
	I_.SetInertia(ixx_,iyy_,izz_,ixy_,ixz_,iyz_);
	
	Inertia I_b = I_.Transform(Inv(T_ref_));	// I_b = the generalized inertia w.r.t. {body}

	// I += I_b
	for (int i=0; i<6; i++) { I._I[i] += I_b._I[i]; }
	for (int i=0; i<3; i++) { I._r[i] += I_b._r[i]; }
	I._m += I_b._m;
}

void GBody::extractMass(const gReal &mass_, const gReal &ixx_, const gReal &iyy_, const gReal &izz_, const gReal &ixy_, const gReal &ixz_, const gReal &iyz_, const SE3 &T_ref_)
{
	Inertia I_;									// I_ = the generalized inertia w.r.t. {ref}
	I_.SetMass(mass_);
	I_.SetInertia(ixx_,iyy_,izz_,ixy_,ixz_,iyz_);

	Inertia I_b = I_.Transform(Inv(T_ref_));	// I_b = the generalized inertia w.r.t. {body}

	// I -= I_b
	for (int i=0; i<6; i++) { I._I[i] -= I_b._I[i]; }
	for (int i=0; i<3; i++) { I._r[i] -= I_b._r[i]; }
	I._m -= I_b._m;
}

void GBody::moveMass(const SE3 &T_ref_new_)
{
	Inertia I_ref_new = I;
	I = I_ref_new.Transform(Inv(T_ref_new_));
}

Vec3 GBody::getPositionGlobal()
{
	return T_global.GetPosition();
}

Vec3 GBody::getPositionGlobal(const Vec3 &p_)
{
	return T_global * p_;
}

SO3 GBody::getOrientationGlobal()
{
	return T_global.GetRotation();
}

SO3 GBody::getOrientationGlobal(const SO3 &R_)
{
	return T_global.GetRotation() * R_;
}

SE3 GBody::getPoseGlobal(const SE3 &T_)
{
	return T_global * T_;
}

Vec3 GBody::getVelocityLinearGlobal()
{
	return T_global.GetRotation() * V.GetV();
}

Vec3 GBody::getVelocityLinearGlobal(const Vec3 &p_)
{
	return T_global.GetRotation() * (Cross(V.GetW(),p_) + V.GetV());
}

Vec3 GBody::getVelocityAngularGlobal()
{
	return T_global.GetRotation() * V.GetW();
}

se3 GBody::getVelocityGlobal()
{
	SO3 R(T_global.GetRotation());
	return se3(R*V.GetW(), R*V.GetV());
}

se3 GBody::getVelocityGlobal(const Vec3 &p_)
{
	// Vp = Ad(SE3(R,0), Ad(Inv(SE3(Eye,p_)), V))
	//    = Ad(SE3(R,0) * Inv(SE3(Eye,p_)), V)
	//    = Ad(SE3(R, -R*p_), V)
	//    = se3(R*w, R*([w]*p + v))
	// where (w,v) = V, R = orientation of {body} w.r.t. {global}, and p_ = a position vector w.r.t. {body}

	SO3 R(T_global.GetRotation());
	Vec3 w(V.GetW()), v(V.GetV());
	return se3(R*w, R*(Cross(w,p_)+v));
}

Vec3 GBody::getAccelerationLinearGlobal()
{
	return T_global.GetRotation() * (Cross(V.GetW(), V.GetV()) + dV.GetV());
}

Vec3 GBody::getAccelerationLinearGlobal(const Vec3 &p_)
{
	SO3 R(T_global.GetRotation());
	Vec3 w(V.GetW()), v(V.GetV()), dw(dV.GetW()), dv(dV.GetV());
	return R*(Cross(w,Cross(w,p_)) + Cross(w,v) + Cross(dw,p_) + dv);
}

Vec3 GBody::getAccelerationAngularGlobal()
{
	return T_global.GetRotation() * dV.GetW();
}

se3 GBody::getAccelerationGlobal()
{
	SO3 R(T_global.GetRotation());
	Vec3 w(V.GetW()), v(V.GetV()), dw(dV.GetW()), dv(dV.GetV());
	return se3(R*dw, R*(Cross(w,v)+dv));
}

se3 GBody::getAccelerationGlobal(const Vec3 &p_)
{
	// Vp = Ad(T, V) where T = SE3(R, -R*p_)
	// d(Vp)/dt = ad(dT/dt*Inv(T), Ad(T, V)) + Ad(T, dV/dt) where dV/dt = GBody::dV
	//          = se3( R*dw, R*([w]*[w]*p_ + [w]*v + [dw]*p_ + dv ) where (w,v) = V, (dw,dv) = dV/dt

	SO3 R(T_global.GetRotation());
	Vec3 w(V.GetW()), v(V.GetV()), dw(dV.GetW()), dv(dV.GetV());
	return se3(R*dw, R*(Cross(w,Cross(w,p_)) + Cross(w,v) + Cross(dw,p_) + dv));
}

gReal GBody::getMass()
{
	return I.GetMass();
}

Vec3 GBody::getPositionCOM()
{
	if ( getMass() > 1E-8 ) {
		return ((gReal)1./getMass()) * I.GetOffDiag();
	} else {
		return Vec3(0,0,0);
	}
}

Vec3 GBody::getPositionCOMGlobal() 
{ 
	return T_global * getPositionCOM(); 
}

Vec3 GBody::getVelocityCOMGlobal()
{
	return getVelocityLinearGlobal(getPositionCOM());
}

Vec3 GBody::getAccelerationCOMGlobal()
{
	return getAccelerationLinearGlobal(getPositionCOM());
}

dse3 GBody::getMomentum()
{
	return  I*V;
}

dse3 GBody::getMomentumGlobal()
{
	return dAd(Inv(T_global), I*V);
}

Vec3 GBody::getDerivative_PositionCOMGlobal_Dq(GCoordinate *pCoordinate_)
{
	if ( find(fwdJointChain.pCoordinates.begin(), fwdJointChain.pCoordinates.end(), pCoordinate_) == fwdJointChain.pCoordinates.end() ) return Vec3(0,0,0);
	se3 Ji; fwdJointChain.get_J(Ji.GetArray(), pCoordinate_); //se3 Ji = convert_to_se3(fwdJointChain.get_J(pCoordinate_));
	return T_global.GetRotation() * ( Cross(Ji.GetW(), getPositionCOM()) + Ji.GetV() );
}

dse3 GBody::getDerivative_MomentumGlobal_Dq(GCoordinate *pCoordinate_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	se3 DVDqi, Ji;

	if ( find(fwdJointChain.pCoordinates.begin(), fwdJointChain.pCoordinates.end(), pCoordinate_) == fwdJointChain.pCoordinates.end() ) return dse3(0,0,0,0,0,0);

	fwdJointChain.get_J(Ji.GetArray(), pCoordinate_); // Ji = convert_to_se3(fwdJointChain.get_J(pCoordinate_));

	DVDqi.SetZero();
	for (iter_pcoord = fwdJointChain.pCoordinates.begin(); iter_pcoord != fwdJointChain.pCoordinates.end(); iter_pcoord++) {
		if ( *iter_pcoord == pCoordinate_ ) break;
		//DVDqi += ad( convert_to_se3(fwdJointChain.get_J(*iter_pcoord)), convert_to_se3(fwdJointChain.get_J(pCoordinate_)) ) * (*iter_pcoord)->dq;
		DVDqi += ad( fwdJointChain.getPtr_J(*iter_pcoord), fwdJointChain.getPtr_J(pCoordinate_) ) * (*iter_pcoord)->dq;
	}

	return dAd(Inv(T_global), I*DVDqi - dad(Ji, I*V));
}

dse3 GBody::getDerivative_MomentumGlobal_Ddq(GCoordinate *pCoordinate_)
{
	//return dAd(Inv(T_global), I * convert_to_se3(fwdJointChain.get_J(pCoordinate_)));
	return dAd(Inv(T_global), I * fwdJointChain.getPtr_J(pCoordinate_));
}

void GBody::setExternalForceLocally(const dse3 &Fe_local_)
{
	Fe = Fe_local_;
}

void GBody::setExternalForceGlobally(const dse3 &Fe_global_)
{
	Fe = dAd(T_global, Fe_global_);
}

void GBody::setExternalForceGlobally(const Vec3 &p_, const Vec3 &fg_)
{
	Vec3 f = ~getOrientationGlobal() * fg_;
	Fe = dse3(Cross(p_, f), f);
}

void GBody::addExternalForceLocally(const dse3 &Fe_local_)
{
	Fe += Fe_local_;
}

void GBody::addExternalForceGlobally(const dse3 &Fe_global_)
{
	Fe += dAd(T_global, Fe_global_);
}

void GBody::addExternalForceGlobally(const Vec3 &p_, const Vec3 &fg_)
{
	Vec3 f = ~getOrientationGlobal() * fg_;
	Fe += dse3(Cross(p_, f), f);
}

void GBody::neDynaRecursion_a()
{
	update_base_joint_info();

	update_T();
	update_V();				
	update_eta();
	update_dV(false);
}

void GBody::neDynaRecursion_b()
{
	update_F();
	update_tau();
}

void GBody::fsDynaRecursion_a()
{
	update_base_joint_info();

	update_T();
	update_V();				
	update_eta();
}

void GBody::fsDynaRecursion_b()
{
	update_aI();
	update_aB();
	update_Psi();
	update_Pi();
	update_beta();
}

void GBody::fsDynaRecursion_c()
{
	if ( pBaseJoint->isPrescribed() ) {
		update_dV(false);
		update_F_fs();
		update_tau();
	} else {
		update_ddq();
		update_dV(true);
		update_F_fs();
	}
}

void GBody::neDynaRecursion_DaDp()
{
	update_DVDp();				
	update_DetaDp();
	update_DdVDp();
}

void GBody::neDynaRecursion_DbDp()
{
	update_DFDp();
	update_DtauDp();
}

void GBody::fsDynaRecursion_DaDp()
{
	update_DVDp();				
	update_DetaDp();
}

void GBody::fsDynaRecursion_DbDp()
{
	update_DaIDp();
	update_DaBDp();
	update_DPsiDp();
	update_DPiDp();
	update_DbetaDp();
}

void GBody::fsDynaRecursion_DcDp()
{
	if ( pBaseJoint->isPrescribed() ) {
		update_DdVDp();
		update_DFDp_fs();
		update_DtauDp();
	} else {
		update_DddqDp();
		update_DdVDp();
		update_DFDp_fs();
	}
}

dse3 GBody::getTransformed_F()
{
	return dAd(invT, F);
}

AInertia GBody::getTransformed_aI()
{
	return Pi.Transform(invT);
}

dse3 GBody::getTransformed_aB()
{
	return dAd(invT, beta);
}

dse3 GBody::getTransformed_DFDp()
{
	dse3 tmp(DFDp);

	if ( bnzDhDp ) {
		tmp -= dad(DhDp, F);
	}
	
	return dAd(invT, tmp);
}

AInertia GBody::getTransformed_DaIDp()
{
	AInertia tmp(DPiDp);

	if ( bnzDhDp ) {
		tmp.SubstractTransform_ad(Pi, DhDp); // tmp -= Pi*ad(DhDp) + ~(Pi*ad(DhDp))
	}

	return tmp.Transform(invT);
}

dse3 GBody::getTransformed_DaBDp()
{
	dse3 tmp(DbetaDp);

	if ( bnzDhDp ) {
		tmp -= dad(DhDp, beta);
	}

	return dAd(invT, tmp);
}

void GBody::setDifferentiatingVariable_Dq(GCoordinate *pCoordinate_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (iter_pcoord = pBaseJoint->pCoordinates.begin(); iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DqDp = 0.0;
		(*iter_pcoord)->DdqDp = 0.0;
		if ( pBaseJoint->isPrescribed() ) {
			(*iter_pcoord)->DddqDp = 0.0;
		} else {
			(*iter_pcoord)->DtauDp = 0.0;	////////////// in case of spring, this is wrong!
		}
	}
	set_bnzAll(false);

	// if pCoordinate_ is included in pBaseJoint->pCoordinates
	if ( isIncluded(pCoordinate_) ) {

		pCoordinate_->DqDp = 1.0;
		set_bnzDqDp(true);

		if ( !(pBaseJoint->isConstantScrew()) ) {
			set_DSDp(get_DSDq(pCoordinate_));
			set_DdSDp(get_DdSDq(pCoordinate_));
			set_bnzDSDp(true);
			set_bnzDdSDp(true);
		}

		set_DhDp(get_S(pCoordinate_));
		set_bnzDhDp(true);		// DhDp = S*DqDp
	}
}

void GBody::setDifferentiatingVariable_Ddq(GCoordinate *pCoordinate_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (iter_pcoord = pBaseJoint->pCoordinates.begin(); iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DqDp = 0.0;
		(*iter_pcoord)->DdqDp = 0.0;
		if ( pBaseJoint->isPrescribed() ) {
			(*iter_pcoord)->DddqDp = 0.0;
		} else {
			(*iter_pcoord)->DtauDp = 0.0;	////////////// in case of spring, this is wrong!
		}
	}
	set_bnzAll(false);

	// if pCoordinate_ is included in pBaseJoint->pCoordinates
	if ( isIncluded(pCoordinate_) ) {
		pCoordinate_->DdqDp = 1.0;
		set_bnzDdqDp(true);

		if ( !(pBaseJoint->isConstantScrew()) ) {
			set_DdSDp(get_DdSDdq(pCoordinate_));
			set_bnzDdSDp(true);
		}
	}
}

void GBody::setDifferentiatingVariable_Dddq(GCoordinate *pCoordinate_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (iter_pcoord = pBaseJoint->pCoordinates.begin(); iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DqDp = 0.0;
		(*iter_pcoord)->DdqDp = 0.0;
		if ( pBaseJoint->isPrescribed() ) {
			(*iter_pcoord)->DddqDp = 0.0;
		} else {
			(*iter_pcoord)->DtauDp = 0.0;
		}
	}
	set_bnzAll(false);

	set_bDpAlien(true);

	// if pCoordinate_ is included in pBaseJoint->pCoordinates
	if ( isIncluded(pCoordinate_) ) {
		if ( pBaseJoint->isPrescribed() ) {
			pCoordinate_->DddqDp = 1.0;
		}
	}
}

void GBody::setDifferentiatingVariable_Dtau(GCoordinate *pCoordinate_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (iter_pcoord = pBaseJoint->pCoordinates.begin(); iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DqDp = 0.0;
		(*iter_pcoord)->DdqDp = 0.0;
		if ( pBaseJoint->isPrescribed() ) {
			(*iter_pcoord)->DddqDp = 0.0;
		} else {
			(*iter_pcoord)->DtauDp = 0.0;
		}
	}
	set_bnzAll(false);

	set_bDpAlien(true);

	// if pCoordinate_ is included in pBaseJoint->pCoordinates
	if ( isIncluded(pCoordinate_) ) {
		if ( !pBaseJoint->isPrescribed() ) {
			pCoordinate_->DtauDp = 1.0;
		}
	}
}

void GBody::setDifferentiatingVariable_DFe(GBody *pBody_, int idx_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (iter_pcoord = pBaseJoint->pCoordinates.begin(); iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DqDp = 0.0;
		(*iter_pcoord)->DdqDp = 0.0;
		if ( pBaseJoint->isPrescribed() ) {
			(*iter_pcoord)->DddqDp = 0.0;
		} else {
			(*iter_pcoord)->DtauDp = 0.0;
		}
	}
	set_bnzAll(false);

	set_bDpAlien(true);

	if ( pBody_ == this && idx_ >= 0 && idx_ < 6 ) {
		bnzDFeDp = true;
		DFeDp.SetZero();
		DFeDp[idx_] = 1.0;
	}
}

void GBody::load_base_joint_info()
{
	Sdq.set_Ad(pBaseJoint->T_right, pBaseJoint->Sdq);
	dSdq.set_Ad(pBaseJoint->T_right, pBaseJoint->dSdq);
	Sddq.set_Ad(pBaseJoint->T_right, pBaseJoint->Sddq);
	DSdqDt = Sddq;
	DSdqDt += dSdq;
	matSet_Ad(S.GetPtr(), pBaseJoint->T_right, pBaseJoint->S.GetPtr(), pBaseJoint->getDOF());
	matSet_Ad(dS.GetPtr(), pBaseJoint->T_right, pBaseJoint->dS.GetPtr(), pBaseJoint->getDOF());
}

void GBody::update_base_joint_info()
{
	Sdq.set_Ad(pBaseJoint->T_right, pBaseJoint->Sdq);
	dSdq.set_Ad(pBaseJoint->T_right, pBaseJoint->dSdq);
	Sddq.set_Ad(pBaseJoint->T_right, pBaseJoint->Sddq);
	DSdqDt = Sddq;
	DSdqDt += dSdq;

	if ( !(pBaseJoint->isConstantScrew()) ) {
		matSet_Ad(S.GetPtr(), pBaseJoint->T_right, pBaseJoint->S.GetPtr(), pBaseJoint->getDOF());
		matSet_Ad(dS.GetPtr(), pBaseJoint->T_right, pBaseJoint->dS.GetPtr(), pBaseJoint->getDOF());
	}
}

void GBody::update_T()
{
	T = pBaseJoint->T_left;
	T *= pBaseJoint->T;
	T *= pBaseJoint->inv_T_right;

	invT.SetInvOf(T);
	
	T_global = pParentBody->T_global * T; // update global location of the body
}

void GBody::update_V()
{
	if ( bjDOF <= 0 ) {
		V.set_Ad(invT, pParentBody->V);
	} else {
		V.set_Ad(invT, pParentBody->V);
		V += Sdq;
	}
}

void GBody::update_eta()
{
	if ( bjDOF <= 0 ) {
		eta.SetZero();
	} else {
		eta.set_ad(V, Sdq);
		eta += dSdq;
	}
}

void GBody::update_dV(bool b_update_)
{
	gReal ddq_[6], Sddq_[6];

	if ( bjDOF <= 0 ) {
		dV.set_Ad(invT, pParentBody->dV);
	} else {
		dV.set_Ad(invT, pParentBody->dV);
		if ( b_update_ ) { // use new joint acceleration 
			pBaseJoint->get_ddq(ddq_);
			matSet_multAB(Sddq_, S.GetPtr(), ddq_, 6, bjDOF, bjDOF, 1);
			dV += Sddq_;
		} else {
			dV += Sddq;
		}
		dV += eta;
	}
}

void GBody::update_F()
{
	std::list<GBody *>::iterator iter_pbody_child;

	F = I * dV;			// inertial force
	F -= dad(V, I * V);	// Coriolis force
	F -= Fe;			// external force

	// force from child bodies
	for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
		F += (*iter_pbody_child)->getTransformed_F();
	}
}

void GBody::update_F_fs()
{
	F = aI * dV + aB;
}

void GBody::update_aI()
{
	std::list<GBody *>::iterator iter_pbody_child;

	aI = I; // aI = AInertia(I);

	for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
		aI += (*iter_pbody_child)->getTransformed_aI();
		//aI.AddTransform(((GBody*)(*iter_pbody_child))->Pi, ((GBody*)(*iter_pbody_child))->invT); 
	}
}

void GBody::update_aB()
{
	std::list<GBody *>::iterator iter_pbody_child;

	aB = -dad(V, I*V);
	aB -= Fe;

	for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
		aB += (*iter_pbody_child)->getTransformed_aB();
	}
}

void GBody::update_aB_zeroV_zeroeta()
{
	aB.SetZero();
	for (std::list<GBody *>::iterator iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
		aB += (*iter_pbody_child)->getTransformed_aB();
	}
}

void GBody::update_Psi()
{
	gReal tmp[36];

	if ( pBaseJoint->isPrescribed() ) {
		;
	} else {
		if ( bjDOF > 0 ) {
			// Psi = Inv(~S * aI * S)
			set_Mult_AInertia_se3(aI_S.GetPtr(), aI, S.GetPtr(), bjDOF); // aI_S = aI * S
			matSet_multAtB(tmp, S.GetPtr(), aI_S.GetPtr(), 6, bjDOF, 6, bjDOF);	// tmp = ~S * aI_S
			// inverse of ~S * aI * S
			// Note: If matrix inverse fails (bsuccess == false), then Psi will not be updated. (Psi keeps its previous values.)
			bool bsuccess = true;
			switch ( bjDOF ) {
				case 1:
					if ( fabs(tmp[0]) > EPS_INV_SaIS ) {
						Psi[0] = (gReal)1./tmp[0];
					} else {
						bsuccess = false;
					}
					break;
				case 2:
					bsuccess = matSet_inv22sym(Psi.GetPtr(), tmp);
					break;
				case 3:
					bsuccess = matSet_inv33sym(Psi.GetPtr(), tmp);
					break;
				case 4:
					bsuccess = matSet_inv44sym(Psi.GetPtr(), tmp);
					break;
				default:
					bsuccess = matSet_invNNfast(Psi.GetPtr(), tmp, bjDOF);
			}
		}
	}
}

void GBody::update_Pi()
{
	gReal tmp1[36], tmp2[36];

	if ( pBaseJoint->isPrescribed() ) {
		Pi = aI;
	} else {
		// Pi = aI - aI*S*Psi*~S*~aI, aI is symmetric.
		Pi = aI;

		if ( bjDOF <= 3 ) {
			for (int i=0; i<bjDOF; i++) {
				for (int j=0; j<bjDOF; j++) {
					Pi.SubstractAlphaXYt(Psi[i+bjDOF*j], &aI_S[6*i], &aI_S[6*j]);
				}
			}
		} else {
			matSet_multABt(tmp1, Psi.GetPtr(), aI_S.GetPtr(), bjDOF, bjDOF, 6, bjDOF);
			matSet_multAB(tmp2, aI_S.GetPtr(), tmp1, 6, bjDOF, bjDOF, 6);
			Pi -= tmp2;
		}
	}
}

void GBody::update_beta()
{
	gReal tau_[6], Psi_tau_[6], S_Psi_tau_[6]; // with asuming bjDOF <= 6

	if ( pBaseJoint->isPrescribed() ) {
		if ( bjDOF <= 0 ) {
			beta = aB + aI * eta;
		} else {
			beta = aB + aI * (eta + Sddq);
		}
	} else {
		if ( bjDOF <= 0 ) {
			beta = aB + aI * eta;
		} else {
			//beta = aB + aI * ( eta + convert_to_se3( S * Psi * ( pBaseJoint->get_tau() - ~S * convert_to_RMatrix( aI * eta + aB ) ) ) );
			dse3 aI_eta_aB_ = aI * eta + aB;
			pBaseJoint->get_tau(tau_);
			for (int i=0; i<bjDOF; i++) { tau_[i] -= aI_eta_aB_.InnerProductWith(&S[6*i]); } //this is faster than "matSet_multAB(tau2_, aI_eta_aB_.GetArray(), S.GetPtr(), 1, 6, 6, bjDOF); for (int i=0; i<bjDOF; i++) { tau_[i] -= tau2_[i]; }"
			matSet_multAB(Psi_tau_, Psi.GetPtr(), tau_, bjDOF, bjDOF, bjDOF, 1);
			matSet_multAB(S_Psi_tau_, S.GetPtr(), Psi_tau_, 6, bjDOF, bjDOF, 1);
			beta = aI_eta_aB_ + aI * se3(S_Psi_tau_);
		}
	}
}

void GBody::update_beta_zeroeta()
{
	gReal tau_[6], Psi_tau_[6], S_Psi_tau_[6]; // with asuming bjDOF <= 6

	if ( pBaseJoint->isPrescribed() ) {
		if ( bjDOF <= 0 ) {
			beta = aB;
		} else {
			beta = aB + aI * Sddq;
		}
	} else {
		if ( bjDOF <= 0 ) {
			beta = aB;
		} else {
			//beta = aB + aI * ( convert_to_se3( S * Psi * ( pBaseJoint->get_tau() - ~S * convert_to_RMatrix( aB ) ) ) );
			pBaseJoint->get_tau(tau_);
			for (int i=0; i<bjDOF; i++) { tau_[i] -= aB.InnerProductWith(&S[6*i]); } //this is faster than "matSet_multAB(tau2_, aB.GetArray(), S.GetPtr(), 1, 6, 6, bjDOF); for (int i=0; i<bjDOF; i++) { tau_[i] -= tau2_[i]; }"
			matSet_multAB(Psi_tau_, Psi.GetPtr(), tau_, bjDOF, bjDOF, bjDOF, 1);
			matSet_multAB(S_Psi_tau_, S.GetPtr(), Psi_tau_, 6, bjDOF, bjDOF, 1);
			beta = aB + aI * se3(S_Psi_tau_);
		}
	}
}

void GBody::update_tau()
{
	gReal tau_[6];

	int i;
	std::list<GCoordinate *>::iterator iter_pcoord;

	if ( bjDOF <= 0 ) return;

	//RMatrix tau_ = ~S * convert_to_RMatrix(F);
	for (int i=0; i<bjDOF; i++) { tau_[i] = F.InnerProductWith(&S[6*i]); }

	for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
		(*iter_pcoord)->tau = tau_[i];
	}
}

void GBody::update_ddq()
{
	gReal tau_[6], ddq_[6];

	int i;
	std::list<GCoordinate *>::iterator iter_pcoord;

	if ( bjDOF <= 0 ) return;

	//RMatrix ddq_ = Psi * ( pBaseJoint->get_tau() - ~S * convert_to_RMatrix( aI * ( Ad(invT, pParentBody->dV) + eta ) ) - ~S * convert_to_RMatrix(aB) );
	se3 eta_ = Ad(invT, pParentBody->dV) + eta;
	dse3 aB_ = aI * eta_ + aB;
	pBaseJoint->get_tau(tau_);
	for (int i=0; i<bjDOF; i++) { tau_[i] -= aB_.InnerProductWith(&S[6*i]); }
	matSet_multAB(ddq_, Psi.GetPtr(), tau_, bjDOF, bjDOF, bjDOF, 1);

	for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
		(*iter_pcoord)->ddq = ddq_[i];
	}
}

void GBody::set_eta_zero()
{
	eta.SetZero();
}

void GBody::update_DVDp()
{
	gReal dq[6], DSDp_dq[6];
	gReal DdqDp[6], S_DdqDp[6];

	if ( bDpAlien ) {
		DVDp.SetZero();
		return;
	}

	DVDp.set_Ad(invT, pParentBody->DVDp);
	
	if ( bnzDhDp ) {
		DVDp -= ad(DhDp, Ad(invT, pParentBody->V));
	}

	if ( bnzDSDp ) {
		//DVDp += convert_to_se3(DSDp * pBaseJoint->get_dq());
		pBaseJoint->get_dq(dq);
		matSet_multAB(DSDp_dq, DSDp.GetPtr(), dq, 6, bjDOF, bjDOF, 1);
		DVDp += DSDp_dq;
	}

	if ( bnzDdqDp ) {
		//DVDp += convert_to_se3(S * pBaseJoint->get_DdqDp());
		pBaseJoint->get_DdqDp(DdqDp);
		matSet_multAB(S_DdqDp, S.GetPtr(), DdqDp, 6, bjDOF, bjDOF, 1);
		DVDp += S_DdqDp;
	}
}

void GBody::update_DetaDp()
{
	gReal dq[6], DdqDp[6], tmp1[6], tmp2[6];
	pBaseJoint->get_dq(dq);
	pBaseJoint->get_DdqDp(DdqDp);
	//RMatrix tmp = Zeros(6,1);
	for (int i=0; i<6; i++) { tmp2[i] = 0.0; }

	if ( bDpAlien ) {
		DetaDp.SetZero();
		return;
	}

	DetaDp.set_ad(DVDp, Sdq);

	if ( bnzDSDp && bnzDdqDp ) {
		//tmp += ad(V, DSDp * pBaseJoint->get_dq() + S * pBaseJoint->get_DdqDp());
		matSet_multAB(tmp1, DSDp.GetPtr(), dq, 6, bjDOF, bjDOF, 1);	// tmp1 = DSDp * dq
		matAdd_multAB(tmp1, S.GetPtr(), DdqDp, 6, bjDOF, bjDOF, 1);	// tmp1 += S * DdqDp
		matAdd_ad(tmp2, V.GetArray(), tmp1);						// tmp2 += ad(V,tmp1)
	} else if ( bnzDSDp && !bnzDdqDp ) {
		//tmp += ad(V, DSDp * pBaseJoint->get_dq());
		matSet_multAB(tmp1, DSDp.GetPtr(), dq, 6, bjDOF, bjDOF, 1);	// tmp1 = DSDp * dq
		matAdd_ad(tmp2, V.GetArray(), tmp1);						// tmp2 += ad(V,tmp1)
	} else if ( !bnzDSDp && bnzDdqDp ) {
		//tmp += ad(V, S * pBaseJoint->get_DdqDp());
		matSet_multAB(tmp1, S.GetPtr(), DdqDp, 6, bjDOF, bjDOF, 1);	// tmp1 = S * DdqDp
		matAdd_ad(tmp2, V.GetArray(), tmp1);						// tmp2 += ad(V,tmp1)
	} else {
		;
	}

	if ( bnzDdSDp ) {
		//tmp += DdSDp * pBaseJoint->get_dq();
		matAdd_multAB(tmp2, DdSDp.GetPtr(), dq, 6, bjDOF, bjDOF, 1);	// tmp2 += DdSDp * dq
	}

	if ( bnzDdqDp ) {
		//tmp += dS * pBaseJoint->get_DdqDp();
		matAdd_multAB(tmp2, dS.GetPtr(), DdqDp, 6, bjDOF, bjDOF, 1);	// tmp2 += dS * DdqDp
	}

	//DetaDp += convert_to_se3(tmp);
	DetaDp += tmp2;
}

void GBody::update_DdVDp()
{
	gReal DddqDp[6], ddq[6];

	if ( bDpAlien ) {
		DdVDp.set_Ad(invT, pParentBody->DdVDp);
		//DdVDp += convert_to_se3(S * pBaseJoint->get_DddqDp());
		pBaseJoint->get_DddqDp(DddqDp);
		matAdd_multAB(DdVDp.GetArray(), S.GetPtr(), DddqDp, 6, bjDOF, bjDOF, 1);
		return;
	}

	DdVDp.set_Ad(invT, pParentBody->DdVDp);

	if ( bnzDhDp ) {
		DdVDp -= ad(DhDp, Ad(invT, pParentBody->dV));
	}

	if ( bnzDSDp ) {
		//DdVDp += convert_to_se3(DSDp * pBaseJoint->get_ddq());
		pBaseJoint->get_ddq(ddq);
		matAdd_multAB(DdVDp.GetArray(), DSDp.GetPtr(), ddq, 6, bjDOF, bjDOF, 1);
	}

	//DdVDp += convert_to_se3(S * pBaseJoint->get_DddqDp());
	pBaseJoint->get_DddqDp(DddqDp);
	matAdd_multAB(DdVDp.GetArray(), S.GetPtr(), DddqDp, 6, bjDOF, bjDOF, 1);

	DdVDp += DetaDp;
}

void GBody::update_DFDp()
{
	std::list<GBody *>::iterator iter_pbody_child;

	if ( bDpAlien ) {
		DFDp = I * DdVDp;
		if ( bnzDFeDp ) { DFDp -= DFeDp; }
		for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
			DFDp += (*iter_pbody_child)->getTransformed_DFDp();
		}
		return;
	}

	dse3 tmp = I * DVDp;

	DFDp = I * DdVDp - dad(DVDp, I*V);

	if ( bnzDIDp ) {
		DFDp += DIDp * dV;
		tmp += DIDp * V;
	}

	DFDp -= dad(V, tmp);

	if ( bnzDFeDp ) {
		DFDp -= DFeDp;
	}

	for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
		DFDp += (*iter_pbody_child)->getTransformed_DFDp();
	}
}

void GBody::update_DFDp_fs()
{
	if ( bDpAlien ) {
		DFDp = aI * DdVDp;
		DFDp += DaBDp;	
		return;
	}

	DFDp = DaIDp * dV;
	DFDp += aI * DdVDp;
	DFDp += DaBDp;	
}

void GBody::update_DaIDp()
{
	if ( bDpAlien ) {
		DaIDp.SetZero();
		return;
	}

	std::list<GBody *>::iterator iter_pbody_child;

	DaIDp.SetZero();

	if ( bnzDIDp ) {
		DaIDp += DIDp;
	}

	for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
		DaIDp += (*iter_pbody_child)->getTransformed_DaIDp();
	}
}

void GBody::update_DaBDp()
{
	std::list<GBody *>::iterator iter_pbody_child;

	if ( bDpAlien ) {
		DaBDp.SetZero();
		if ( bnzDFeDp ) { 
			DaBDp -= DFeDp; 
		}
		for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
			DaBDp += (*iter_pbody_child)->getTransformed_DaBDp();
		}
		return;
	}

	dse3 tmp = I*DVDp;

	DaBDp = -dad(DVDp, I*V);

	if ( bnzDIDp ) {
		tmp += DIDp * V;
	}

	DaBDp -= dad(V, tmp);

	if ( bnzDFeDp ) {
		DaBDp -= DFeDp;
	}

	for (iter_pbody_child = pChildBodies.begin(); iter_pbody_child != pChildBodies.end(); iter_pbody_child++) {
		DaBDp += (*iter_pbody_child)->getTransformed_DaBDp();
	}
}

void GBody::update_DPsiDp()
{
	// DPsiDp = - Psi * ( ~S * DaIDp * S + (~DSDp * aI_S) + ~(~DSDp * aI_S) ) * Psi

	gReal tmp1[36], tmp2[36], tmp3[36], tmp4[36], tmp5[36];

	if ( bDpAlien ) {
		DPsiDp.SetZero(0,0);
		return;
	}

	if ( pBaseJoint->isPrescribed() ) {
		;
	} else {
		if ( bjDOF > 0 ) {
			// tmp2 = ~S * DaIDp * S + (~DSDp * aI_S) + ~(~DSDp * aI_S)
			set_Mult_AInertia_se3(tmp1, DaIDp, S.GetPtr(), bjDOF); // tmp1 = DaIDp * S
			matSet_multAtB(tmp2, S.GetPtr(), tmp1, 6, bjDOF, 6, bjDOF);	// tmp2 = ~S * (DaIDp * S)
			if ( bnzDSDp ) {
				matSet_multAtB(tmp3, DSDp.GetPtr(), aI_S.GetPtr(), 6, bjDOF, 6, bjDOF);	// tmp3 = ~DSDp * aI_S
				matSet_transpose(tmp4, tmp3, bjDOF, bjDOF);	// tmp4 = ~(~DSDp * aI_S)
				matAdd(tmp2, tmp3, bjDOF, bjDOF);
				matAdd(tmp2, tmp4, bjDOF, bjDOF);
			}
			matSet_multAB(tmp5, tmp2, Psi.GetPtr(), bjDOF, bjDOF, bjDOF, bjDOF);	// tmp5 = tmp2 * Psi
			matSet_multAB(DPsiDp.GetPtr(), Psi.GetPtr(), tmp5, bjDOF, bjDOF, bjDOF, bjDOF); // DPsiDp = Psi * tmp5
			matMult(DPsiDp.GetPtr(), -1.0, bjDOF, bjDOF); // DPsiDp *= -1
		}
		//if ( bjDOF <= 0 ) {
		//	DPsiDp.SetZero(0,0);
		//} else {
		//	RMatrix tmp = ~S * convert_to_RMatrix(DaIDp) * S;
		//	if ( bnzDSDp ) {
		//		RMatrix tmp2 = ~DSDp * convert_to_RMatrix(aI) * S;
		//		tmp += tmp2 + ~tmp2;
		//	}
		//	DPsiDp = -Psi*tmp*Psi;
		//}
	}
}

void GBody::update_DPiDp()
{
	// DPiDp = DaIDp 
	//		 - aI_S * DPsiDp * ~aI_S 
	//		 - ( (DaIDp * S * Psi * ~aI_S) + ~(DaIDp * S * Psi * ~aI_S) )
	//		 - ( (aI * DSDp * Psi * ~aI_S) + ~(aI * DSDp * Psi * ~aI_S) )

	gReal tmp1[36], tmp2[36], tmp3[36], tmp4[36];

	if ( bDpAlien ) {
		DPiDp.SetZero();
		return;
	}

	if ( pBaseJoint->isPrescribed() ) {
		DPiDp = DaIDp;
	} else {
		if ( bjDOF <= 0 ) {
			DPiDp = DaIDp;
		} else {
			// DPiDp = DaIDp
			DPiDp = DaIDp;
			// DPiDp -= aI_S * DPsiDp * ~aI_S
			matSet_multABt(tmp1, DPsiDp.GetPtr(), aI_S.GetPtr(), bjDOF, bjDOF, 6, bjDOF); // tmp1 = DPsiDp * ~aI_S
			matSet_multAB(tmp2, aI_S.GetPtr(), tmp1, 6, bjDOF, bjDOF, 6); // tmp2 = aI_S * (DPsiDp * ~aI_S)
			DPiDp -= tmp2;
			// DPiDp -= (DaIDp * S * Psi * ~aI_S) + ~(DaIDp * S * Psi * ~aI_S)
			matSet_multABt(tmp1, Psi.GetPtr(), aI_S.GetPtr(), bjDOF, bjDOF, 6, bjDOF); // tmp1 = Psi * ~aI_S (this will be used later again in the if statement)
			set_Mult_AInertia_se3(tmp2, DaIDp, S.GetPtr(), bjDOF);	// tmp2 = DaIDp * S
			matSet_multAB(tmp3, tmp2, tmp1, 6, bjDOF, bjDOF, 6); // tmp3 = (DaIDp * S) * (Psi * ~aI_S)
			matSet_transpose(tmp4, tmp3, 6, 6);	// tmp4 = ~tmp3
			DPiDp -= tmp3; DPiDp -= tmp4; // DPiDp -= (tmp3 + ~tmp3)
			if ( bnzDSDp ) {
				// DPiDp -= (aI * DSDp * Psi * ~aI_S) + ~(aI * DSDp * Psi * ~aI_S)
				set_Mult_AInertia_se3(tmp2, aI, DSDp.GetPtr(), bjDOF); // tmp2 = aI * DSDp
				matSet_multAB(tmp3, tmp2, tmp1, 6, bjDOF, bjDOF, 6); // tmp3 = (aI * DSDp) * (Psi * ~aI_S)
				matSet_transpose(tmp4, tmp3, 6, 6);	// tmp4 = ~tmp3
				DPiDp -= tmp3; DPiDp -= tmp4; // DPiDp -= (tmp3 + ~tmp3)
			}

			//RMatrix tmp = ( convert_to_RMatrix(DaIDp) * S ) * Psi * ~aI_S;
			//DPiDp = DaIDp;
			//DPiDp -= AInertia((aI_S * DPsiDp * ~aI_S).GetPtr());
			//DPiDp -= AInertia((tmp + ~tmp).GetPtr());
			//if ( bnzDSDp ) {
			//	RMatrix tmp2 = ( convert_to_RMatrix(aI) * DSDp ) * Psi * ~aI_S;
			//	DPiDp -= AInertia((tmp2 + ~tmp2).GetPtr());
			//}
		}
	}
}

void GBody::update_DbetaDp()
{
	gReal vtmp1[6], vtmp2[6], vtmp3[6], vtmp4[6], vtmp5[6];
	gReal stmp1[6], stmp2[6];

	if ( bDpAlien ) {
		if ( pBaseJoint->isPrescribed() ) {
			DbetaDp = DaBDp;
			if ( bjDOF > 0 ) { 
				//DbetaDp += aI * convert_to_se3(S * pBaseJoint->get_DddqDp());
				pBaseJoint->get_DddqDp(vtmp1); // vtmp1 = DddqDp
				matSet_multAB(stmp1, S.GetPtr(), vtmp1, 6, bjDOF, bjDOF, 1); // stmp1 = S * DddqDp
				add_Mult_AInertia_se3(DbetaDp, aI, stmp1); // DbetaDp += aI * (S * DddqDp)
			}
		} else {
			DbetaDp = DaBDp;
			if ( bjDOF > 0 ) {
				//RMatrix S_Psi = S * Psi;
				//RMatrix tmp3 = pBaseJoint->get_DtauDp() - ~S * convert_to_RMatrix(DaBDp);
				//DbetaDp += aI * convert_to_se3(S_Psi * tmp3);
				pBaseJoint->get_DtauDp(vtmp1); // vtmp1 = DtauDp
				for (int i=0; i<bjDOF; i++) { vtmp1[i] -= DaBDp.InnerProductWith(&S[6*i]); } // vtmp1 -= ~S * DaBDp
				matSet_multAB(vtmp2, Psi.GetPtr(), vtmp1, bjDOF, bjDOF, bjDOF, 1); // vtmp2 = Psi * (DtauDp - ~S * DaBDp)
				matSet_multAB(stmp1, S.GetPtr(), vtmp2, 6, bjDOF, bjDOF, 1); // stmp1 = S * Psi * (DtauDp - ~S * DaBDp)
				add_Mult_AInertia_se3(DbetaDp, aI, stmp1); // DbetaDp += aI * S * Psi * (DtauDp - ~S * DaBDp)
			}
		}
		return;
	}

	if ( pBaseJoint->isPrescribed() ) {
		if ( bjDOF <= 0 ) {
			//DbetaDp = DaBDp;
			//DbetaDp += DaIDp * eta;
			//DbetaDp += aI * DetaDp;
			DbetaDp = DaBDp;
			add_Mult_AInertia_se3(DbetaDp, DaIDp, eta);
			add_Mult_AInertia_se3(DbetaDp, aI, DetaDp);
		} else {
			//se3 tmp(DetaDp);
			//tmp += convert_to_se3(S * pBaseJoint->get_DddqDp());
			//if ( bnzDSDp ) {
			//	tmp += convert_to_se3(DSDp * pBaseJoint->get_ddq());
			//}
			//DbetaDp = DaBDp;
			//DbetaDp += DaIDp * (eta + Sddq);
			//DbetaDp += aI * tmp;
			pBaseJoint->get_DddqDp(vtmp1); // vtmp1 = DddqDp
			matSet(stmp1, DetaDp.GetArray(), 6); // stmp1 = DetaDp
			matAdd_multAB(stmp1, S.GetPtr(), vtmp1, 6, bjDOF, bjDOF, 1); // stmp1 += S * DddqDp
			if ( bnzDSDp ) {
				pBaseJoint->get_ddq(vtmp2); // vtmp2 = ddq
				matAdd_multAB(stmp1, DSDp.GetPtr(), vtmp2, 6, bjDOF, bjDOF, 1); // stmp1 += DSDp * ddq
			}
			DbetaDp = DaBDp;
			DbetaDp += DaIDp * (eta + Sddq);
			add_Mult_AInertia_se3(DbetaDp, aI, stmp1); // DbetaDp += aI * stmp1
		}
	} else {
		if ( bjDOF <= 0 ) {
			//DbetaDp = DaBDp;
			//DbetaDp += DaIDp * eta;
			//DbetaDp += aI * DetaDp;
			DbetaDp = DaBDp;
			add_Mult_AInertia_se3(DbetaDp, DaIDp, eta);
			add_Mult_AInertia_se3(DbetaDp, aI, DetaDp);
		} else {
			//dse3 aI_eta_aB = aI * eta + aB;
			//RMatrix S_Psi = S * Psi;
			//RMatrix S_DPsiDp = S * DPsiDp;
			//RMatrix tmp = pBaseJoint->get_tau() - ~S * convert_to_RMatrix(aI_eta_aB);
			//RMatrix tmp2 = S_DPsiDp;
			//RMatrix tmp3 = pBaseJoint->get_DtauDp() - ~S * convert_to_RMatrix(DaIDp * eta + aI * DetaDp + DaBDp);
			//if ( bnzDSDp ) {
			//	tmp2 += DSDp * Psi;
			//	tmp3 -= ~DSDp * convert_to_RMatrix(aI_eta_aB);
			//}
			//DbetaDp = DaBDp;
			//DbetaDp += DaIDp * (eta + convert_to_se3(S_Psi * tmp));
			//DbetaDp += aI * (DetaDp + convert_to_se3(tmp2 * tmp + S_Psi * tmp3));

			pBaseJoint->get_tau(vtmp1); // vtmp1 = tau
			pBaseJoint->get_DtauDp(vtmp2); // vtmp2 = DtauDp
			dse3 dse3tmp1 = aI * eta + aB;
			dse3 dse3tmp2 = DaIDp * eta + aI * DetaDp + DaBDp;
			for (int i=0; i<bjDOF; i++) { 
				vtmp1[i] -= dse3tmp1.InnerProductWith(&S[6*i]); // vtmp1 = tau - ~S * (aI * eta + aB)
				vtmp2[i] -= dse3tmp2.InnerProductWith(&S[6*i]); // vtmp2 = DtauDp - ~S * (DaIDp * eta + aI * DetaDp + DaBDp)
			} 
			matSet_multAB(vtmp3, Psi.GetPtr(), vtmp1, bjDOF, bjDOF, bjDOF, 1); // vtmp3 = Psi * ( tau - ~S * (aI * eta + aB) )
			matSet_multAB(vtmp4, DPsiDp.GetPtr(), vtmp1, bjDOF, bjDOF, bjDOF, 1); // vtmp4 = DPsiDp * ( tau - ~S * (aI * eta + aB) )
			matSet(stmp1, eta.GetArray(), 6); // stmp1 = eta
			matSet(stmp2, DetaDp.GetArray(), 6); // stmp2 = DetaDp
			matAdd_multAB(stmp1, S.GetPtr(), vtmp3, 6, bjDOF, bjDOF, 1); // stmp1 = eta + S * Psi * ( tau - ~S * (aI * eta + aB) )
			matAdd_multAB(stmp2, S.GetPtr(), vtmp4, 6, bjDOF, bjDOF, 1); // stmp2 = DetaDp + S * DPsiDp * ( tau - ~S * (aI * eta + aB) )
			if ( bnzDSDp ) {
				matAdd_multAB(stmp2, DSDp.GetPtr(), vtmp3, 6, bjDOF, bjDOF, 1); // stmp2 += DSDp * Psi * ( tau - ~S * (aI * eta + aB) )
				for (int i=0; i<bjDOF; i++) { 
					vtmp2[i] -= dse3tmp1.InnerProductWith(&DSDp[6*i]); // vtmp2 -= ~DSDp * (aI * eta + aB)
				}
			}
			matSet_multAB(vtmp5, Psi.GetPtr(), vtmp2, bjDOF, bjDOF, bjDOF, 1); // vtmp5 = Psi * ( DtauDp - ~S * (DaIDp * eta + aI * DetaDp + DaBDp) - ~DSDp * (aI * eta + aB) )
			matAdd_multAB(stmp2, S.GetPtr(), vtmp5, 6, bjDOF, bjDOF, 1); // stmp2 += S * Psi * ( DtauDp - ~S * (DaIDp * eta + aI * DetaDp + DaBDp) - ~DSDp * (aI * eta + aB) )
			DbetaDp = DaBDp;
			add_Mult_AInertia_se3(DbetaDp, DaIDp, stmp1);
			add_Mult_AInertia_se3(DbetaDp, aI, stmp2);
		}
	}
}

void GBody::update_DtauDp()
{
	gReal vtmp1[6];

	int i;
	std::list<GCoordinate *>::iterator iter_pcoord;

	if ( bjDOF <= 0 ) return;

	//RMatrix _DtauDp = ~S * convert_to_RMatrix(DFDp);
	//if ( bnzDSDp ) {
	//	_DtauDp += ~DSDp * convert_to_RMatrix(F);
	//}
	//for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
	//	(*iter_pcoord)->DtauDp = _DtauDp[i];
	//}

	for (int i=0; i<bjDOF; i++) {
		vtmp1[i] = DFDp.InnerProductWith(&S[6*i]); // vtmp = ~S * DFDp
	}
	if ( bnzDSDp ) {
		for (int i=0; i<bjDOF; i++) {
			vtmp1[i] += F.InnerProductWith(&DSDp[6*i]); // vtmp += ~DSDp * F
		}
	}
	for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
		(*iter_pcoord)->DtauDp = vtmp1[i];
	}
}

void GBody::update_DddqDp()
{
	int i;
	std::list<GCoordinate *>::iterator iter_pcoord;
	gReal vtmp1[6], vtmp2[6], vtmp3[6];

	if ( bjDOF <= 0 ) return;

	if ( bDpAlien ) {
		//RMatrix tmp1 = pBaseJoint->get_DtauDp() - ~S * convert_to_RMatrix(DaBDp);
		//se3 tmp3; tmp3.set_Ad(invT, pParentBody->DdVDp());
		//RMatrix _DddqDp = Psi * ( tmp1 - ~S * convert_to_RMatrix(aI * tmp3) );
		//for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
		//	(*iter_pcoord)->DddqDp = _DddqDp[i];
		//}

		pBaseJoint->get_DtauDp(vtmp1); // vtmp1 = DtauDp
		for (i=0; i<bjDOF; i++) {
			vtmp1[i] -= DaBDp.InnerProductWith(&S[6*i]); // vtmp1 = DtauDp - ~S * DaBDp
		}
		se3 Ad_DdVDp; Ad_DdVDp.set_Ad(invT, pParentBody->DdVDp);
		dse3 aI_Ad_DdVDp = aI * Ad_DdVDp;
		for (i=0; i<bjDOF; i++) {
			vtmp1[i] -= aI_Ad_DdVDp.InnerProductWith(&S[6*i]); // vtmp1 = DtauDp - ~S * DaBDp - ~S * aI * Ad(DdVDp)
		}
		matSet_multAB(vtmp2, Psi.GetPtr(), vtmp1, bjDOF, bjDOF, bjDOF, 1); // vtmp2 = Psi * vtmp1
		for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
			(*iter_pcoord)->DddqDp = vtmp2[i];
		}

		return;
	}

	//se3 Ad_dV; Ad_dV.set_Ad(invT, pParentBody->dV);
	//se3 Ad_dV_eta = Ad_dV + eta;

	//// tmp1 = DtauDp - ~DSDp * aB - ~S * DaBDp;
	//// tmp2 = (~DSDp * aI + ~S * DaIDp) * (Ad(invT, pParentBody->dV) + eta);
	//// tmp3 = Ad(invT, pParentBody->DdVDp()) - ad(DhDp, Ad(invT, pParentBody->dV)) + DetaDp;
	//RMatrix tmp1 = pBaseJoint->get_DtauDp() - ~S * convert_to_RMatrix(DaBDp);
	//RMatrix tmp2 = ~S * convert_to_RMatrix(DaIDp * Ad_dV_eta);
	//se3 tmp3; tmp3.set_Ad(invT, pParentBody->DdVDp()); tmp3 += DetaDp;
	//if ( bnzDSDp ) {
	//	tmp1 -= ~DSDp * convert_to_RMatrix(aB);
	//	tmp2 += ~DSDp * convert_to_RMatrix(aI * Ad_dV_eta);
	//}
	//if ( bnzDhDp ) {
	//	tmp3 -= ad(DhDp, Ad_dV);
	//}

	//// DddqDp
	//RMatrix _DddqDp = DPsiDp * ( pBaseJoint->get_tau() - ~S * convert_to_RMatrix(aI * Ad_dV_eta + aB) );
	//_DddqDp += Psi * ( tmp1 - tmp2 - ~S * convert_to_RMatrix(aI * tmp3) );

	//for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
	//	(*iter_pcoord)->DddqDp = _DddqDp[i];
	//}

	se3 Ad_dV; Ad_dV.set_Ad(invT, pParentBody->dV); // Ad_dV = Ad(dV)
	se3 Ad_dV_eta(Ad_dV); Ad_dV_eta += eta; // Ad_dV_eta = Ad(dV) + eta
	pBaseJoint->get_tau(vtmp1); // vtmp1 = tau
	dse3 dse3tmp1 = aI * Ad_dV_eta + aB;
	for (i=0; i<bjDOF; i++) {
		vtmp1[i] -= dse3tmp1.InnerProductWith(&S[6*i]); // vtmp1 -= ~S * ( aI * (Ad(dV) + eta ) + aB )
	}
	pBaseJoint->get_DtauDp(vtmp2); // vtmp2 = DtauDp
	if ( bnzDhDp ) {
		// vtmp2 -= ~S * ( DaIDp * (Ad(dV)+eta) + aI * (Ad(DdVDp) + DetaDp - ad(DhDp, Ad(dV))) + DaBDp )
		dse3 dse3tmp2 = DaIDp * Ad_dV_eta + aI * ( Ad(invT, pParentBody->DdVDp) + DetaDp - ad(DhDp, Ad_dV) ) + DaBDp;
		for (i=0; i<bjDOF; i++) {
			vtmp2[i] -= dse3tmp2.InnerProductWith(&S[6*i]); 
		}
	} else {
		// vtmp2 -= ~S * ( DaIDp * (Ad(dV)+eta) + aI * (Ad(DdVDp) + DetaDp) + DaBDp )
		dse3 dse3tmp2 = DaIDp * Ad_dV_eta + aI * ( Ad(invT, pParentBody->DdVDp) + DetaDp ) + DaBDp;
		for (i=0; i<bjDOF; i++) {
			vtmp2[i] -= dse3tmp2.InnerProductWith(&S[6*i]); 
		}
	}
	if ( bnzDSDp ) {
		// vtmp2 -= ~DSDp * ( aI * (Ad(dV)+eta) + aB )
		dse3 dse3tmp3 = aI * Ad_dV_eta + aB;
		for (i=0; i<bjDOF; i++) {
			vtmp2[i] -= dse3tmp3.InnerProductWith(&DSDp[6*i]); 
		}
	}
	// vtmp3 = DPsiDp * vtmp1 + Psi * vtmp2
	matSet_multAB(vtmp3, DPsiDp.GetPtr(), vtmp1, bjDOF, bjDOF, bjDOF, 1);
	matAdd_multAB(vtmp3, Psi.GetPtr(), vtmp2, bjDOF, bjDOF, bjDOF, 1);
	for (iter_pcoord = pBaseJoint->pCoordinates.begin(), i=0; iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
		(*iter_pcoord)->DddqDp = vtmp3[i];
	}
}

se3 GBody::get_S(GCoordinate *pCoordinate_)
{
	int i;
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (i=0, iter_pcoord = pBaseJoint->pCoordinates.begin(); iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
		if ( *iter_pcoord == pCoordinate_ ) { return se3(S(0,i), S(1,i), S(2,i), S(3,i), S(4,i), S(5,i)); }
	}

	return se3(0,0,0,0,0,0);
}

se3 GBody::get_dS(GCoordinate *pCoordinate_)
{
	int i;
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (i=0, iter_pcoord = pBaseJoint->pCoordinates.begin(); iter_pcoord != pBaseJoint->pCoordinates.end(); iter_pcoord++, i++) {
		if ( *iter_pcoord == pCoordinate_ ) { return se3(dS(0,i), dS(1,i), dS(2,i), dS(3,i), dS(4,i), dS(5,i)); }
	}

	return se3(0,0,0,0,0,0);
}

RMatrix GBody::get_DSDq(GCoordinate *pCoordinate_)
{
	if ( isIncluded(pCoordinate_) ) {
		return Ad(pBaseJoint->T_right, pBaseJoint->get_DSDq(pCoordinate_));
	} else {
		return Zeros(6, bjDOF);
	}
}

RMatrix GBody::get_DdSDq(GCoordinate *pCoordinate_)
{
	if ( isIncluded(pCoordinate_) ) {
		return Ad(pBaseJoint->T_right, pBaseJoint->get_DdSDq(pCoordinate_));
	} else {
		return Zeros(6, bjDOF);
	}
}

RMatrix GBody::get_DdSDdq(GCoordinate *pCoordinate_)
{
	return get_DSDq(pCoordinate_);
}

bool GBody::isIncluded(GCoordinate *pCoordinate_)
{
	if ( find(pBaseJoint->pCoordinates.begin(), pBaseJoint->pCoordinates.end(), pCoordinate_) != pBaseJoint->pCoordinates.end() ) 
		return true;
	else
		return false;
}

std::string GBody::getInfoStr()
{
	std::stringstream sstr;
	std::list<GJoint *>::iterator iter_pjoint;
	std::list<GBody *>::iterator iter_pbody;

	sstr << GElement::getInfoStr();
	sstr << "GBody:: " << std::endl;
	sstr << "number of joints attached = " << int(pJoints.size()) << std::endl;
	sstr << "joints attached = ("; for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) { sstr << (*iter_pjoint)->getName() << ", "; } sstr << ")" << std::endl;
	sstr << "base joint = " << pBaseJoint->getName() << std::endl;
	sstr << "parent body = " << pParentBody->getName() << std::endl;
	sstr << "child bodies = ("; for (iter_pbody = pChildBodies.begin(); iter_pbody != pChildBodies.end(); iter_pbody++) { sstr << (*iter_pbody)->getName() << ", "; } sstr << ")" << std::endl;
	sstr << "forward joint chain = ("; for (iter_pjoint = fwdJointChain.pJoints.begin(); iter_pjoint != fwdJointChain.pJoints.end(); iter_pjoint++) { sstr << (*iter_pjoint)->getName() << ", "; } sstr << ")" << std::endl;
	sstr << "mass = " << getMass() << std::endl;
	sstr << "center of mass (local)  = " << getPositionCOM();
	sstr << "center of mass (global) = " << getPositionCOMGlobal();
	sstr << "generalized momentum (local)  = " << getMomentum();
	sstr << "generalized momentum (global) = " << getMomentumGlobal();
	sstr << "T_global = " << T_global;
	sstr << "I = " << I;
	sstr << "Fe = " << Fe;
	sstr << "T = " << T;
	sstr << "invT = " << invT;
	sstr << "V = " << V;
	sstr << "dV = " << dV;
	sstr << "F = " << F;
	sstr << "aI = " << aI;
	sstr << "aB = " << aB;
	sstr << "eta = " << eta;
	sstr << "Psi = " << Psi;
	sstr << "Pi = " << Pi;
	sstr << "beta = " << beta;

	return sstr.str();
}
