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
#include "gspringdamper.h"

//=============================================================
//                 GSpringDamperBody
//=============================================================

GSpringDamperBody::GSpringDamperBody()
{ 
	forceType = GFORCE_SPRINGDAMPER_BODY;
	pLeftBody = pRightBody = NULL;
	T_left.SetIdentity();
	T_right.SetIdentity();
	inv_T_left.SetIdentity();
	inv_T_right.SetIdentity();
	K[0] = K[1] = K[2] = K[3] = K[4] = K[5] = 0;
	C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0;
}

GSpringDamperBody::~GSpringDamperBody() 
{ 
}

void GSpringDamperBody::clear()
{
	pLeftBody = pRightBody = NULL;
	T_left.SetIdentity();
	T_right.SetIdentity();
	inv_T_left.SetIdentity();
	inv_T_right.SetIdentity();
	K[0] = K[1] = K[2] = K[3] = K[4] = K[5] = 0;
	C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0;
}

bool GSpringDamperBody::connectBodies(GBody *pLeftBody_, GBody *pRightBody_)
{
	if ( pLeftBody_ == NULL || pRightBody_ == NULL ) return false;

	pLeftBody = pLeftBody_;
	pRightBody = pRightBody_;

	pLeftBody->pForces.push_back(this);
	pRightBody->pForces.push_back(this);

	return true;
}

bool GSpringDamperBody::disconnectBodies()
{
	if ( pLeftBody == NULL || pRightBody == NULL ) return false;

	pLeftBody->pForces.remove(this);
	pRightBody->pForces.remove(this);

	pLeftBody = NULL;
	pRightBody = NULL;

	return true;
}

void GSpringDamperBody::setPosition(const Vec3 &pL_, const Vec3 &pR_)
{
	T_left.SetPosition(pL_);
	T_right.SetPosition(pR_);
	inv_T_left = Inv(T_left);
	inv_T_right = Inv(T_right);
}

void GSpringDamperBody::setOrientation(const SO3 &RL_, const SO3 &RR_)
{
	T_left.SetRotation(RL_);
	T_right.SetRotation(RR_);
	inv_T_left = Inv(T_left);
	inv_T_right = Inv(T_right);
}

void GSpringDamperBody::setPositionAndOrientation(const SE3 &TL_, const SE3 &TR_)
{
	T_left = TL_;
	T_right = TR_;
	inv_T_left = Inv(T_left);
	inv_T_right = Inv(T_right);
}

void GSpringDamperBody::setSpring(const Vec3 &K_rotation_, const Vec3 &K_position_)
{
	K[0] = K_rotation_[0];	K[1] = K_rotation_[1];	K[2] = K_rotation_[2];
	K[3] = K_position_[0];	K[4] = K_position_[1];	K[5] = K_position_[2];
}

void GSpringDamperBody::setDamper(const Vec3 &C_rotation_, const Vec3 &C_position_)
{
	C[0] = C_rotation_[0];	C[1] = C_rotation_[1];	C[2] = C_rotation_[2];
	C[3] = C_position_[0];	C[4] = C_position_[1];	C[5] = C_position_[2];
}

bool GSpringDamperBody::applyForce(bool badd_)
{
	if ( pLeftBody == NULL || pRightBody == NULL ) return false;

	SE3 T = inv_T_left * Inv(pLeftBody->T_global) * pRightBody->T_global * T_right;

	//se3 x(Log(T.GetRotation()), T.GetPosition());	// relative location of the right body from the left body
	se3 x(Log(T));	// relative location of the right body from the left body

	se3 V_left = Ad(inv_T_left, pLeftBody->V) - Ad(T * inv_T_right, pRightBody->V);		// relative velocity of the left body

	dse3 F_left;	// force to be acting on the left body
	for (int i=0; i<6; i++) {
		F_left[i] = K[i] * x[i] - C[i] * V_left[i]; 
	}

	if ( badd_ ) {
		pLeftBody->Fe += dAd(inv_T_left, F_left);
		pRightBody->Fe += dAd(T * inv_T_right, -F_left);
	} else {
		pLeftBody->Fe -= dAd(inv_T_left, F_left);
		pRightBody->Fe -= dAd(T * inv_T_right, -F_left);
	}

	return true;
}

GSpringDamperJoint::GSpringDamperJoint()
{
	forceType = GFORCE_SPRINGDAMPER_JOINT;
	pJoint = NULL;
	q_neutral = NULL;
	K = C = NULL;
}

GSpringDamperJoint::~GSpringDamperJoint()
{
	delete [] q_neutral;
	delete [] K;
	delete [] C;
}

void GSpringDamperJoint::clear()
{
	delete [] q_neutral;
	delete [] K;
	delete [] C;

	pJoint = NULL;
	q_neutral = NULL;
	K = C = NULL;
}

bool GSpringDamperJoint::connectJoint(GJoint *pjoint_)
{
	if ( pjoint_ == NULL ) return false;

	int n = pjoint_->getDOF();

	pJoint = pjoint_;
	pJoint->pForces.push_back(this);

	q_neutral = new gReal[n];
	K = new gReal[n];
	C = new gReal[n];

	for (int i=0; i<n; i++) {
		q_neutral[i] = K[i] = C[i] = 0;
	}

	return true;
}

bool GSpringDamperJoint::disconnectJoint()
{
	if ( pJoint == NULL ) return false;
	pJoint->pForces.remove(this);
	pJoint = NULL;
	return true;
}

bool GSpringDamperJoint::setNeutralPosition(gReal *q_neutral_)
{
	if ( pJoint == NULL ) return false;

	if ( q_neutral_ == NULL ) {
		RMatrix q = pJoint->get_q();
		for (int i=0; i<pJoint->getDOF(); i++) {
			q_neutral[i] = q[i];
		}
	} else {
		for (int i=0; i<pJoint->getDOF(); i++) {
			q_neutral[i] = q_neutral_[i];
		}
	}

	return true;
}

bool GSpringDamperJoint::setNeutralPosition(gReal q_neutral_)
{
	if ( pJoint == NULL ) return false;

	for (int i=0; i<pJoint->getDOF(); i++) {
		q_neutral[i] = q_neutral_;
	}

	return true;
}

bool GSpringDamperJoint::setNeutralPosition(RMatrix q_neutral_)
{
	if ( q_neutral_.RowSize() * q_neutral_.ColSize() != pJoint->getDOF() ) return false;
	return setNeutralPosition(q_neutral_.GetPtr());
}

bool GSpringDamperJoint::setSpring(gReal *K_)
{
	if ( pJoint == NULL ) return false;
	if ( K_ == NULL ) return false;

	for (int i=0; i<pJoint->getDOF(); i++) {
		K[i] = K_[i];
	}

	return true;
}

bool GSpringDamperJoint::setSpring(gReal K_)
{
	if ( pJoint == NULL ) return false;

	for (int i=0; i<pJoint->getDOF(); i++) {
		K[i] = K_;
	}

	return true;
}

bool GSpringDamperJoint::setSpring(RMatrix K_)
{
	if ( K_.RowSize() * K_.ColSize() != pJoint->getDOF() ) return false;
	return setSpring(K_.GetPtr());
}

bool GSpringDamperJoint::setDamper(gReal *C_)
{
	if ( pJoint == NULL ) return false;
	if ( C_ == NULL ) return false;

	for (int i=0; i<pJoint->getDOF(); i++) {
		C[i] = C_[i];
	}
	
	return true;
}

bool GSpringDamperJoint::setDamper(gReal C_)
{
	if ( pJoint == NULL ) return false;

	for (int i=0; i<pJoint->getDOF(); i++) {
		C[i] = C_;
	}

	return true;
}

bool GSpringDamperJoint::setDamper(RMatrix C_)
{
	if ( C_.RowSize() * C_.ColSize() != pJoint->getDOF() ) return false;
	return setDamper(C_.GetPtr());
}

bool GSpringDamperJoint::applyForce(bool badd_)
{
	if ( pJoint == NULL ) return false;

	int i;
	std::list<GCoordinate*>::iterator iter_pcoord;

	if ( badd_ ) {
		for (i=0, iter_pcoord = pJoint->pCoordinates.begin(); iter_pcoord != pJoint->pCoordinates.end(); iter_pcoord++, i++) {
			(*iter_pcoord)->tau -= K[i] * ( (*iter_pcoord)->q - q_neutral[i] ) + C[i] * (*iter_pcoord)->dq;
		}
	} else {
		for (i=0, iter_pcoord = pJoint->pCoordinates.begin(); iter_pcoord != pJoint->pCoordinates.end(); iter_pcoord++, i++) {
			(*iter_pcoord)->tau += K[i] * ( (*iter_pcoord)->q - q_neutral[i] ) + C[i] * (*iter_pcoord)->dq;
		}
	}

	return true;
}