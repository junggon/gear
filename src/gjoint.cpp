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
#include <string>
#include <sstream>
#include "gjoint.h"
#include "gelement.h"
#include "gbody.h"
#include "gcoordinate.h"
#include "liegroup.h"

//=============================================================
//                 GJoint
//=============================================================
GJoint::GJoint()
{
	jointType = GJOINT_NULL;
	pLeftBody = pRightBody = NULL;
	T_left.SetIdentity();
	T_right.SetIdentity();
	inv_T_left.SetIdentity();
	inv_T_right.SetIdentity();
	T_global.SetIdentity();
	T.SetIdentity();
	inv_T.SetIdentity();
	Sdq.SetZero();
	dSdq.SetZero();
	Sddq.SetZero();
	DSdqDt.SetZero();
	bReversed = false;
	bCut = false;
	setPrescribed(false);
}

bool GJoint::connectBodies(GBody *pLeftBody_, GBody *pRightBody_)
{
	if ( pLeftBody_ == NULL || pRightBody_ == NULL ) return false;

	pLeftBody = pLeftBody_;
	pRightBody = pRightBody_;
	bReversed = false;
	bCut = false;
	pLeftBody->pJoints.push_back(this);
	pRightBody->pJoints.push_back(this);
	return true;
}

void GJoint::disconnectBodies()
{
	if ( pLeftBody != NULL ) {
		pLeftBody->pJoints.remove(this);
		pLeftBody = NULL;
	}
	if ( pRightBody != NULL ) {
		pRightBody->pJoints.remove(this);
		pRightBody = NULL;
	}
	bReversed = false;
	bCut = false;
}

void GJoint::setPosition(const Vec3 &pL_, const Vec3 &pR_)
{
	T_left.SetPosition(pL_);
	T_right.SetPosition(pR_);
	inv_T_left.SetInvOf(T_left);
	inv_T_right.SetInvOf(T_right);
}

void GJoint::setOrientation(const SO3 &RL_, const SO3 &RR_)
{
	T_left.SetRotation(RL_);
	T_right.SetRotation(RR_);
	inv_T_left.SetInvOf(T_left);
	inv_T_right.SetInvOf(T_right);
}

void GJoint::setPositionAndOrientation(const SE3 &TL_, const SE3 &TR_)
{
	T_left = TL_;
	T_right = TR_;
	inv_T_left.SetInvOf(T_left);
	inv_T_right.SetInvOf(T_right);
}

// set pCoordinates[]->bPrescribed
void GJoint::setPrescribed(bool b_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;

	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->bPrescribed = b_;
	}

	bPrescribed = b_;
}

bool GJoint::reverse()
{
	GBody *pbody_tmp;
	SE3 T_tmp;

	// swap pLeftBody and pRightBody
	pbody_tmp = pLeftBody;
	pLeftBody = pRightBody;
	pRightBody = pbody_tmp;

	// swap T_left and T_right
	T_tmp = T_left;
	T_left = T_right;
	T_right = T_tmp;

	// swap inv_T_left and inv_T_right
	T_tmp = inv_T_left;
	inv_T_left = inv_T_right;
	inv_T_right = T_tmp;

	// set flag
	bReversed = (bReversed == true ? false : true );

	return true;
}

std::string GJoint::getInfoStr()
{
	std::stringstream sstr;
	std::list<GCoordinate *>::iterator iter_pcoord;
	std::string body_name_1, body_name_2;
	if ( pLeftBody == NULL ) {
		body_name_1 = "NULL";
	} else {
		body_name_1 = pLeftBody->getName();
	}
	if ( pRightBody == NULL ) {
		body_name_2 = "NULL";
	} else {
		body_name_2 = pRightBody->getName();
	}

	sstr << GElement::getInfoStr();
	sstr << "GJoint:: " << std::endl;
	sstr << "    bodies connected = (" << body_name_1 << ", " << body_name_2 << ")" << std::endl;
	sstr << "    is reversed? " << (isReversed() ? "yes" : "no") << std::endl;
	sstr << "    is cut? " << (isCut() ? "yes" : "no") << std::endl;
	sstr << "    is prescribed? " << (isPrescribed() ? "yes" : "no") << std::endl;
	sstr << "    d.o.f. = " << getDOF() << std::endl;
	sstr << "    q = {";
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		sstr << (*iter_pcoord)->q << ", ";
	}
	sstr << "}" << std::endl;
	sstr << "    dq = {";
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		sstr << (*iter_pcoord)->dq << ", ";
	}
	sstr << "}" << std::endl;
	sstr << "    ddq = {";
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		sstr << (*iter_pcoord)->ddq << ", ";
	}
	sstr << "}" << std::endl;
	sstr << "    tau = {";
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		sstr << (*iter_pcoord)->tau << ", ";
	}
	sstr << "}" << std::endl;
	sstr << "T_global = " << T_global;
	sstr << "T_left = " << T_left;
	sstr << "T_right = " << T_right;
	sstr << "T = " << T << std::endl;
	sstr << "inv_T = " << inv_T << std::endl;
	sstr << "S = " << S << std::endl;
	sstr << "dS = " << dS << std::endl;
	sstr << "Sdq = " << Sdq << std::endl;
	sstr << "dSdq = " << dSdq << std::endl;
	sstr << "Sddq = " << Sddq << std::endl;
	sstr << "DSdqDt = " << DSdqDt << std::endl;
	sstr << std::endl;

	return sstr.str();
}

se3 GJoint::get_S(int idx_)
{
	if ( idx_ < 0 || idx_ >= getDOF() ) {
		return se3(0,0,0,0,0,0);
	} else {
		//return se3(S(0,idx_),S(1,idx_),S(2,idx_),S(3,idx_),S(4,idx_),S(5,idx_));
		return se3(&S[6*idx_]);
	}
}

se3 GJoint::get_dS(int idx_)
{
	if ( idx_ < 0 || idx_ >= getDOF() ) {
		return se3(0,0,0,0,0,0);
	} else {
		//return se3(dS(0,idx_),dS(1,idx_),dS(2,idx_),dS(3,idx_),dS(4,idx_),dS(5,idx_));
		return se3(&dS[6*idx_]);
	}
}

se3 GJoint::get_S(GCoordinate *pCoordinate_)
{
	return get_S(getIndexOfCoordinate(pCoordinate_));
}

se3 GJoint::get_dS(GCoordinate *pCoordinate_)
{
	return get_dS(getIndexOfCoordinate(pCoordinate_));
}

void GJoint::allocate_memory(int n_)
{
	S.SetZero(6,n_);
	dS.SetZero(6,n_);
}
