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

#include "gconstraint_jointloop.h"
#include "gconstraint.h"
#include "gjoint.h"
#include "gcoordinate.h"

#include "liegroup.h"
#include "rmatrix3j.h"


//=============================================================
//                 GConstraintJointLoop
//=============================================================
GConstraintJointLoop::GConstraintJointLoop()
{
	M1.SetIdentity();
	M2.SetIdentity();
	T.SetIdentity();

	num_coord = 0;
	num_coord2 = 0;

	jointLoopConstraintType = JOINTLOOP_ORIENTATION_POSITION;

	jacobian.SetZero(0,0);
}

void GConstraintJointLoop::setJointLoopConstraintType(JOINTLOOP_CONSTRAINT_TYPE jointLoopConstraintType_)
{
	jointLoopConstraintType = jointLoopConstraintType_;

	switch ( jointLoopConstraintType )
	{
	case JOINTLOOP_ORIENTATION_POSITION:
		constrNum = 6;
		break;
	case JOINTLOOP_ORIENTATION_ONLY:
		constrNum = 3;
		break;
	case JOINTLOOP_POSITION_ONLY:
		constrNum = 3;
		break;
	}
}

bool GConstraintJointLoop::setJoints(std::list<GJoint *> pjoints)
{
	std::list<GJoint *>::iterator iter_pjoint;
	std::list<GCoordinate *>::iterator iter_pcoord;

	pCoordinates.clear();
	pJoints.clear();
	pJoints2.clear();
	num_coord = 0;
	num_coord2 = 0;

	for (iter_pjoint = pjoints.begin(); iter_pjoint != pjoints.end(); iter_pjoint++) {
		pJoints.push_back(*iter_pjoint);
		for (iter_pcoord = (*iter_pjoint)->pCoordinates.begin(); iter_pcoord != (*iter_pjoint)->pCoordinates.end(); iter_pcoord++) {
			pCoordinates.push_back(*iter_pcoord);
		}
		num_coord += (*iter_pjoint)->getDOF();
	}
	
	jointLoopConstraintType = JOINTLOOP_ORIENTATION_POSITION;
	constrNum = 6;

	return true;
}

bool GConstraintJointLoop::setJoints(std::list<GJoint *> pjoints, std::list<GJoint *> pjoints2)
{
	std::list<GJoint *>::iterator iter_pjoint;
	std::list<GCoordinate *>::iterator iter_pcoord;

	pCoordinates.clear();
	pJoints.clear();
	pJoints2.clear();
	num_coord = 0;
	num_coord2 = 0;

	for (iter_pjoint = pjoints.begin(); iter_pjoint != pjoints.end(); iter_pjoint++) {
		pJoints.push_back(*iter_pjoint);
		for (iter_pcoord = (*iter_pjoint)->pCoordinates.begin(); iter_pcoord != (*iter_pjoint)->pCoordinates.end(); iter_pcoord++) {
			pCoordinates.push_back(*iter_pcoord);
		}
		num_coord += (*iter_pjoint)->getDOF();
	}

	for (iter_pjoint = pjoints2.begin(); iter_pjoint != pjoints2.end(); iter_pjoint++) {
		pJoints2.push_back(*iter_pjoint);
		for (iter_pcoord = (*iter_pjoint)->pCoordinates.begin(); iter_pcoord != (*iter_pjoint)->pCoordinates.end(); iter_pcoord++) {
			pCoordinates.push_back(*iter_pcoord);
		}
		num_coord2 += (*iter_pjoint)->getDOF();
	}

	if ( getNumCoordinates() != num_coord + num_coord2 ) return false;

	jointLoopConstraintType = JOINTLOOP_ORIENTATION_POSITION;
	constrNum = 6;

	return true;
}

SE3 GConstraintJointLoop::getLoopSE3()
{
	std::list<GJoint *>::iterator iter_pjoint;
	SE3 re;

	re.SetIdentity();
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		re *= (*iter_pjoint)->T_left;
		re *= (*iter_pjoint)->T;
		re *= (*iter_pjoint)->inv_T_right;
	}
	re *= M1;
	
	return re;
}

void GConstraintJointLoop::update_C()
{
	std::list<GJoint *>::iterator iter_pjoint;
	SE3 T_loop_left, T_loop_right; 
	gReal re[6];
	
	// left side of the loop
	T_loop_left.SetIdentity();
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		(*iter_pjoint)->update_short();
		T_loop_left *= (*iter_pjoint)->T_left;
		T_loop_left *= (*iter_pjoint)->T;
		T_loop_left *= (*iter_pjoint)->inv_T_right;
	}
	T_loop_left *= M1;

	// right side of the loop
	T_loop_right.SetIdentity();
	for (iter_pjoint = pJoints2.begin(); iter_pjoint != pJoints2.end(); iter_pjoint++) {
		(*iter_pjoint)->update_short();
		T_loop_right *= (*iter_pjoint)->T_right;
		T_loop_right *= (*iter_pjoint)->T;
		T_loop_right *= (*iter_pjoint)->inv_T_right;
	}
	T_loop_right *= M2;

	//put_se3_to_matrix(re, InvSkew(Inv(T_loop_right)*T_loop_left), 0);		// InvSkew(G) = unskew(G-I), where G-I is assumed to be a 4x4 se3.
	matSet(re, InvSkew(Inv(T_loop_right)*T_loop_left).GetArray(), 6);		// InvSkew(G) = unskew(G-I), where G-I is assumed to be a 4x4 se3.
																// ** Modified in 2007.05.22: 
																//    InvSkew(Inv(T_right)*T_left) --> InvSkew(Inv(T_left)*T_right)
																// ** Modification canceled in 2007.05.22, i.e.,
																//    InvSkew(Inv(T_left)*T_right) --> InvSkew(Inv(T_right)*T_left)

	C.SetZero(constrNum, 1);
	switch ( jointLoopConstraintType )
	{
	case JOINTLOOP_ORIENTATION_POSITION:
		matSet(C.GetPtr(), re, 6); // C = re
		break;
	case JOINTLOOP_ORIENTATION_ONLY:
		matSet(C.GetPtr(), re, 3); // C = re[0:2]
		break;
	case JOINTLOOP_POSITION_ONLY:
		matSet(C.GetPtr(), &re[3], 3); // C = re[3:5]
		break;
	}
}

// body J of the joint loop constraint
// case f=g   : J = [ inv(f)*dot(f), -inv(g)*dot(g) ]
// case f*M=T : J = [ inv(f*M)*dot(f*M) ] = Ad( inv(M), [ inv(f)*dot(f) ] )
void GConstraintJointLoop::update_J()
{
	int idx; 
	SE3 Ti;
	std::list<GJoint *>::reverse_iterator riter_pjoint;
	
	jacobian.SetZero(6, getNumCoordinates());

	idx = num_coord;
	Ti.SetInvOf(M1);
	for (riter_pjoint = pJoints.rbegin(); riter_pjoint != pJoints.rend(); riter_pjoint++) {
		idx -= (*riter_pjoint)->getDOF();
		Ti *= (*riter_pjoint)->T_right;
		//jacobian.Push(0, idx, Ad(Ti, (*riter_pjoint)->S));
		matSet_Ad(&jacobian[6*idx], Ti, (*riter_pjoint)->S.GetPtr(), (*riter_pjoint)->getDOF());
		Ti *= (*riter_pjoint)->inv_T;
		Ti *= (*riter_pjoint)->inv_T_left;
	}

	idx = num_coord + num_coord2;
	Ti.SetInvOf(M2);
	for (riter_pjoint = pJoints2.rbegin(); riter_pjoint != pJoints2.rend(); riter_pjoint++) {
		idx -= (*riter_pjoint)->getDOF();
		Ti *= (*riter_pjoint)->T_right;
		//jacobian.Push(0, idx, -Ad(Ti, (*riter_pjoint)->S));
		matSet_Ad_minus(&jacobian[6*idx], Ti, (*riter_pjoint)->S.GetPtr(), (*riter_pjoint)->getDOF());
		Ti *= (*riter_pjoint)->inv_T;
		Ti *= (*riter_pjoint)->inv_T_left;
	}
	
	J.SetZero(constrNum, getNumCoordinates());
	switch ( jointLoopConstraintType )
	{
	case JOINTLOOP_ORIENTATION_POSITION:
		J = jacobian;
		break;
	case JOINTLOOP_ORIENTATION_ONLY:
		J = jacobian.Sub(0, 2, 0, getNumCoordinates()-1);
		break;
	case JOINTLOOP_POSITION_ONLY:
		J = jacobian.Sub(3, 5, 0, getNumCoordinates()-1);
		break;
	}
}

void GConstraintJointLoop::update_dJdt()
{
	int i, j;
	se3 Jsum, Jtmp; // RMatrix Jsum;
	std::list<GCoordinate *>::iterator iter_pcoord_i, iter_pcoord_j;

	RMatrix dotjacobian(6, getNumCoordinates());
	dotjacobian.SetZero();

	i=0;
	iter_pcoord_i = pCoordinates.begin();

	for ( ; i<num_coord; i++)
	{
		Jsum.SetZero(); // Jsum = Zeros(6,1);
		iter_pcoord_j = iter_pcoord_i;
		iter_pcoord_j++;
		for (j=i+1; j<num_coord; j++)
		{
			// Jsum += ad(get_jacobian(i), get_jacobian(j)) * (*iter_pcoord_j)->dq;
			Jtmp.set_ad(&jacobian[6*i], &jacobian[6*j]); // &jacobian[6*k] --> pointer to the k-th column of jacobian
			Jtmp *= (*iter_pcoord_j)->dq;
			Jsum += Jtmp;

			iter_pcoord_j++;
		}
		matSet(&dotjacobian[6*i], Jsum.GetArray(), 6); // dotjacobian.Push(0, i, Jsum);
		iter_pcoord_i++;
	}

	for ( ; i<num_coord+num_coord2; i++)
	{
		Jsum.SetZero(); // Jsum = Zeros(6,1);
		iter_pcoord_j = iter_pcoord_i;
		iter_pcoord_j++;
		for (j=i+1; j<num_coord+num_coord2; j++)
		{
			// Jsum += ad(get_jacobian(i), get_jacobian(j)) * (*iter_pcoord_j)->dq;
			Jtmp.set_ad(&jacobian[6*i], &jacobian[6*j]);
			Jtmp *= (*iter_pcoord_j)->dq;
			Jsum -= Jtmp; // cautious! '-=' used here 

			iter_pcoord_j++;
		}
		// dotjacobian.Push(0, i, -Jsum);
		matSet(&dotjacobian[6*i], Jsum.GetArray(), 6); // cautious! see above where '-=' used!
		iter_pcoord_i++;
	}

	dJdt.SetZero(constrNum, getNumCoordinates());
	switch ( jointLoopConstraintType )
	{
	case JOINTLOOP_ORIENTATION_POSITION:
		dJdt = dotjacobian;
		break;
	case JOINTLOOP_ORIENTATION_ONLY:
		dJdt = dotjacobian.Sub(0, 2, 0, getNumCoordinates()-1);
		break;
	case JOINTLOOP_POSITION_ONLY:
		dJdt = dotjacobian.Sub(3, 5, 0, getNumCoordinates()-1);
		break;
	}
}

std::string GConstraintJointLoop::getInfoStr()
{
	std::stringstream sstr;
	std::list<GJoint*>::iterator iter_pjoint;

	sstr << GElement::getInfoStr();
	sstr << "jointLoopConstraintType = " << jointLoopConstraintType << std::endl;
	sstr << "pJoints = ";
	for (iter_pjoint = pJoints.begin(); iter_pjoint != pJoints.end(); iter_pjoint++) {
		sstr << (*iter_pjoint)->getName() << ", ";
	}
	sstr << std::endl;
	sstr << "pJoints2 = ";
	for (iter_pjoint = pJoints2.begin(); iter_pjoint != pJoints2.end(); iter_pjoint++) {
		sstr << (*iter_pjoint)->getName() << ", ";
	}
	sstr << std::endl;
	sstr << "M1 = " << M1 << "M2 = " << M2 << "T = " << T << std::endl;

	return sstr.str();
}