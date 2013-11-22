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

//================================================================================
//         GConstraintJointLoop: class for joint-loop constraints
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_CONSTRAINT_JOINTLOOP_
#define _GEAR_CONSTRAINT_JOINTLOOP_

#include <list>
#include "gconstraint.h"
#include "liegroup.h"

class GJoint;

//=============================================================
//                 GConstraintJointLoop
//=============================================================
class GConstraintJointLoop: public GConstraint
{
public:
	enum JOINTLOOP_CONSTRAINT_TYPE
	{
		JOINTLOOP_ORIENTATION_POSITION,		// both position and orientation considered
		JOINTLOOP_ORIENTATION_ONLY,			// only orientation considered
		JOINTLOOP_POSITION_ONLY,			// only position considered
	};

	// constraint: f(joints)*M1 = f(joints2)*M2, if pJoints2.size() > 0
	//             f(joints)*M1 = T, if pJoints2.size() = 0
	//             Note: Joints in pJoints and pJoints2 must be unique.
	//                   Cut-joint = pJoints[pJoints.size()-1]
	//					 M1 = M2 = T = Eye by default.
	std::list<GJoint *> pJoints, pJoints2;
	SE3 M1, M2;
	SE3 T;

	int num_coord, num_coord2;

	JOINTLOOP_CONSTRAINT_TYPE jointLoopConstraintType;	// JOINTLOOP_ORIENTATION_POSITION(default) | JOINTLOOP_ORIENTATION_ONLY | JOINTLOOP_POSITION_ONLY

	RMatrix jacobian;	// temporary Jacobian matrix (update_dJdt() needs this.)

public:
	GConstraintJointLoop();
	~GConstraintJointLoop() {}

public:
	bool setJoints(std::list<GJoint *> pjoints);
	bool setJoints(std::list<GJoint *> pjoints, std::list<GJoint *> pjoints2);
	bool setM(const SE3 &M1_, const SE3 &M2_ = SE3()) { M1 = M1_; M2 = M2_; return true; }
	bool setT(const SE3 &T_) { M2 = T_; return true; }
	void setJointLoopConstraintType(JOINTLOOP_CONSTRAINT_TYPE jointLoopConstraintType_);

	SE3 getLoopSE3();	// return f(joints)*M1

	// virtual functions
	void update_C();		// In addition to C, pJoints[]->(T,inv_T,S) will be updated in this function.
	void update_J();		// pJoints[]->(inv_T, S) need to be updated before calling this function.
	void update_dJdt();		// update_J() need to be called before calling this function.
	std::string getInfoStr();
};



#endif

