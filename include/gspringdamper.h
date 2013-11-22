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
//         GSpringDamperBody: spring-damper connecting two rigid bodies
//         GSpringDamperJoint: spring-damper for joints
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_SPRING_DAMPER_
#define _GEAR_SPRING_DAMPER_

#include "gforce.h"
#include "gbody.h"
#include "gjoint.h"


//=============================================================
//                 GSpringDamperBody
//=============================================================
class GSpringDamperBody : public GForce
{
public:
	GBody *pLeftBody, *pRightBody;		// two rigid bodies connected by the spring-damper

	SE3 T_left, T_right;				// locations of both ends of the spring-damper w.r.t. {pLeftBody} and {pRightBody}
	SE3 inv_T_left, inv_T_right;		// inverse of T_left and T_right
	
	gReal K[6], C[6];					// spring and damping coefficients on each screw component

public:
	GSpringDamperBody();
	~GSpringDamperBody();

public:
	virtual void clear();

	bool connectBodies(GBody *pLeftBody_, GBody *pRightBody_);
	bool disconnectBodies();

	void setPosition(const Vec3 &pL_, const Vec3 &pR_);
	void setOrientation(const SO3 &RL_, const SO3 &RR_);
	void setPositionAndOrientation(const SE3 &TL_, const SE3 &TR_);

	void setSpring(const Vec3 &K_rotation_, const Vec3 &K_position_);
	void setDamper(const Vec3 &C_rotation_, const Vec3 &C_position_);

public:
	bool applyForce(bool badd_ = false);		// If badd = true, adds spring-damper force to pLeftBody->Fe and pRightBody->Fe. If badd = false, subtracts the force.
};

//=============================================================
//                 GSpringDamperJoint
//=============================================================
class GSpringDamperJoint : public GForce
{
public:
	GJoint *pJoint;

	gReal *q_neutral;
	gReal *K, *C;

public:
	GSpringDamperJoint();
	~GSpringDamperJoint();

public:
	virtual void clear();

	bool connectJoint(GJoint *pjoint_);
	bool disconnectJoint();
	
	bool setNeutralPosition(gReal *q_neutral_ = NULL);	// if q_neutral_ == NULL, set current pJoint->pCoordinates[]->q as the neutral position.
	bool setNeutralPosition(RMatrix q_neutral_);
	bool setNeutralPosition(gReal q_neutral_);

	bool setSpring(gReal *K_);
	bool setSpring(RMatrix K_);
	bool setSpring(gReal K_);

	bool setDamper(gReal *C_);
	bool setDamper(gReal C_);
	bool setDamper(RMatrix C_);

public:
	bool applyForce(bool badd_ = false);		// If badd = true, adds spring-damper force (or torque) to joint torques. If badd = false, subtracts the force.
};



#endif

