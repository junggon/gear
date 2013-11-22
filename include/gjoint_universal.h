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
//         GJointUniversal: class for universal joints (X-Y type)
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_JOINT_UNIVERSAL_
#define _GEAR_JOINT_UNIVERSAL_

#include "gcoordinate.h"
#include "gjoint.h"
#include "gjoint_composite.h"
#include "gjoint_revolute.h"
#include "liegroup.h"
#include "rmatrix3j.h"

class GJointUniversalXY;
class GJointUniversalGeneral;

typedef GJointUniversalXY GJointUniversal;

//=============================================================
//                 GJointUniversalXY
//=============================================================
//
//  T = SE3(R(x,q0) * R(y,q1), Vec3(0,0,0)) 
//  where q0 = coordinates[0].q and q1 = coordinates[1].q 
//
class GJointUniversalXY: public GJoint
{
public:
	GCoordinate coordinates[2];	

public:
	GJointUniversalXY();
	~GJointUniversalXY() {}

public:
	void getAxes(Vec3 &axis1_, Vec3 &axis2_);	// return the x, y axes w.r.t. {joint left}
												// i.e., axis1_ = (1,0,0), axis2_ = R(x,q0)*Vec3(0,1,0)

public:
	bool isConstantScrew() { return false; }

	void update_short();
	void update();

	RMatrix get_DSDq(GCoordinate *pCoordinate_);
	RMatrix get_DdSDq(GCoordinate *pCoordinate_);
};

//=============================================================
//                 GJointUniversalGeneral
//=============================================================
//
//  T = SE3(R(axis1,q1) * R(axis2,q2), Vec3(0,0,0)) 
//
class GJointUniversalGeneral: public GJointComposite
{
public:
	GJointRevolute rjoint1, rjoint2;

public:
	GJointUniversalGeneral() { rjoint1.setAxis(1,0,0); rjoint2.setAxis(0,1,0); compose(&rjoint1, &rjoint2); jointType = GJOINT_UNIVERSAL_GENERAL; }
	GJointUniversalGeneral(const Vec3 &axis1, const Vec3 &axis2) { rjoint1.setAxis(axis1); rjoint2.setAxis(axis2); compose(&rjoint1, &rjoint2); jointType = GJOINT_UNIVERSAL_GENERAL; }
	~GJointUniversalGeneral() {}

	void setAxes(const Vec3 &axis1, const Vec3 &axis2) { rjoint1.setAxis(axis1); rjoint2.setAxis(axis2); }
	void getAxes(Vec3 &axis1, Vec3 &axis2) { axis1 = rjoint1.getAxis(); axis2 = rjoint2.getAxis(); }
};



#endif

