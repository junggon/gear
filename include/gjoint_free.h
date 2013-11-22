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
//         GJointFree: class for free joints
//           (variants: GJointFreeST(=GJointFree), GJointFreeTS)
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_JOINT_FREE_
#define _GEAR_JOINT_FREE_

#include "gjoint_composite.h"
#include "gjoint_spherical.h"
#include "gjoint_translational.h"
#include "liegroup.h"
#include "rmatrix3j.h"

class GJointFreeST; // spherical + translational
class GJointFreeTS; // translational + spherical

typedef GJointFreeST GJointFree;

//=============================================================
//                 GJointFreeST
//=============================================================
class GJointFreeST: public GJoint
{
public:
	GJointSpherical spherical_joint;
	GJointTranslational translational_joint;

public:
	GJointFreeST();
	~GJointFreeST() {}

public:
	bool isConstantScrew() { return false; }

	void update_short();
	void update();

	void _update_short_for_reversed_joint();	// modify T, inv_T, S for reversed joint
	void _update_for_reversed_joint();			// modify T, inv_T, S, dS, Sdq, dSdq, Sddq, DSdqDt for reversed joint

	RMatrix get_DSDq(GCoordinate *pCoordinate_);
	RMatrix get_DdSDq(GCoordinate *pCoordinate_);

	int getCoordinateChartForRotation() { return spherical_joint.getCoordinateChart(); }
	void setFixedCoordinateChartForRotation(GJointSpherical::CoordinateChartForSphericalJoint cc_) { spherical_joint.setFixedCoordinateChart(cc_); }
	void resetCoordinateChartForRotation() { spherical_joint.resetCoordinateChart(); }

	// set coordinates[]->(q,dq,ddq) with given set of (T, dT, ddT) or (T, V, dV)
	// T = SE(3): {joint left} --> {joint right}
	// dT = time derivative of T (4x4 matrix)
	// ddT = time derivative of dT (4x4 matrix)
	// V = Inv(T)*dT/dt = generalized velocity of {joint right} w.r.t. {joint left} viewed in {joint right}
	// dV = dV/dt
	void setMotion(const SE3 &T, const RMatrix &dT, const RMatrix &ddT);
	void setMotion(const SE3 &T, const se3 &V, const se3 &dV);
};

//=============================================================
//                 GJointFreeTS
//=============================================================
class GJointFreeTS: public GJointComposite
{
public:
	GJointTranslational translational_joint;
	GJointSpherical spherical_joint;

public:
	GJointFreeTS() { compose(&translational_joint, &spherical_joint); jointType = GJOINT_FREE_TS;}
	~GJointFreeTS() {}

	int getCoordinateChartForRotation() { return spherical_joint.getCoordinateChart(); }
	void setFixedCoordinateChartForRotation(GJointSpherical::CoordinateChartForSphericalJoint cc_) { spherical_joint.setFixedCoordinateChart(cc_); }
	void resetCoordinateChartForRotation() { spherical_joint.resetCoordinateChart(); }

	// set coordinates[]->(q,dq,ddq) with given set of (T, dT, ddT) or (T, V, dV)
	// T = SE(3): {joint left} --> {joint right}
	// dT = time derivative of T (4x4 matrix)
	// ddT = time derivative of dT (4x4 matrix)
	// V = Inv(T)*dT/dt = generalized velocity of {joint right} w.r.t. {joint left} viewed in {joint right}
	// dV = dV/dt
	void setMotion(const SE3 &T, const RMatrix &dT, const RMatrix &ddT);
	void setMotion(const SE3 &T, const se3 &V, const se3 &dV);
};

//=============================================================
//                 GJointFreeST2 (Just for testing GJointComposite. Use GJointFree (=GJointFreeST) instead.)
//=============================================================
class GJointFreeST2: public GJointComposite
{
public:
	GJointSpherical spherical_joint;
	GJointTranslational translational_joint;

public:
	GJointFreeST2() { compose(&spherical_joint, &translational_joint); }
	~GJointFreeST2() {}
};


#endif

