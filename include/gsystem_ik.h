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
//         GSystemIK: class for systems with IK
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_SYSTEM_WITH_INVERSE_KINEMATICS_
#define _GEAR_SYSTEM_WITH_INVERSE_KINEMATICS_

#include <vector>
#include "gsystem.h"
#include "gbody.h"
#include "liegroup.h"



//=============================================================
//                 GSystemIK
//=============================================================
class GSystemIK: public GSystem
{
public:
	GSystemIK() {}
	~GSystemIK() {}

public:

	// solve q = pCoordinates[]->q satisfying IK constraints with Newton-Raphson method using analytic Jacobian
	//  - pbodies[i] = pointers to the bodies considered
	//  - p[i] = location of the origin of a local frame attached to pbodies[i] 
	//           ** the local frame can be regarded as an end-effector frame
	//           ** the local frame can be different from {pbodies[i]}.
	//           ** the orientation of the frame is assumed to be aligned to {pbodies[i]}
	//  - T_target[i] = target transform of the local frame of pbodies[i] w.r.t. {global}
	//  - idxC[i] = constraint setting for pbodies[i]
	//	         e.g., by setting idxC[i] = [0, 1, 2, 5], all three components of the orientation
	//		     and the z-axis location of the local body frame will be constrained by T_target[i].
	//  - tol = tolerance for constraint violation
	//  - max_iter_num = maximum number of iteration allowed
	//bool solveIK_q(
	//	RMatrix &q,
	//	std::vector<GBody*> pbodies,
	//	std::vector<Vec3> p,
	//	std::vector<SE3> T_target,
	//	std::vector< std::vector<int> > idxC,
	//	gReal tol = 1E-6, int max_iter_num = 50);


	// solve dq = pCoordinates[]->dq satisfying the primary goal and minimizing the error on the secondary goal.
	//	- primary goal: V_primary[i] = pbodies_primary[i]'s (angular velocity, linear velocity at p_primary[i]) w.r.t. {global}
	//	- secondary goal: V_secondary[i] = pbodies_secondary[i]'s (angular velocity, linear velocity at p_secondary[i]) w.r.t. {global}
	//	- p_primary[i] and p_secondary[i] are position vectors w.r.t. {pbodies_primary[i]} and {pbodies_secondary[i]} respectively.
	//  - idxC_primary[i]: constraint setting for the i-th primary body
	//	- idxC_secondary[i]: constraint setting for the i-th secondary body
	//	    e.g., by setting idxC_primary[i] = [0, 1, 2, 5], all three components of the angular velocity 
	//		and the z-axis component of the linear velocity of pbodies_primary[i] will be constrained.
	//	- pcoords_prescribed: the coordinates whose velocities(dq) are already prescribed. (pcoords_prescribed[]->dq will be remained in the solution)
	//	- alpha: a parameter for singularity-robust inverse, srInv()
	bool solveIK_dq(
		RMatrix &dq, 
		std::vector<GBody*> pbodies_primary, std::vector<GBody*> pbodies_secondary, 
		std::vector<Vec3> p_primary, std::vector<Vec3> p_secondary,
		std::vector<se3> V_primary, std::vector<se3> V_secondary, 
		std::vector< std::vector<int> > idxC_primary, std::vector< std::vector<int> > idxC_secondary, 
		gReal alpha_primary = 0, gReal alpha_secondary = 0.001);
	bool solveIK_dq(
		RMatrix &dq, 
		std::vector<GBody*> pbodies_primary, std::vector<GBody*> pbodies_secondary, 
		std::vector<Vec3> p_primary, std::vector<Vec3> p_secondary,
		std::vector<se3> V_primary, std::vector<se3> V_secondary, 
		std::vector< std::vector<int> > idxC_primary, std::vector< std::vector<int> > idxC_secondary, 
		std::vector<GCoordinate*> pcoords_prescribed,
		gReal alpha_primary = 0, gReal alpha_secondary = 0.001);
	bool solveIK_dq(
		RMatrix &dq, 
		std::vector<GBody*> pbodies_primary, std::vector<GBody*> pbodies_secondary, 
		std::vector<Vec3> p_primary, std::vector<Vec3> p_secondary,
		std::vector<se3> V_primary, std::vector<se3> V_secondary, 
		std::vector< std::vector<int> > idxC_primary, std::vector< std::vector<int> > idxC_secondary, 
		std::vector<GCoordinate*> pcoords_prescribed,
		std::ofstream *pfout,
		gReal alpha_primary = 0, gReal alpha_secondary = 0.001);

	// build IK constraints J*dq = V where dq = pCoordinates[]->dq
	//	- IK constraints: V_target[i] = pbodies[i]'s (angular velocity, linear velocity at pos[i]) w.r.t. {global}
	//	- pos[i] is a position std::vector w.r.t. {pbodies[i]}
	//	- idxC[i]: index of the active constraints on the i-th body
	//	    e.g., by setting idxC[i] = [0, 1, 2, 5], all three components of the angular velocity 
	//		and the z-axis component of the linear velocity of pbodies[i] will be constrained.
	bool buildConstrIK_dq(
		RMatrix &J, RMatrix &V,
		std::vector<GBody*> pbodies, std::vector<Vec3> pos, std::vector<se3> V_target, std::vector< std::vector<int> > idxC);

};



#endif

