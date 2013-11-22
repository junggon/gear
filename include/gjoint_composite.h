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
//         GJointComposite: class for composite-type joints
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_JOINT_COMPOSITE_
#define _GEAR_JOINT_COMPOSITE_

#include "gjoint.h"
#include "gcoordinate.h"
#include "liegroup.h"
#include "rmatrix3j.h"



//=============================================================
//                 GJointComposite
//=============================================================
//
// T = pJoint1->T * pJoint2->T
//
// ** Computational speed of GJointComposite may be slower than that of an optimized single joint implementation.
// For e.g., a spherical joint can be implemented by GJointComposite(GJointComposite(GJointRevolute, GJointRevolute), GJointRevolute),
// and the calculation speed for update() is more than two times slower than that of GJointSpherical.
//
class GJointComposite: public GJoint
{
public:
	GJoint *pJoint1, *pJoint2;
	
public:
	GJointComposite();
	GJointComposite(GJoint *pjoint1_, GJoint *pjoint2_);
	~GJointComposite() {}

public:
	bool compose(GJoint *pjoint1_, GJoint *pjoint2_);	// compose two joints to make a composite joint
	
public:
	bool isConstantScrew() { return false; }

	void update_short();
	void update();

	void _update_short_for_reversed_joint();			// modify T, inv_T, S for reversed joint
	void _update_for_reversed_joint();					// modify T, inv_T, S, dS, Sdq, dSdq, Sddq, DSdqDt for reversed joint

	RMatrix get_DSDq(GCoordinate *pCoordinate_);
	RMatrix get_DdSDq(GCoordinate *pCoordinate_);
};



#endif

