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
#include "gjoint_prismatic.h"
#include "gjoint.h"
#include "liegroup.h"
#include "rmatrix3j.h"



//=============================================================
//                 GJointPrismatic
//=============================================================
GJointPrismatic::GJointPrismatic()
{
	jointType = GJOINT_PRISMATIC;
	pCoordinates.push_back(&coordinate);		// set coordinate pointer
	axis = Vec3(0,0,1);							// default axis
	allocate_memory(1);
}

void GJointPrismatic::setAxis(gReal x_, gReal y_, gReal z_)
{
	axis[0] = x_;
	axis[1] = y_;
	axis[2] = z_;
}

void GJointPrismatic::update_short()
{
	if ( bReversed ) {
		T = SE3(-axis*coordinate.q);
		inv_T = SE3(-T.GetPosition());
		S[3] = -axis[0]; S[4] = -axis[1]; S[5] = -axis[2];
	} else {
		T = SE3(axis*coordinate.q);
		inv_T = SE3(-T.GetPosition());
		S[3] = axis[0]; S[4] = axis[1]; S[5] = axis[2];
	}
}

void GJointPrismatic::update()
{
	//dS, dSdq are still zeros.

	if ( bReversed ) {
		T = SE3(-axis*coordinate.q);
		inv_T = SE3(-T.GetPosition());
		Sdq = se3(Vec3(0,0,0), -axis*coordinate.dq);
		Sddq = se3(Vec3(0,0,0), -axis*coordinate.ddq);
		DSdqDt = Sddq;
		S[3] = -axis[0]; S[4] = -axis[1]; S[5] = -axis[2];
	} else {
		T = SE3(axis*coordinate.q);
		inv_T = SE3(-T.GetPosition());
		Sdq = se3(Vec3(0,0,0), axis*coordinate.dq);
		Sddq = se3(Vec3(0,0,0), axis*coordinate.ddq);
		DSdqDt = Sddq;
		S[3] = axis[0]; S[4] = axis[1]; S[5] = axis[2];
	}
}

RMatrix GJointPrismatic::get_DSDq(GCoordinate *pCoordinate_)
{
	return Zeros(6,1);
}

RMatrix GJointPrismatic::get_DdSDq(GCoordinate *pCoordinate_)
{
	return Zeros(6,1);
}

