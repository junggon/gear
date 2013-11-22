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
#include "gjoint_revolute.h"
#include "gjoint.h"
#include "liegroup.h"
#include "rmatrix3j.h"



//=============================================================
//                 GJointRevolute
//=============================================================
GJointRevolute::GJointRevolute()
{
	jointType = GJOINT_REVOLUTE;
	axis = Vec3(0,0,1);							// default axis
	pCoordinates.push_back(&coordinate);		// set coordinate pointer
	allocate_memory(1);
}

void GJointRevolute::setAxis(gReal x_, gReal y_, gReal z_)
{
	axis[0] = x_;
	axis[1] = y_;
	axis[2] = z_;
}

void GJointRevolute::update_short()
{
	if ( bReversed ) {
		T = SE3(Exp(-axis*coordinate.q), Vec3(0,0,0));
		inv_T = SE3(~T.GetRotation());
		S[0] = -axis[0]; S[1] = -axis[1]; S[2] = -axis[2];
	} else {
		T = SE3(Exp(axis*coordinate.q), Vec3(0,0,0));
		inv_T = SE3(~T.GetRotation());
		S[0] = axis[0]; S[1] = axis[1]; S[2] = axis[2];
	}
}

void GJointRevolute::update()
{
	//dS, dSdq are still zeros.

	if ( bReversed ) {
		T = SE3(Exp(-axis*coordinate.q), Vec3(0,0,0));
		inv_T = SE3(~T.GetRotation());
		Sdq = se3(-axis*coordinate.dq, Vec3(0,0,0));
		Sddq = se3(-axis*coordinate.ddq, Vec3(0,0,0));
		DSdqDt = Sddq;
		S[0] = -axis[0]; S[1] = -axis[1]; S[2] = -axis[2];
	} else {
		T = SE3(Exp(axis*coordinate.q), Vec3(0,0,0));
		inv_T = SE3(~T.GetRotation());
		Sdq = se3(axis*coordinate.dq, Vec3(0,0,0));
		Sddq = se3(axis*coordinate.ddq, Vec3(0,0,0));
		DSdqDt = Sddq;
		S[0] = axis[0]; S[1] = axis[1]; S[2] = axis[2];
	}
}

RMatrix GJointRevolute::get_DSDq(GCoordinate *pCoordinate_)
{
	return Zeros(6,1);
}

RMatrix GJointRevolute::get_DdSDq(GCoordinate *pCoordinate_)
{
	return Zeros(6,1);
}

