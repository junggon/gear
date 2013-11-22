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
#include "gjoint_planar.h"
#include "gjoint.h"
#include "liegroup.h"
#include "rmatrix3j.h"



//=============================================================
//                 GJointPlanar
//=============================================================
GJointPlanar::GJointPlanar()
{
	jointType = GJOINT_PLANAR;
	pCoordinates.push_back(&coordinates[0]);
	pCoordinates.push_back(&coordinates[1]);
	allocate_memory(2);
}

void GJointPlanar::update_short()
{
	if ( bReversed ) {
		T = SE3(1,0,0,0,1,0,0,0,1, -coordinates[0].q, -coordinates[1].q, 0);
		inv_T = SE3(1,0,0,0,1,0,0,0,1, coordinates[0].q, coordinates[1].q, 0);
		S[3] = S[10] = -1.0;
	} else {
		T = SE3(1,0,0,0,1,0,0,0,1, coordinates[0].q, coordinates[1].q, 0);
		inv_T = SE3(1,0,0,0,1,0,0,0,1, -coordinates[0].q, -coordinates[1].q, 0);
		S[3] = S[10] = 1.0;
	}
}

void GJointPlanar::update()
{
	//dS, dSdq are still zeros.

	if ( bReversed ) {
		T = SE3(1,0,0,0,1,0,0,0,1, -coordinates[0].q, -coordinates[1].q, 0);
		inv_T = SE3(1,0,0,0,1,0,0,0,1, coordinates[0].q, coordinates[1].q, 0);
		Sdq = se3(0, 0, 0, -coordinates[0].dq, -coordinates[1].dq, 0);
		Sddq = se3(0, 0, 0, -coordinates[0].ddq, -coordinates[1].ddq, 0);
		DSdqDt = Sddq;
		S[3] = S[10] = -1.0;
	} else {
		T = SE3(1,0,0,0,1,0,0,0,1, coordinates[0].q, coordinates[1].q, 0);
		inv_T = SE3(1,0,0,0,1,0,0,0,1, -coordinates[0].q, -coordinates[1].q, 0);
		Sdq = se3(0, 0, 0, coordinates[0].dq, coordinates[1].dq, 0);
		Sddq = se3(0, 0, 0, coordinates[0].ddq, coordinates[1].ddq, 0);
		DSdqDt = Sddq;
		S[3] = S[10] = 1.0;
	}
}

RMatrix GJointPlanar::get_DSDq(GCoordinate *pCoordinate_)
{
	return Zeros(6,2);
}

RMatrix GJointPlanar::get_DdSDq(GCoordinate *pCoordinate_)
{
	return Zeros(6,2);
}
