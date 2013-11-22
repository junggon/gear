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
//         GCoordinate: class for coordinates
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_COORDINATE_
#define _GEAR_COORDINATE_

#include "greal.h"

//=============================================================
//                 GCoordinate
//=============================================================
class GCoordinate
{
public:
	gReal q;								// displacement
	gReal dq;								// velocity
	gReal ddq;								// acceleration
	gReal tau;								// torque or force

	gReal DqDp, DdqDp, DddqDp, DtauDp;		// derivatives w.r.t. an arbitrary scalar variable p

	bool bPrescribed;						// set bPrescribed = true if ddq are prescribed
											// ** Do NOT directly set this! (Only use GJoint::setPrescribed(bool b_).)

	gReal qLL, qUL;							// lower and upper limits of q
	gReal dqLL, dqUL;						// lower and upper limits of dq
	gReal ddqLL, ddqUL;						// lower and upper limits of ddq
	gReal tauLL, tauUL;						// lower and upper limits of tau

	gReal aux;								// an auxiliary quantity

public:
	GCoordinate() : q(0), dq(0), ddq(0), tau(0), DqDp(0), DdqDp(0), DddqDp(0), DtauDp(0), bPrescribed(false)
		, qLL((gReal)-1E20), qUL((gReal)1E20), dqLL((gReal)-1E20), dqUL((gReal)1E20), ddqLL((gReal)-1E20), ddqUL((gReal)1E20), tauLL((gReal)-1E20), tauUL((gReal)1E20)
		, aux(0)	{}
	~GCoordinate() {}

	void init() {
		q = dq = ddq = tau = DqDp = DdqDp = DddqDp = DtauDp = (gReal)0.0;
		qLL = dqLL = ddqLL = tauLL = (gReal)-1E20;
		qUL = dqUL = ddqUL = tauUL = (gReal)1E20;
		aux = (gReal)0.0;
	}
};



#endif

