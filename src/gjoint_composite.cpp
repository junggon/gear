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
#include <algorithm>
#include "gjoint_composite.h"
#include "gjoint.h"
#include "liegroup.h"
#include "rmatrix3j.h"
#include "liegroup_rmatrix3_ext.h"



//=============================================================
//                 GJointComposite
//=============================================================
GJointComposite::GJointComposite()
{
	jointType = GJOINT_COMPOSITE;
	pJoint1 = pJoint2 = NULL;
}

GJointComposite::GJointComposite(GJoint *pjoint1_, GJoint *pjoint2_)
{
	if ( pjoint1_ == NULL || pjoint2_ == NULL ) {
		pJoint1 = pJoint2 = NULL;
	} else {
		compose(pjoint1_, pjoint2_);
	}
}

bool GJointComposite::compose(GJoint *pjoint1_, GJoint *pjoint2_)
{
	std::list<GCoordinate*>::iterator iter_pcoord;

	if ( pjoint1_ == NULL || pjoint2_ == NULL ) return false;

	pJoint1 = pjoint1_;
	pJoint2 = pjoint2_;
	
	pCoordinates.clear();
	for (iter_pcoord = pJoint1->pCoordinates.begin(); iter_pcoord != pJoint1->pCoordinates.end(); iter_pcoord++) {
		pCoordinates.push_back(*iter_pcoord);
	}
	for (iter_pcoord = pJoint2->pCoordinates.begin(); iter_pcoord != pJoint2->pCoordinates.end(); iter_pcoord++) {
		pCoordinates.push_back(*iter_pcoord);
	}

	allocate_memory(getDOF());

	return true;
}

void GJointComposite::update_short()
{
	pJoint1->update_short();
	pJoint2->update_short();

	SE3 inv_T2 = pJoint2->inv_T;

	T = pJoint1->T;
	T *= pJoint2->T;

	inv_T.SetInvOf(T);

	// S = [Ad(inv_T2, pJoint1->S), pJoint2->S]
	S.SetZero(6, getDOF());
	matSet_Ad(S.GetPtr(), inv_T2, pJoint1->S.GetPtr(), pJoint1->getDOF()); //S.Push(0, 0, Ad(inv_T2, pJoint1->S));
	matSet(&S[6*pJoint1->getDOF()], pJoint2->S.GetPtr(), 6*pJoint2->getDOF()); //S.Push(0, pJoint1->getDOF(), pJoint2->S);

	if ( bReversed ) { _update_short_for_reversed_joint(); }
}

void GJointComposite::update()
{
	pJoint1->update();
	pJoint2->update();

	SE3 inv_T2 = pJoint2->inv_T;

	T = pJoint1->T;
	T *= pJoint2->T;

	inv_T.SetInvOf(T);

	Sdq.set_Ad(inv_T2, pJoint1->Sdq);
	Sdq += pJoint2->Sdq;

	dSdq.set_ad(-pJoint2->Sdq, Ad(inv_T2, pJoint1->Sdq));
	dSdq += Ad(inv_T2, pJoint1->dSdq);
	dSdq += pJoint2->dSdq;

	Sddq.set_Ad(inv_T2, pJoint1->Sddq);
	Sddq += pJoint2->Sddq;

	DSdqDt = Sddq;
	DSdqDt += dSdq;

	// S = [Ad(inv_T2, pJoint1->S), pJoint2->S]
	S.SetZero(6, getDOF());
	matSet_Ad(S.GetPtr(), inv_T2, pJoint1->S.GetPtr(), pJoint1->getDOF()); // S.Push(0, 0, Ad(inv_T2, pJoint1->S));
	matSet(&S[6*pJoint1->getDOF()], pJoint2->S.GetPtr(), 6*pJoint2->getDOF()); // S.Push(0, pJoint1->getDOF(), pJoint2->S);

	dS.SetZero(6, getDOF());
	dS.Push(0, 0, - ad(pJoint2->Sdq, Ad(inv_T2, pJoint1->S)) + Ad(inv_T2, pJoint1->dS));
	dS.Push(0, pJoint1->getDOF(), pJoint2->dS);

	if ( bReversed ) { _update_for_reversed_joint(); }
}

RMatrix GJointComposite::get_DSDq(GCoordinate *pCoordinate_)
{
	if ( find(pCoordinates.begin(), pCoordinates.end(), pCoordinate_) == pCoordinates.end() ) return Zeros(6, getDOF());

	RMatrix DSDq, DS1Dq, DS2Dq;

	DS1Dq = -ad(pJoint2->get_S(pCoordinate_), Ad(pJoint2->inv_T, pJoint1->S));
	DS1Dq += Ad(pJoint2->inv_T, pJoint1->get_DSDq(pCoordinate_));

	DS2Dq = pJoint2->get_DSDq(pCoordinate_);

	DSDq.SetZero(6, getDOF());
	DSDq.Push(0, 0, DS1Dq);
	DSDq.Push(0, pJoint1->getDOF(), DS2Dq);
	
	if ( bReversed ) {
		DSDq = -Ad(inv_T, DSDq);
		DSDq -= ad(get_S(pCoordinate_), S);
	}

	return DSDq;
}

RMatrix GJointComposite::get_DdSDq(GCoordinate *pCoordinate_)
{
	if ( find(pCoordinates.begin(), pCoordinates.end(), pCoordinate_) == pCoordinates.end() ) return Zeros(6, getDOF());

	RMatrix DdSDq, DdS1Dq, DdS2Dq, Ss;
	RMatrix dq2(pJoint2->getDOF(),1); pJoint2->get_dq(dq2.GetPtr());

	DdS1Dq = -ad(pJoint2->get_DSDq(pCoordinate_) * dq2, Ad(pJoint2->inv_T, pJoint1->S));
	DdS1Dq += ad(pJoint2->Sdq, ad(pJoint2->get_S(pCoordinate_), Ad(pJoint2->inv_T, pJoint1->S)));
	DdS1Dq -= ad(pJoint2->Sdq, Ad(pJoint2->inv_T, pJoint1->get_DSDq(pCoordinate_)));
	DdS1Dq -= ad(pJoint2->get_S(pCoordinate_), Ad(pJoint2->inv_T, pJoint1->dS));
	DdS1Dq += Ad(pJoint2->inv_T, pJoint1->get_DdSDq(pCoordinate_));

	DdS2Dq = pJoint2->get_DdSDq(pCoordinate_);

	DdSDq.SetZero(6, getDOF());
	DdSDq.Push(0, 0, DdS1Dq);
	DdSDq.Push(0, pJoint1->getDOF(), DdS2Dq);

	if ( bReversed ) {
		RMatrix DSDq = get_DSDq(pCoordinate_);
		RMatrix dq(getDOF(),1); get_dq(dq.GetPtr());
		DdSDq = -Ad(inv_T, DdSDq);
		DdSDq -= ad(get_S(pCoordinate_), dS + ad(Sdq, S));
		DdSDq -= ad(DSDq*dq, S);
		DdSDq -= ad(Sdq, DSDq);
	}

	return DdSDq;
}

void GJointComposite::_update_short_for_reversed_joint()
{
	SE3 T_tmp = T;
	T = inv_T;
	inv_T = T_tmp;

	S = -Ad(inv_T, S);
}

void GJointComposite::_update_for_reversed_joint()
{
	SE3 T_tmp = T;
	T = inv_T;
	inv_T = T_tmp;

	Sdq = -Ad(inv_T, Sdq);
	dSdq = -Ad(inv_T, dSdq);
	Sddq = -Ad(inv_T, Sddq);
	DSdqDt = Sddq + dSdq;
	S = -Ad(inv_T, S);
	dS = -Ad(inv_T, dS) - ad(Sdq, S);
}
