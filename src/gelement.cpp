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
#include <string>
#include <sstream>
#include "gelement.h"
#include "gcoordinate.h"
#include "liegroup.h"

//=============================================================
//                 GElement
//=============================================================
std::string GElement::getInfoStr()
{
	std::stringstream sstr;
	sstr << "GElement:: name = " << name << ", id = " << id << std::endl;
	return sstr.str();
}

//=============================================================
//                 GElementWithCoordinates
//=============================================================
void GElementWithCoordinates::clear()
{
	GElement::clear();
	pCoordinates.clear();
}

int GElementWithCoordinates::getIndexOfCoordinate(GCoordinate *pcoord_)
{
	std::list<GCoordinate *>::iterator iter_pcoord, begin = pCoordinates.begin(), end = pCoordinates.end();
	iter_pcoord = std::find(begin, end, pcoord_);
	if ( iter_pcoord == end ) {
		return -1;
	} else {
		return (int)std::distance(begin, iter_pcoord);
	}
}

std::vector<GCoordinate*> GElementWithCoordinates::getPrescribedCoordinates()
{
	int i;
	std::vector<GCoordinate*> pcoords;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++, i++) {
		if ( (*iter_pcoord)->bPrescribed ) {
			pcoords.push_back(*iter_pcoord);
		}
	}
	return pcoords;
}

std::vector<GCoordinate*> GElementWithCoordinates::getUnprescribedCoordinates()
{
	int i;
	std::vector<GCoordinate*> pcoords;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++, i++) {
		if ( !(*iter_pcoord)->bPrescribed ) {
			pcoords.push_back(*iter_pcoord);
		}
	}
	return pcoords;
}

std::vector<int> GElementWithCoordinates::getIndexOfPrescribedCoordinates()
{
	int i;
	std::vector<int> idx;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++, i++) {
		if ( (*iter_pcoord)->bPrescribed ) {
			idx.push_back(i);
		}
	}
	return idx;
}

std::vector<int> GElementWithCoordinates::getIndexOfUnprescribedCoordinates()
{
	int i;
	std::vector<int> idx;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (i=0, iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++, i++) {
		if ( !(*iter_pcoord)->bPrescribed ) {
			idx.push_back(i);
		}
	}
	return idx;
}

void GElementWithCoordinates::initCoordinates()
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->init();
	}
}

void GElementWithCoordinates::set_q(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->q = x_;
	}
}

void GElementWithCoordinates::set_dq(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dq = x_;
	}
}

void GElementWithCoordinates::set_ddq(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->ddq = x_;
	}
}

void GElementWithCoordinates::set_tau(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->tau = x_;
	}
}

void GElementWithCoordinates::set_qLL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->qLL = x_;
	}
}

void GElementWithCoordinates::set_dqLL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dqLL = x_;
	}
}

void GElementWithCoordinates::set_ddqLL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->ddqLL = x_;
	}
}

void GElementWithCoordinates::set_tauLL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->tauLL = x_;
	}
}

void GElementWithCoordinates::set_qUL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->qUL = x_;
	}
}

void GElementWithCoordinates::set_dqUL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dqUL = x_;
	}
}

void GElementWithCoordinates::set_ddqUL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->ddqUL = x_;
	}
}

void GElementWithCoordinates::set_tauUL(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->tauUL = x_;
	}
}

void GElementWithCoordinates::set_aux(gReal x_)
{
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->aux = x_;
	}
}

void GElementWithCoordinates::set_q(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->q = x_[i++];
	}
}

void GElementWithCoordinates::set_dq(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dq = x_[i++];
	}
}

void GElementWithCoordinates::set_ddq(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->ddq = x_[i++];
	}
}

void GElementWithCoordinates::set_tau(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->tau = x_[i++];
	}
}

void GElementWithCoordinates::set_DqDp(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DqDp = x_[i++];
	}
}

void GElementWithCoordinates::set_DdqDp(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DdqDp = x_[i++];
	}
}

void GElementWithCoordinates::set_DddqDp(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DddqDp = x_[i++];
	}
}

void GElementWithCoordinates::set_DtauDp(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->DtauDp = x_[i++];
	}
}

void GElementWithCoordinates::set_qLL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->qLL = x_[i++];
	}
}

void GElementWithCoordinates::set_dqLL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dqLL = x_[i++];
	}
}

void GElementWithCoordinates::set_ddqLL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->ddqLL = x_[i++];
	}
}

void GElementWithCoordinates::set_tauLL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->tauLL = x_[i++];
	}
}

void GElementWithCoordinates::set_qUL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->qUL = x_[i++];
	}
}

void GElementWithCoordinates::set_dqUL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->dqUL = x_[i++];
	}
}

void GElementWithCoordinates::set_ddqUL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->ddqUL = x_[i++];
	}
}

void GElementWithCoordinates::set_tauUL(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->tauUL = x_[i++];
	}
}

void GElementWithCoordinates::set_aux(const gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		(*iter_pcoord)->aux = x_[i++];
	}
}

void GElementWithCoordinates::get_q(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->q;
	}
}

void GElementWithCoordinates::get_dq(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->dq;
	}
}

void GElementWithCoordinates::get_ddq(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->ddq;
	}
}

void GElementWithCoordinates::get_tau(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->tau;
	}
}

void GElementWithCoordinates::get_DqDp(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->DqDp;
	}
}

void GElementWithCoordinates::get_DdqDp(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->DdqDp;
	}
}

void GElementWithCoordinates::get_DddqDp(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->DddqDp;
	}
}

void GElementWithCoordinates::get_DtauDp(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->DtauDp;
	}
}

void GElementWithCoordinates::get_qLL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->qLL;
	}
}

void GElementWithCoordinates::get_dqLL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->dqLL;
	}
}

void GElementWithCoordinates::get_ddqLL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->ddqLL;
	}
}

void GElementWithCoordinates::get_tauLL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->tauLL;
	}
}

void GElementWithCoordinates::get_qUL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->qUL;
	}
}

void GElementWithCoordinates::get_dqUL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->dqUL;
	}
}

void GElementWithCoordinates::get_ddqUL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->ddqUL;
	}
}

void GElementWithCoordinates::get_tauUL(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->tauUL;
	}
}

void GElementWithCoordinates::get_aux(gReal *x_)
{
	int i = 0;
	std::list<GCoordinate *>::iterator iter_pcoord;
	for (iter_pcoord = pCoordinates.begin(); iter_pcoord != pCoordinates.end(); iter_pcoord++) {
		x_[i++] = (*iter_pcoord)->aux;
	}
}


