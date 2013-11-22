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
//         GElement: base class for GEAR elements
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_ELEMENT_
#define _GEAR_ELEMENT_

#include <string>
#include <list>
#include <vector>
#include "rmatrix3j.h"


class GCoordinate;

//=============================================================
//                 GElement
//=============================================================
class GElement
{
protected:
	std::string name;						// name
	int id;									// id

public:
	GElement() : name(""), id(0) {}
	virtual ~GElement() {}

	void setName(std::string name_) { name = name_; }
	void setID(int id_) { id = id_; }

	std::string getName() { return name; }
	int getID() { return id; }

	virtual void clear() { id=0; name=""; } // clear all allocated memories (such as vector and list) and initialize variables
	virtual bool getReady() { return true; }
	virtual std::string getInfoStr();
	virtual void render() {}
};

//=============================================================
//                 GElementWithCoordinates
//=============================================================
class GElementWithCoordinates: public GElement
{
public:
	std::list<GCoordinate *> pCoordinates;	// pointer to coordinates
	
public:
	GElementWithCoordinates() {}
	~GElementWithCoordinates() {}

	virtual void clear();

	int getNumCoordinates() { return (int)pCoordinates.size(); }
	std::vector<GCoordinate*> getCoordinates() { return std::vector<GCoordinate*>(pCoordinates.begin(), pCoordinates.end()); }

	int getIndexOfCoordinate(GCoordinate *pcoord_);	// return the index of the coordinate (return -1 if not exists)

	std::vector<GCoordinate*> getPrescribedCoordinates();
	std::vector<GCoordinate*> getUnprescribedCoordinates();

	std::vector<int> getIndexOfPrescribedCoordinates();
	std::vector<int> getIndexOfUnprescribedCoordinates();

	void initCoordinates();

	// set coordinate values with a scalar value
	void set_q(const gReal x_);
	void set_dq(const gReal x_);
	void set_ddq(const gReal x_);
	void set_tau(const gReal x_);

	void set_qLL(const gReal x_);
	void set_dqLL(const gReal x_);
	void set_ddqLL(const gReal x_);
	void set_tauLL(const gReal x_);

	void set_qUL(const gReal x_);
	void set_dqUL(const gReal x_);
	void set_ddqUL(const gReal x_);
	void set_tauUL(const gReal x_);

	void set_aux(const gReal x_);

	// set coordinate values with an array
	void set_q(const gReal *x_);
	void set_dq(const gReal *x_);
	void set_ddq(const gReal *x_);
	void set_tau(const gReal *x_);

	void set_DqDp(const gReal *x_);
	void set_DdqDp(const gReal *x_);
	void set_DddqDp(const gReal *x_);
	void set_DtauDp(const gReal *x_);

	void set_qLL(const gReal *x_);
	void set_dqLL(const gReal *x_);
	void set_ddqLL(const gReal *x_);
	void set_tauLL(const gReal *x_);

	void set_qUL(const gReal *x_);
	void set_dqUL(const gReal *x_);
	void set_ddqUL(const gReal *x_);
	void set_tauUL(const gReal *x_);

	void set_aux(const gReal *x_);

	// get coordinate values to an array
	void get_q(gReal *x_);
	void get_dq(gReal *x_);
	void get_ddq(gReal *x_);
	void get_tau(gReal *x_);

	void get_DqDp(gReal *x_);
	void get_DdqDp(gReal *x_);
	void get_DddqDp(gReal *x_);
	void get_DtauDp(gReal *x_);

	void get_qLL(gReal *x_);
	void get_dqLL(gReal *x_);
	void get_ddqLL(gReal *x_);
	void get_tauLL(gReal *x_);

	void get_qUL(gReal *x_);
	void get_dqUL(gReal *x_);
	void get_ddqUL(gReal *x_);
	void get_tauUL(gReal *x_);

	void get_aux(gReal *x_);

	// set coordinate values with RMatrix
	bool set_q(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_q(in_.GetPtr()); return true; }
	bool set_dq(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_dq(in_.GetPtr()); return true; }
	bool set_ddq(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_ddq(in_.GetPtr()); return true; }
	bool set_tau(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_tau(in_.GetPtr()); return true; }

	bool set_DqDp(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_DqDp(in_.GetPtr()); return true; }
	bool set_DdqDp(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_DdqDp(in_.GetPtr()); return true; }
	bool set_DddqDp(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_DddqDp(in_.GetPtr()); return true; }
	bool set_DtauDp(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_DtauDp(in_.GetPtr()); return true; }

	bool set_qLL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_qLL(in_.GetPtr()); return true; }
	bool set_dqLL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_dqLL(in_.GetPtr()); return true; }
	bool set_ddqLL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_ddqLL(in_.GetPtr()); return true; }
	bool set_tauLL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_tauLL(in_.GetPtr()); return true; }

	bool set_qUL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_qUL(in_.GetPtr()); return true; }
	bool set_dqUL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_dqUL(in_.GetPtr()); return true; }
	bool set_ddqUL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_ddqUL(in_.GetPtr()); return true; }
	bool set_tauUL(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_tauUL(in_.GetPtr()); return true; }

	bool set_aux(const RMatrix &in_) { if ( in_.RowSize() * in_.ColSize() != getNumCoordinates() ) return false; set_aux(in_.GetPtr()); return true; }

	// get coordinate values to RMatrix
	RMatrix get_q() { RMatrix re(getNumCoordinates(), 1); get_q(re.GetPtr()); return re; }
	RMatrix get_dq() { RMatrix re(getNumCoordinates(), 1); get_dq(re.GetPtr()); return re; }
	RMatrix get_ddq() { RMatrix re(getNumCoordinates(), 1); get_ddq(re.GetPtr()); return re; }
	RMatrix get_tau() { RMatrix re(getNumCoordinates(), 1); get_tau(re.GetPtr()); return re; }

	RMatrix get_DqDp() { RMatrix re(getNumCoordinates(), 1); get_DqDp(re.GetPtr()); return re; }
	RMatrix get_DdqDp() { RMatrix re(getNumCoordinates(), 1); get_DdqDp(re.GetPtr()); return re; }
	RMatrix get_DddqDp() { RMatrix re(getNumCoordinates(), 1); get_DddqDp(re.GetPtr()); return re; }
	RMatrix get_DtauDp() { RMatrix re(getNumCoordinates(), 1); get_DtauDp(re.GetPtr()); return re; }

	RMatrix get_qLL() { RMatrix re(getNumCoordinates(), 1); get_qLL(re.GetPtr()); return re; }
	RMatrix get_dqLL() { RMatrix re(getNumCoordinates(), 1); get_dqLL(re.GetPtr()); return re; }
	RMatrix get_ddqLL() { RMatrix re(getNumCoordinates(), 1); get_ddqLL(re.GetPtr()); return re; }
	RMatrix get_tauLL() { RMatrix re(getNumCoordinates(), 1); get_tauLL(re.GetPtr()); return re; }

	RMatrix get_qUL() { RMatrix re(getNumCoordinates(), 1); get_qUL(re.GetPtr()); return re; }
	RMatrix get_dqUL() { RMatrix re(getNumCoordinates(), 1); get_dqUL(re.GetPtr()); return re; }
	RMatrix get_ddqUL() { RMatrix re(getNumCoordinates(), 1); get_ddqUL(re.GetPtr()); return re; }
	RMatrix get_tauUL() { RMatrix re(getNumCoordinates(), 1); get_tauUL(re.GetPtr()); return re; }

	RMatrix get_aux() { RMatrix re(getNumCoordinates(), 1); get_aux(re.GetPtr()); return re; }
};



#endif

