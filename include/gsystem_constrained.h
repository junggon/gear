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
//         GSystemConstrained: class for constrained articulated rigid-body systems
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_SYSTEM_CONSTRAINED_
#define _GEAR_SYSTEM_CONSTRAINED_

#include <list>
#include <vector>
#include "rmatrix3j.h"
#include "liegroup.h"
#include "gsystem.h"
#include "gconstraint.h"
#include "gconstraint_jointloop.h"



//=============================================================
//                 GSystemConstrained
//=============================================================
class GSystemConstrained: public GSystem
{
public:	
	std::list<GConstraint *> pConstraints;	// pointer to constraints

	// setting for constraint solving
	gReal tolerance_C;		// Newton-Raphson iteration stops if fabs(C[i]) < tolerance_C for all i where C[i] means the i-th constraint function value
	gReal tolerance_Cmax;
	gReal tolerance_delq;	// perturbs joint displacement if Norm(del_q) < tolerance_delq and fabs(C[i]) > tolerance_Cmax 
	gReal tolerance_pinv;	// tolerance for matrix pseudo-inverse when solving A*x=b
	gReal perturb_amount;	// amount of joint displacement perturbation
	int maxIterNum;			// Newton-Raphson iteration stops if the iteration number exceeds maxIterNum

	// automatic handling of closed loops
	std::list<GConstraintJointLoop> closedJointLoopConstraints;	// embedded closed joint-loop constraints
	std::list<GBody *> pCutBodies;								// pointers to the dummy bodies attached to the cut-joints (The dummy bodies are automatically created and must be deleted in the destructor.)

	// current constraint forces (result of calcDynamics(), calcForwardDynamics() and calcInverseDynamics())
	RMatrix _lambda;

public:
	GSystemConstrained() : tolerance_C(1E-9), tolerance_Cmax(1E-3), tolerance_delq(1E-9), tolerance_pinv(1E-9), perturb_amount(1E-3), maxIterNum(50) {}
	~GSystemConstrained() { _deleteMemoryForAllCutBodies(); }

	// virtual function override
	bool buildSystem(GBody *pGround_, bool b_scan_fwd_joint_chain_ = true);
	bool stepSimulation(gReal h_);
	void calcDynamics();
	void calcInverseDynamics();	// computes minimum joint torques with assuming all coordinates are prescribed (Note: calcDynamics() chooses a solution with minimum constraint forces when there exist multiple (infinite number of) solutions.)
	void diffDynamics() { std::cerr << "error:: GSystemConstrained::diffDynamics() not implemented yet!" << std::endl; }

	// number of constraints
	int getNumConstraints() { int nc=0; for (std::list<GConstraint*>::iterator iter = pConstraints.begin(); iter != pConstraints.end(); iter++) { nc += (*iter)->getNumConstraints(); } return nc; }

	// add/remove constraints
	bool addConstraint(GConstraint *pConstraint_);
	bool removeConstraint(GConstraint *pConstraint_);
	void removeAllConstraints();

	// get constraint function (C), Jacobian (J), and the time derivative of the Jacobian (dJdt) where C=0 is the system constraints, dC/dt = J * dq where dq = pCoordinates[]->dq, and dJdt = dJ/dt
 	void getConstraintFunction(RMatrix &C);					
	void getConstraintJacobian(RMatrix &J) { return getConstraintJacobian(J, pCoordinates);	}
	void getConstraintJacobianDerivative(RMatrix &dJdt) { return getConstraintJacobianDerivative(dJdt, pCoordinates); }

	// get Jacobian (J) and the time derivative of the Jacobian (dJdt) whose columns are corresponding to the given coordinates pcoords
	void getConstraintJacobian(RMatrix &J, std::vector<GCoordinate*> pcoords);	
	void getConstraintJacobian(RMatrix &J, std::list<GCoordinate*> pcoords) { getConstraintJacobian(J, std::vector<GCoordinate*>(pcoords.begin(),pcoords.end())); }
	void getConstraintJacobianDerivative(RMatrix &dJdt, std::vector<GCoordinate*> pcoords);
	void getConstraintJacobianDerivative(RMatrix &dJdt, std::list<GCoordinate*> pcoords) { getConstraintJacobianDerivative(dJdt, std::vector<GCoordinate*>(pcoords.begin(),pcoords.end())); }

	// get/set parameters for constraint solving
	gReal getTolerance_C() { return tolerance_C; }					// get torlerance for solving the constraints
	gReal getTolerance_Cmax() { return tolerance_Cmax; }			// get torlerance for solving the constraints
	gReal getTolerance_delq() { return tolerance_delq; }			// get torlerance for solving the constraints
	gReal getTolerance_pinv() { return tolerance_pinv; }			// get torlerance for solving the constraints
	int getMaxIterNum() { return maxIterNum; }						// get maximum iteration number for solving the constraints
	void setTolerance_C(gReal tol) { tolerance_C = tol; }			// set tolerance for solving the constraints
	void setTolerance_Cmax(gReal tol) { tolerance_Cmax = tol; }		// set tolerance for solving the constraints
	void setTolerance_delq(gReal tol) { tolerance_delq = tol; }		// set tolerance for solving the constraints
	void setTolerance_pinv(gReal tol) { tolerance_pinv = tol; }		// set tolerance for solving the constraints
	void setMaxIterNum(int n) { maxIterNum = n; }					// set maximum iteration number for solving the constraints

	// check if the constraints are satisfied
	bool checkConstrainedKinematics(gReal eps);	// return true if the current pCoordinates[]->(q,dq,ddq) satisfy the kinematic constraints (C=0, J*dq=0, J*ddq+dJdt*dq=0)
	bool checkConstrainedDynamics(gReal eps);	// return true if the current pCoordinates[]->(q,dq,ddq|tau) and the constraint forces (_lambda) satisfy the constrained dynamics equations: 
												//    M*ddq+b=tau+~J*_lambda, J*ddq+dJdt*dq=0 where M = M(q), b = b(q,dq), J = J(q), dJdt = dJ/dt.
												//    (Note: Call this function right after calling calcDynamics().)

	// enforce the constraints by adjusting displacement(q), velocity(dq) and acceleration(ddq)
	bool enforceConstraints() { return enforceConstraints(getUnprescribedCoordinates()); } // enforce the constraints by adjusting unprescribed coordinates
	bool enforceConstraints(std::vector<GCoordinate*> pcoords);	// enforce the constraints by adjusting pcoords[]->(q,dq,ddq)
	bool enforceConstraints(std::list<GCoordinate*> pcoords) { return enforceConstraints(std::vector<GCoordinate*>(pcoords.begin(),pcoords.end())); }

	//// check constraint solving
	//void checkConstraintSolving() { return checkConstraintSolving(getUnprescribedCoordinates()); }
	//void checkConstraintSolving(std::vector<GCoordinate*> pcoords);	// enforce the constraints by adjusting pcoords[]->(q,dq,ddq)
	//void checkConstraintSolving(std::list<GCoordinate*> pcoords) { return checkConstraintSolving(std::vector<GCoordinate*>(pcoords.begin(),pcoords.end())); }

public:	
	// sub-functions for buildSystem()
	bool _findClosedJointLoopConstraints();
	bool _findClosedJointLoop(GJoint *pCutJoint_, std::list<GJoint *> &loopjoints_);
	bool _findJointLoop(GJoint *pEndJoint_, std::list<GJoint *> &loopJoints_);
	GBody* _createCutBody(GJoint *pcutjoint);
	void _deleteMemoryForAllCutBodies() { for (std::list<GBody*>::iterator iter = pCutBodies.begin(); iter != pCutBodies.end(); iter++) { delete *iter; } }

	// sub-functions for adjusting unprescribed coordinates to satisfy the constraints
	bool _adjust_qv(std::vector<GCoordinate*> pcoords);	// adjust qv = pcoords[]->q to satisfy the constraints
	bool _adjust_qv(std::list<GCoordinate*> pcoords) { return _adjust_qv(std::vector<GCoordinate*>(pcoords.begin(),pcoords.end())); }			
	bool _adjust_dqv(std::vector<GCoordinate*> pcoords, bool bupdate_J = true); // adjust dqv = pcoords[]->dq to satisfy the constraints
	bool _adjust_dqv(std::list<GCoordinate*> pcoords, bool bupdate_J = true) { return _adjust_dqv(std::vector<GCoordinate*>(pcoords.begin(),pcoords.end()), bupdate_J); }			
	bool _adjust_ddqv(std::vector<GCoordinate*> pcoords, bool bupdate_J = true, bool bupdate_dJdt = true); // adjust ddqv = pcoords[]->ddq to satisfy the constraints
	bool _adjust_ddqv(std::list<GCoordinate*> pcoords, bool bupdate_J = true, bool bupdate_dJdt = true) { return _adjust_ddqv(std::vector<GCoordinate*>(pcoords.begin(),pcoords.end()), bupdate_J, bupdate_dJdt); }	

	// sub-functions for evaluating constraint function, Jacobian, and the time derivative of the Jacobian
	void _update_C();													// update the function values of all constraints
	void _update_J();													// update the Jacobian matrices of all constraints
	void _update_dJdt();												// update the time derivative of the Jacobian matrices of all constraints
	void _copy_C(RMatrix &C);											// copy the function values of each constraint to build the system constraint function C
	void _copy_J(RMatrix &J, std::vector<GCoordinate*> pcoords);		// copy Jacobian of each constraint to build a system constraint Jacobian matrix J whose columns are corresponding to the given coordinates pcoords
	void _copy_J(RMatrix &J, std::list<GCoordinate*> pcoords) { _copy_J(J, std::vector<GCoordinate*>(pcoords.begin(),pcoords.end())); }
	void _copy_dJdt(RMatrix &dJdt, std::vector<GCoordinate*> pcoords);	// copy the time derivative of Jacobian of each constraint to build a system constraint Jacobian derivative dJdt whose columns are corresponding to the given coordinates	
	void _copy_dJdt(RMatrix &dJdt, std::list<GCoordinate*> pcoords) { _copy_dJdt(dJdt, std::vector<GCoordinate*>(pcoords.begin(),pcoords.end())); }

	//// sub-functions for dynamics
	//// inefficient!! we don't have to consider coordinates not related to the constraints.
	//bool _calculateEquivalentIndependentCoordinatesForce(std::list<gReal> &tauu);	
};



#endif

