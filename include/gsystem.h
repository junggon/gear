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
//         GSystem: class for articulated rigid-body systems
// 
//                                                               junggon@gmail.com
//================================================================================

#ifndef _GEAR_SYSTEM_
#define _GEAR_SYSTEM_

#include <list>
#include <vector>
#include <string>
#include "gelement.h"
#include "gcoordinate.h"
#include "gbody.h"
#include "gjoint.h"
#include "gforce.h"
#include "rmatrix3j.h"
#include "liegroup.h"

//=============================================================
//                 GSystem
//=============================================================
class GSystem: public GElementWithCoordinates
{
public:	
	GBody *pGround;							// pointer to the ground body

	std::list<GBody *> pBodies;				// pointer to bodies (pGround is not included in pBodies! See buildSystem().)
	std::list<GJoint *> pJoints;			// pointer to joints
	std::list<GForce *> pForces;			// pointer to force elements (e.g., spring-dampers)

	std::vector<RMatrix> Jacobian;			// Jacobian[i] = Jacobian of pBodies[i]

public:
	GSystem();
	virtual ~GSystem() {}

public:
	// ---------- building system ----------------------------------------------------

	virtual bool buildSystem(GBody *pGround_, bool b_scan_fwd_joint_chain_ = true);	
											// build system by scanning from the ground
											// All the bodies having a direct or indirect network connection to the ground will be scanned.
											// Set b_scan_fwd_joint_chain_ = true to scan forward joint chain for each body.
											
	virtual bool buildSystemWith(GBody *pGround_, std::vector<GBody *> pFirstBodiesIn_, bool b_scan_fwd_joint_chain_ = true);
											// build system by scanning from the child bodies of the ground which are included in pFirstBodiesIn_.
											// Set b_scan_fwd_joint_chain_ = true to scan forward joint chain for each body.
			
	virtual bool buildSystemWithout(GBody *pGround_, std::vector<GBody *> pFirstBodiesOut_, bool b_scan_fwd_joint_chain_ = true);
											// build system by scanning from the child bodies of the ground which are NOT included in pFirstBodiesOut_.
											// Set b_scan_fwd_joint_chain_ = true to scan forward joint chain for each body.


	// ---------- gravity ------------------------------------------------------------

	bool setGravity(Vec3 g_);				// set gravity e.g. setGravity(Vec3(0,0,-9.81))
	Vec3 getGravity();

	// ---------- access to system elements ------------------------------------------

	GBody* getGround() { return pGround; }	// return pointer to the ground

	std::vector<GBody*> getBodies();		// return pointers to the bodies
	std::vector<GJoint*> getJoints();		// return pointers to the joints

	GBody* getBody(std::string name_);		// return pointer to a body whose name is name_ (return NULL if there is no such a body)
	GBody* getBody(int idx_);				// return pointer to a body whose index is idx_ (zero-based index assumed in pBodies) (return NULL if there is no such a joint)

	GJoint* getJoint(std::string name_);	// return pointer to a joint whose name is name_ (return NULL if there is no such a joint)
	GJoint* getJoint(int idx_);				// return pointer to a joint whose index is idx_ (zero-based index assumed in pJoints) (return NULL if there is no such a joint)

	int getIndexBody(std::string name_);	// return index of a body whose name is name_
	int getIndexJoint(std::string name_);	// return index of a joint whose name is name_

	// ---------- basic functions for dynamics simulation -------------------------

	virtual void initBodyForces();			// initializes external body forces

	virtual void initJointTorques();		// initializes joint torques
	
	virtual void initBodyForcesAndJointTorques(); // initializes external body forces and joint torques

	virtual bool stepSimulation(gReal h_);	// step simulation with the step size (h_)
											// (internally, this calculates joint acceleration|torque, integrates the accelerations, and updates position/velocity of the bodies.)
											// ** Simulation inputs must be set before calling this. 
											//    (inputs = joint torque|acceleration + external forces acting on the bodies)
											// ** This is just an example implementation with mixed Euler integration for showing how to write a simulation code. 
											//    (Users may need to write their own implementation with their favorite integrator.)
											// ** A typical implementation of stepSimulation()
											//        0. (NOT INCLUDED HERE) initBodyForcesAndJointTorques()
											//        1. (NOT INCLUDED HERE) apply external body forces (e.g., penality-based contact forces)
											//        2. (NOT INCLUDED HERE) apply joint input (e.g., joint torques or acceleration from a controller)
											//        3. apply force elements to the bodies and joints
											//        4. calculate joint accelerations (for unprescribed coordinates) and torques (for prescribed coordinates)
											//        5. integrate joint accelerations to obtain displacements/velocities of the joints in the next time step
											//        6. update kinematics (update positions and velocities of the bodies)
											// ** A typical procedure for simulation
											//      set initial system state (updateKinematics() must be called at the end of this.)
											//      while (t<tf) {
											//        stepSimulation(h)
											//        t += h
											//        render the system if needed
											//      }

	// ---------- advanced functions for kinematics and dynamics -------------------------

	virtual void updateKinematics();		// update all kinematic information of bodies, such as position, velocity and acceleration, 
											// with current pCoordinates[]->(q,dq,ddq).

	virtual void updateSystemJacobianOfAllBodies();	// update system Jacobian of all bodies
											// prerequisite: updateKinematics() or calcDynamics()

	virtual void calcInverseDynamics();		// calculate pCoordinates[]->tau from current pCoordinates[]->(q,dq,ddq)
											// (all coordinates will be regarded to be prescribed.)

	virtual void calcForwardDynamics();		// calculate pCoordinates[]->ddq from current pCoordinates[]->(q,dq,tau)
											// (all coordinates will be regarded to be unprescribed.)

	virtual void calcDynamics();			// calculate pCoordinates[]->(ddq|tau) from current pCoordinates[]->(q,dq,tau|ddq)
											// if pCoordinates[]->bPrescribed == true, then pCoordinates[]->tau will be updated.
											// if pCoordinates[]->bPrescribed == false, then pCoordinates[]->ddq will be updated.
											// (forward/inverse/hybrid dynamics can be covered by this function only.)

	void getMassMatrix(RMatrix &M);			// M = mass matrix
	void getEquationsOfMotion(RMatrix &M, RMatrix &b); // M * ddq + b = tau
	bool calcProductOfInvMassAndMatrix(RMatrix &invM_A, const RMatrix &A); // calculate Inv(M)*A where M = mass matrix using O(n^2)
	bool calcProductOfInvMassAndMatrix2(RMatrix &invM_A, const RMatrix &A); // calculate Inv(M)*A where M = mass matrix using O(n^2)


	// ---------- derivatives of dynamics ----------------------------------------------

	virtual void diffDynamics();			// calculate pCoordinates[]->(DddqDp|DtauDp)
											// prerequisite: setting differentiating variable p

	virtual void setDeriv_Dq(GCoordinate *pCoordinate_);	// set pCoordinate_->q as the differentiating variable
	virtual void setDeriv_Ddq(GCoordinate *pCoordinate_);	// set pCoordinate_->dq as the differentiating variable
	virtual void setDeriv_Dddq(GCoordinate *pCoordinate_);	// set pCoordinate_->ddq as the differentiating variable
	virtual void setDeriv_Dtau(GCoordinate *pCoordinate_);	// set pCoordinate_->tau as the differentiating variable
	virtual void setDeriv_DFe(GBody *pBody_, int idx_);		// set idx_-th element of the external force acting on pBody_ as the differentiating variable
	

	// ---------- mass, center of mass, momentum --------------------------------------

	gReal getMass();						// return total mass of the system
	
	Vec3 getPositionCOMGlobal();			// return the position of the center of mass of the system w.r.t. {global}

	Vec3 getVelocityCOMGlobal();			// return the velocity of the center of mass of the system w.r.t. {global}

	Vec3 getAccelerationCOMGlobal();		// return the acceleration of the center of mass of the system w.r.t. {global}
	
	dse3 getMomentumGlobal();				// return Hg = the momentum of the system w.r.t. {global}

	dse3 getMomentumCOM();					// return Hc = the momentum of the system w.r.t. {com}
											//      {com} = a moving coordinate frame whose origin is located at the center of mass of the system
											//              and whose orientation is always aligned with that of {global}.

	gReal getGravitationalPotentialEnergy();	// return gravitational potential energy
	gReal getKineticEnergy();					// return kinetic energy

	// ---------- derivatives of center of mass, momentum -----------------------------
	//
	// ** Local information of all joints needs to be updated before calling functions below.  
	// ** updateKinematics() or calcDynamics() will do this automatically.                     
	// ** To do this manually, use update_joint_local_info_short() or update_joint_local_info().

	void calcDerivative_PositionCOMGlobal_Dq(std::vector<GCoordinate *> pCoordinates_, RMatrix &DpDq_);
											// DpDq_ = the derivative of the position of the center of mass w.r.t. pCoordinates_[]->q
    
	void calcDerivative_PositionCOMGlobal_Dq(RMatrix &DpDq_);
											// DpDq_ = the derivative of the position of the center of mass w.r.t. GSystem::pCoordinates[]->q

	void calcDerivative_PositionCOMGlobal_Dq_2(RMatrix &DpDq_);
											// DpDq_ = the derivative of the position of the center of mass w.r.t. GSystem::pCoordinates[]->q
											// prerequisite: updateSystemJacobianOfAllBodies()

	void calcDerivative_MomentumGlobal_Dq_Ddq(std::vector<GCoordinate *> pCoordinates_, RMatrix &DHgDq_, RMatrix &DHgDdq_);
											// DHgDq_, DHgDdq_ = the derivatives of Hg w.r.t. pCoordinates_[]->q and pCoordinates_[]->dq

	void calcDerivative_MomentumGlobal_Dq_Ddq(RMatrix &DHgDq_, RMatrix &DHgDdq_);
											// DHgDq_, DHgDdq_ = the derivatives of Hg w.r.t. GSystem::pCoordinates[]->q and GSystem::pCoordinates[]->dq

	void calcDerivative_MomentumCOM_Dq_Ddq(std::vector<GCoordinate*> pCoordinates_, RMatrix &DHcDq_, RMatrix &DHcDdq_);
											// DHcDq_, DHcDdq_ = the derivatives of Hc w.r.t. pCoordinates_[]->q and pCoordinates_[]->dq

	void calcDerivative_MomentumCOM_Dq_Ddq(RMatrix &DHcDq_, RMatrix &DHcDdq_);
											// DHcDq_, DHcDdq_ = the derivatives of Hc w.r.t. GSystem::pCoordinates[]->q and GSystem::pCoordinates[]->dq

	// ---------- sub-functions ------------------------------------------------------
	
	// sub-functions for buildSystem()
	bool _scanBodiesAndJoints(GBody *pbody_);	// scan bodies and joints
	bool _scanCoordinates();					// scan coordinates for bodies and joints
	bool _findParentChildRelation();			// find pBodies[]->(pBaseJoint, pParentBody, pChildBodies)
	bool _findFwdJointChainOfBodies();			// find the forward joint chain of each body
	bool _scanForceElements();					// scan spring-damper elements connected to the bodies and joints
	bool _getReady();							// get ready

	// sub-functions for updating joint information
	void update_joint_local_info_short();		// updates part of local information of joints (T, inv_T, S)
	void update_joint_local_info();				// updates local information of joints
	void update_joint_global_info();			// updates global location of joints

	// sub-function for updating body Jacobian
	void update_Jacobian_child_bodies(GBody *pbody_, int idx_, const RMatrix &S_);	
												// update idx_-th column of pbody_->pChildBodies[]->Jacobian
												// S_ = idx_-th column of pbody_->Jacobian

	// sub-functions for Newton-Euler inverse dynamics
	void neFwdRecursion_a();
	void neBwdRecursion_b();

	// sub-functions for calcDynamics()
	void fsFwdRecursion_a();
	void fsBwdRecursion_b();
	void fsFwdRecursion_c();

	// sub-functions for dynamics with zero gravity and zero joint velocity
	void initDynamicsWithZeroGravityAndVelocity();	// calculates intermediate quantities depending pose (joint angle) assuming zero joint velocity
	void calcDynamicsWithZeroGravityAndVelocity();	// calculates hybrid dynamics with current joint input (acceleration or torque) assuming zero gravity and zero joint velocity (the result will be saved in pCoordinates[]->(tau|ddq).)
	
	// sub-functions for inverse dynamics with zero gravity and zero joint velocity
	void initInverseDynamicsWithZeroGravityAndVelocity();
	void calcInverseDynamicsWithZeroGravityAndVelocity();

	// sub-functions for diffDynamics()
	void fsFwdRecursion_DaDp();
	void fsBwdRecursion_DbDp();
	void fsFwdRecursion_DcDp();

	// sub-functions for the derivatives of the position of the center of mass and momentum
	Vec3 getDerivative_PositionCOMGlobal_Dq_2(GCoordinate *pCoordinate_);
	Vec3 getDerivative_PositionCOMGlobal_Dq(GCoordinate *pCoordinate_);
											// return the derivative of the position of the center of mass of the system w.r.t. pCoordinate_->q
											// Prerequisites: prerequisites for pBodies[]->getDerivative_PositionCOMGlobal_Dq(pCoordinate_)
											//                ( Usually, it requires that
											//					 1. pBodies[]->fwdJointChain.M1 = identity and pBodies[]->fwdJointChain.jointLoopConstraintType = JOINTLOOP_ORIENTATION_POSITION.
											//					 2. pBodies[]->fwdJointChain.update_J(); )
	dse3 getDerivative_MomentumGlobal_Dq(GCoordinate *pCoordinate_);
	dse3 getDerivative_MomentumGlobal_Ddq(GCoordinate *pCoordinate_);
											// return the derivatives of the momentum w.r.t. pCoordinate_->q and pCoordinate_->dq
											// Prerequisites: prerequisites for pBodies[]->getMomentumGlobalDerivative_[Dq,Ddq](pCoordinate_)
											//                ( Usually, it requires that
											//					 1. pBodies[]->fwdJointChain.M1 = identity and pBodies[]->fwdJointChain.jointLoopConstraintType = JOINTLOOP_ORIENTATION_POSITION.
											//					 2. pBodies[]->fwdJointChain.update_J(); )

	// --------- auxiliary functions ---------------------------------------------------

	std::string getInfoStr();
	void setAllJointsPrescribed(bool b);
};



#endif

