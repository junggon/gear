#ifdef _WIN32
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#endif

#include "myworld.h"

using namespace std;

void glsubBox(double dx, double dy, double dz)
{
	glBegin(GL_QUADS);
		glNormal3d(1, 0, 0); glVertex3d(dx/2, -dy/2, dz/2); glVertex3d(dx/2, -dy/2, -dz/2); glVertex3d(dx/2, dy/2, -dz/2); glVertex3d(dx/2, dy/2, dz/2);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3d(0, 1, 0); glVertex3d(dx/2, dy/2, dz/2); glVertex3d(dx/2, dy/2, -dz/2); glVertex3d(-dx/2, dy/2, -dz/2); glVertex3d(-dx/2, dy/2, dz/2);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3d(-1, 0, 0); glVertex3d(-dx/2, dy/2, -dz/2); glVertex3d(-dx/2, -dy/2, -dz/2); glVertex3d(-dx/2, -dy/2, dz/2); glVertex3d(-dx/2, dy/2, dz/2);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3d(0, -1, 0); glVertex3d(-dx/2, -dy/2, dz/2); glVertex3d(-dx/2, -dy/2, -dz/2); glVertex3d(dx/2, -dy/2, -dz/2); glVertex3d(dx/2, -dy/2, dz/2);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3d(0, 0, 1); glVertex3d(dx/2, -dy/2, dz/2); glVertex3d(dx/2, dy/2, dz/2); glVertex3d(-dx/2, dy/2, dz/2); glVertex3d(-dx/2, -dy/2, dz/2);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3d(0, 0, -1); glVertex3d(dx/2, -dy/2, -dz/2); glVertex3d(-dx/2, -dy/2, -dz/2); glVertex3d(-dx/2, dy/2, -dz/2); glVertex3d(dx/2, dy/2, -dz/2);
	glEnd();
}

//=============================================================
//                 RigidBody
//=============================================================
void RigidBody::render()
{
	glColor3d(_color[0],_color[1],_color[2]);
	glPushMatrix();
	glMultMatrixd(getPoseGlobal().GetArray()); // move to {body}
	glMultMatrixd(_T0.GetArray()); // move to box's reference frame
	glsubBox(_Lx, _Ly, _Lz);
	glPopMatrix();
}

//=============================================================
//                 MyWorld
//=============================================================
bool MyWorld::create()
{
	//----------------------------------------	
	//  link
	//	-----------------------------
	//	|             z             |
	//	|             |             | (A local coordinate frame {link} is attached to each link.)
	//	|             ---> x        |
	//	|                           |
	//	-----------------------------
	// 
	//  Note: There is no restriction on where the local coordinate frame must be located.
	//        In this illustrative example, {link} is attached to the center of mass, i.e., {link} = {com}. 
	//
	//----------------------------------------
	//	two-link model
	//
	//  Z                     
	//  |                     
	//  O=-->X====O==========     (X,Y,Z) represents the world coordinate frame {global}  ({global} is the same as the coordinate frame of ground.)
	// j1  link1  j2  link2
	//
	//----------------------------------------

	double Lx = 0.3, Ly = 0.05, Lz = 0.05; // link dimension
	double M = 1; // link mass
	double Ixx = 1./12.*M*(Ly*Ly+Lz*Lz), Iyy = 1./12.*M*(Lx*Lx+Lz*Lz), Izz = 1./12.*M*(Lx*Lx+Ly*Ly), Ixy = 0, Ixz = 0, Iyz = 0; // moments of inertia with respect to the mass center
	SE3 T_link_com = SE3(); // transformation from {link} to {com} (In this example, T_link_com = identity because {link} = {com}.)

	// set mass and moments of inertia of each link
	link1.setMass(M, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, T_link_com); // Internally, moments of inertia with respect to {link} are calculated using the moments of inertia with respect to {com} and T_link_com.
	link2.setMass(M, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, T_link_com);

	// set box rendering (optional)
	link1.setBoxRendering(Lx, Ly, Lz, SE3(), Vec3(0,0,1));
	link2.setBoxRendering(Lx, Ly, Lz, SE3(), Vec3(1,0,0));

	// connect the links with joints
	//   Let's think about connecting two links (linkA and linkB) with a joint.
	//   A joint has two reference coordinate frames {joint left} and {joint right} which are attached to linkA and linkB respectively,
	//   and {joint left} and {joint right} are aligned with each other when the joint angle or displacement is zero.
	//   usage:
	//     joint.connectBodies(linkA, linkB) connects linkA and linkB with joint.
	//     joint.setPosition(pA, pB) sets the positions of {joint left} and {joint right} with respect to {linkA} and {linkB} respectively.
	//     joint.setOrientation(Ra, Rb) sets the orientations of {joint left} and {joint right} with respect to {linkA} and {linkB} respectively.
	joint1.connectBodies(&ground, &link1);  joint1.setPosition(Vec3(0,0,0), Vec3(-0.5*Lx,0,0)); // (0,0,0) = position of joint1 in {global}, (-0.5*Lx,0,0) = position of join1 in {link1}
	joint2.connectBodies(&link1, &link2);   joint2.setPosition(Vec3(0.5*Lx,0,0), Vec3(-0.5*Lx,0,0));

	// set rotating axes of the revolute joints
	joint1.setAxis(0,1,0); // rotation axis with respect to {joint left}
	joint2.setAxis(0,1,0);

	// build the two-link system
	system.buildSystem(&ground); // GEAR scans the link-joint structure automatically from ground

	// set gravity
	system.setGravity(Vec3(0,0,-9.81));

	return true;
}

void MyWorld::render()
{
	for (std::list<GBody*>::iterator iter_pbody = system.pBodies.begin(); iter_pbody != system.pBodies.end(); iter_pbody++) {
		(*iter_pbody)->render();
	}
}

void MyWorld::testfunc()
{
	RMatrix q = Rand(2,1), dq = 0.1*Rand(2,1), ddq, tau;

	// getting basic information
	cout << "dof = " << system.getNumCoordinates() << endl;
	cout << "mass = " << system.getMass() << endl;

	// testing kinematics
	{
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_ddq(RMatrix("1 2")); // set joint acceleration (unit: rad)
		system.updateKinematics(); // update position, velocity, and acceleration of the system (links and joints)
		cout << "p1cg = " << link1.getPositionCOMGlobal() << endl; // position of link1's center of mass w.r.t. {global}
		cout << "p1g = " << link1.getPositionGlobal() << endl; // position of link1's local coordinate frame origin w.r.t. {global}
		cout << "v2cg = " << link2.getVelocityCOMGlobal() << endl; // velocity of link2's center of mass w.r.t. {global}
		cout << "p2g = " << link2.getPositionGlobal() << endl; // position of link2's local coordinate frame origin w.r.t. {global}
		se3 V2g = link2.getVelocityGlobal(); // generalized velocity (angular velocity, linear velocity) of link2 w.r.t. {global}
		se3 V2 = link2.V; // generalized velocity of link2 w.r.t. {link2} ({link2} = local body coordinate frame of link2)
		system.updateSystemJacobianOfAllBodies(); // update Jacobian matrices of all links using current system position and velocity
		RMatrix J2 = link2.Jacobian; // Jacobian of link2 mapping the system velocity dq to link2's generalized velocity V2 (w.r.t. {link2})
		cout << "Norm(V2-J2*dq) = " << FNorm(RMatrix(6,1,V2.GetArray()) - J2 * system.get_dq()) << endl; // this must be zeros (V2 = J2*dq)
		SO3 R2 = link2.getOrientationGlobal(); // orientation of {link2} w.r.t. {global}
		RMatrix R2_R2 = Zeros(6,6); R2_R2.Push(0,0,RMatrix(3,3,R2.GetArray())); R2_R2.Push(3,3,RMatrix(3,3,R2.GetArray())); // R2_R2 = [R2, 0; 0 R2]
		RMatrix J2g = R2_R2 * J2; // Jacobian of link2 mapping the system velocity dq to link2's generalized velocity V2g (w.r.t. {global})
		cout << "Norm(V2g-J2g*dq) = " << FNorm(RMatrix(6,1,V2g.GetArray()) - J2g * system.get_dq()) << endl; // this must be zeros (V2g = J2g*dq)
	}
	
	// testing hybrid dynamics: (q,dq,ddq_a,tau_b) --> (tau_a,ddq_b) where 'a' = prescribed joints, 'b' = unprescribed joints
	// (inverse and forward dynamics are the two extreme cases of the hybrid dynamics)
	{
		joint1.setPrescribed(false); // joint 1 is 'unprescribed' (i.e., input = torque)
		joint2.setPrescribed(true); // joint 2 is 'prescribed' (i.e., input = acceleration)
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_ddq(RMatrix("1 2")); // set joint acceleration (only acceleration for prescribed joints will be considered in calcDynamics())
		system.set_tau(RMatrix("0.1 0.2")); // set joint torque (only torque for unprescribed joints will be will be considered in calcDynamics())
		system.calcDynamics(); // compute torque for prescribed joints and acceleration for unprescribed joints
		cout << "(ddq1,tau2) = (" << joint1.coordinate.ddq << ", " << joint2.coordinate.tau << ")" << endl; // print joint1's acceleration and joint2's torque
	}

	// testing inverse dynamics: (q,dq,ddq) --> tau
	{
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_ddq(RMatrix("1 2")); // set joint acceleration
		system.calcInverseDynamics(); // compute joint torque using inverse dynamics (All joints are temporarily regarded as "prescribed" joints.)
		cout << "tau = " << system.get_tau() << endl; // print joint torque
	}

	// testing forward dynamics: (q,dq,tau) --> ddq
	{
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_tau(RMatrix("0.1 0.2")); // set joint torque
		system.calcForwardDynamics(); // compute joint acceleration using forward dynamics (All joints are temporarily regarded as "unprescribed" joints.)
		cout << "ddq = " << system.get_ddq() << endl; // print joint acceleration
	}

	// testing derivatives of dynamics
	{
		joint1.setPrescribed(false); // joint 1 is 'unprescribed'
		joint2.setPrescribed(false); // joint 2 is 'unprescribed'
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_tau(RMatrix("0.1 0.2")); // set joint torque
		system.calcDynamics(); // you must call this first before calling diffDynamics()
		system.setDeriv_Dq(&joint1.coordinate); // set q1 as the differentiating variable
		system.diffDynamics(); // compute D(ddq)/D(q1)
		cout << "D(ddq)/D(q1) = " << system.get_DddqDp() << endl;
		system.setDeriv_Ddq(&joint1.coordinate); // set dq1 as the differentiating variable
		system.diffDynamics(); // compute D(ddq)/D(dq1)
		cout << "D(ddq)/D(dq1) = " << system.get_DddqDp() << endl;
		system.setDeriv_Dtau(&joint1.coordinate); // set tau1 as the differentiating variable
		system.diffDynamics(); // compute D(ddq)/D(tau1)
		cout << "D(ddq)/D(tau1) = " << system.get_DddqDp() << endl;
	}

	cout << "comparing analytical and numerical derivatives..." << endl;
	// compare analytical and numerical derivatives: D(ddq)/D(q)
	{
		joint1.setPrescribed(false); // joint 1 is 'unprescribed'
		joint2.setPrescribed(false); // joint 2 is 'unprescribed'
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_tau(RMatrix("0.1 0.2")); // set joint torque
		system.calcDynamics(); // you must call this first before calling diffDynamics()
		vector<GCoordinate*> pcoords = system.getCoordinates();
		RMatrix DddqDq(2,2); // analytical derivative DddqDq = [D(ddq)/D(q1), D(ddq)/D(q2)]
		for (int i=0; i<pcoords.size(); i++) {
			system.setDeriv_Dq(pcoords[i]);
			system.diffDynamics();
			DddqDq.Push(0,i,system.get_DddqDp());
		}
		RMatrix NDddqDq(2,2); // numerical derivative NDddqDq = [D(ddq)/D(q1), D(ddq)/D(q2)]
		RMatrix ddq0 = system.get_ddq();
		double h = 1E-6;
		for (int i=0; i<pcoords.size(); i++) {
			pcoords[i]->q += h;
			system.calcDynamics();
			NDddqDq.Push(0,i,(1./h)*(system.get_ddq()-ddq0));
			pcoords[i]->q -= h;
		}
		//cout << "DddqDq = " << DddqDq << "NDddqDq = " << NDddqDq << endl;
		cout << "Norm(DddqDq-NDddqDq) = " << FNorm(DddqDq-NDddqDq) << endl;
	}

	// compare analytical and numerical derivatives: D(ddq)/D(dq)
	{
		joint1.setPrescribed(false); // joint 1 is 'unprescribed'
		joint2.setPrescribed(false); // joint 2 is 'unprescribed'
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_tau(RMatrix("0.1 0.2")); // set joint torque
		system.calcDynamics(); // you must call this first before calling diffDynamics()
		vector<GCoordinate*> pcoords = system.getCoordinates();
		RMatrix DddqDdq(2,2); // analytical derivative DddqDdq = [D(ddq)/D(dq1), D(ddq)/D(dq2)]
		for (int i=0; i<pcoords.size(); i++) {
			system.setDeriv_Ddq(pcoords[i]);
			system.diffDynamics();
			DddqDdq.Push(0,i,system.get_DddqDp());
		}
		RMatrix NDddqDdq(2,2); // numerical derivative NDddqDdq = [D(ddq)/D(dq1), D(ddq)/D(dq2)]
		RMatrix ddq0 = system.get_ddq();
		double h = 1E-6;
		for (int i=0; i<pcoords.size(); i++) {
			pcoords[i]->dq += h;
			system.calcDynamics();
			NDddqDdq.Push(0,i,(1./h)*(system.get_ddq()-ddq0));
			pcoords[i]->dq -= h;
		}
		//cout << "DddqDdq = " << DddqDdq << "NDddqDdq = " << NDddqDdq << endl;
		cout << "Norm(DddqDdq-NDddqDdq) = " << FNorm(DddqDdq-NDddqDdq) << endl;
	}

	// compare analytical and numerical derivatives: D(ddq)/D(tau)
	{
		joint1.setPrescribed(false); // joint 1 is 'unprescribed'
		joint2.setPrescribed(false); // joint 2 is 'unprescribed'
		system.set_q(q); // set joint position
		system.set_dq(dq); // set joint velocity
		system.set_tau(RMatrix("0.1 0.2")); // set joint torque
		system.calcDynamics(); // you must call this first before calling diffDynamics()
		vector<GCoordinate*> pcoords = system.getCoordinates();
		RMatrix DddqDtau(2,2); // analytical derivative DddqDtau = [D(ddq)/D(tau1), D(ddq)/D(tau2)]
		for (int i=0; i<pcoords.size(); i++) {
			system.setDeriv_Dtau(pcoords[i]);
			system.diffDynamics();
			DddqDtau.Push(0,i,system.get_DddqDp());
		}
		RMatrix NDddqDtau(2,2); // numerical derivative NDddqDtau = [D(ddq)/D(tau1), D(ddq)/D(tau2)]
		RMatrix ddq0 = system.get_ddq();
		double h = 1E-6;
		for (int i=0; i<pcoords.size(); i++) {
			pcoords[i]->tau += h;
			system.calcDynamics();
			NDddqDtau.Push(0,i,(1./h)*(system.get_ddq()-ddq0));
			pcoords[i]->tau -= h;
		}
		//cout << "DddqDtau = " << DddqDtau << "NDddqDtau = " << NDddqDtau << endl;
		cout << "Norm(DddqDtau-NDddqDtau) = " << FNorm(DddqDtau-NDddqDtau) << endl;
	}
}

