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

	double Lx = 0.3, Ly = 0.05, Lz = 0.05; // link dimension
	double M = 1; // link mass
	double Ixx = 1./12.*M*(Ly*Ly+Lz*Lz), Iyy = 1./12.*M*(Lx*Lx+Lz*Lz), Izz = 1./12.*M*(Lx*Lx+Ly*Ly), Ixy = 0, Ixz = 0, Iyz = 0; // moments of inertia with respect to the mass center
	SE3 T_link_com = SE3(); // transformation from {link} to {center of mass} (In this example, {link} == {center of mass}, so T_link_com = identity.)

	// set mass properties of the links
	link1.setMass(M, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, T_link_com); // set mass and the moments of inertia of the link (with respect to {link1})
	link2.setMass(M, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, T_link_com);
	link3.setMass(M, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, T_link_com);
	link4.setMass(M, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, T_link_com);

	// set box rendering (optional)
	link1.setBoxRendering(0.4*Lx, Ly, Lz, SE3(), Vec3(0,0,1));
	link2.setBoxRendering(Lx, Ly, Lz, SE3(), Vec3(1,0,0));
	link3.setBoxRendering(Lx, Ly, Lz, SE3(), Vec3(0,1,0));
	link4.setBoxRendering(Lx, Ly, Lz, SE3(), Vec3(1,1,0));

	// set link name
	link1.setName("link1");
	link2.setName("link2");
	link3.setName("link3");
	link4.setName("link4");

	// connect the links with joints
	joint1.connectBodies(&ground, &link1);  joint1.setPosition(Vec3(0,0,0), Vec3(-0.2*Lx,0,0));
	joint2.connectBodies(&link1, &link2);   joint2.setPosition(Vec3(0.2*Lx,0,0), Vec3(-0.5*Lx,0,0));
	joint3.connectBodies(&link2, &link3);   joint3.setPosition(Vec3(0.5*Lx,0,0), Vec3(-0.5*Lx,0,0));
	joint4.connectBodies(&link3, &link4);   joint4.setPosition(Vec3(0.5*Lx,0,0), Vec3(-0.5*Lx,0,0));
	joint5.connectBodies(&link4, &ground);  joint5.setPosition(Vec3(0,0,0), Vec3(2*Lx,0,0.5*Lx)); //joint5.setOrientation(Exp(0,0,0),Exp(0,45*3.14159/180.,0));

	// set rotating axes of the revolute joints
	joint1.setAxis(0,1,0);
	joint2.setAxis(0,1,0);
	joint3.setAxis(0,1,0);
	joint4.setAxis(0,1,0);
	joint5.setAxis(1,0,0);

	// set prescribed/unprescribed joints
	joint1.setPrescribed(true);
	joint2.setPrescribed(false);
	joint3.setPrescribed(false);
	joint4.setPrescribed(false);
	joint5.setPrescribed(false);

	// set joint name
	joint1.setName("joint1"); 
	joint2.setName("joint2"); 
	joint3.setName("joint3");
	joint4.setName("joint4");
	joint5.setName("joint5");

	// build system
	if ( !system.buildSystem(&ground) ) return false;
	if ( !system.enforceConstraints() ) {
		cout << "failed in enforcing constraints!" << endl;
	}
	system.updateKinematics();

	// set gravity
	system.setGravity(Vec3(0,0,-9.81));

	// init
	joint1.set_q(0.0);
	joint1.set_dq(2.0);
	joint1.set_ddq(0.0);
	system.enforceConstraints();
	system.updateKinematics();

	//ofstream fout("dd.txt");
	//fout << system.getInfoStr() << endl;
	//fout.close();

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
}

