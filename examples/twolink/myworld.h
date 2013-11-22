#ifndef _MYWORLD_
#define _MYWORLD_

#include <fstream>
#include <vector>
#include "gear.h"

//=============================================================
//                 RigidBody
//=============================================================
class RigidBody: public GBody
{
public:
	RigidBody() : _Lx(0), _Ly(0), _Lz(0), _T0(SE3()), _color(0,0,1) {}
	~RigidBody() {}
	
	// set box rendering: 
	void setBoxRendering(double Lx, double Ly, double Lz, SE3 T0, Vec3 color=Vec3(0,0,1)) { _Lx = Lx; _Ly = Ly; _Lz = Lz; _T0 = T0; _color = color; }

	// virtual function
	void render();
	
public:
	double _Lx, _Ly, _Lz; 		// box size
	SE3 _T0;					// position and orientation of the box w.r.t. {body} ({body} means a local coordinate frame attached to the body)
	Vec3 _color;				// box color
};


//=============================================================
//                 MyWorld
//=============================================================
class MyWorld
{
public:
	MyWorld() {}
	~MyWorld() {}

	// build a two-link system
	bool create();

	// draw the system
	void render();

	// test function
	void testfunc();

public:

	// ground and links
	RigidBody ground, link1, link2;

	// joints
	GJointRevolute joint1, joint2;
	
	// system (for kinematics and dynamics)
	GSystem system;

	// simulation data
	struct SimulData {
		double t;
		RMatrix q, dq;
	};

	// simulated data
	std::vector<SimulData> simuldata;
};

#endif

