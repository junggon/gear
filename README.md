gear
====
Geometric Engine for Articulated Rigid-body simulation

GEAR is a C++ library for kinematics and dynamics of articulated rigid body systems. It uses efficient algorithms based on Lie group and screw theory to solve the following problems:

- Forward dynamics for calculating the resulting joint accelerations when the known joint torques or forces are applied to the mechanical system. 
- Inverse dynamics for calculating the joint torques that are necessary to make the system follow prescribed joint accelerations at the current state (position and velocity). 
- Hybrid dynamics for calculating the torques for prescribed joints, and the accelerations for unprescribed joints where 'prescribed' means joint acceleration is prescribed and 'unprescribed'  (or 'torque-specified') means joint torque is known or given. A joint can be either 'prescribed' or 'unprescribed'. Hybrid dynamics is a generalization of the forward and inverse dynamics, i.e., they can be regarded as the extreme cases of hybrid dynamics when all of the joints are 'unprescribed' and when all of the joints are 'prescribed' respectively. 
- Analytical derivatives of the forward/inverse/hybrid dynamics with respect to an arbitrary system parameter. 

GEAR uses a joint coordinate system to describe the degree of freedom of an articulated rigid body system. Thus, the joint constraints are inherently considered in the algorithms without having to consider them separately. Closed joint loops are automatically detected and handled as constraints. Additional constraints can also be added to the system.

How to compile
==============
CMake and a C++ compiler are necessary to build GEAR.

- On Linux, do make in the terminal to compile the library under gear/build_release.
- On Windows, execute runcmake_win.bat to create a Visual Studio solution file under gear/build. 
- On Mac, execute runcmake_xcode to create an XCode project file under gear/build_xcode. Alternatively, you can do make in the terminal to compile using gcc.

How to use
==========
See gear/examples for how to create your own simulation project. 

