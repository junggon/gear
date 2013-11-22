gear
====
Geometric Engine for Articulated Rigid-body simulation

GEAR is a C++ library for kinematics and dynamics of articulated rigid body systems. It uses effecient algorithms based on Lie group and screw theory to solve the following problems:

- Forward dynamics for calculating the resulting joint accelerations when the known joint torques or forces are applied to the mechanical system. 
- Inverse dynamics for calculating the joint torques that are necessary to make the system follow prescribed joint accelerations at the current state (position and velocity). 
- Hybrid dynamics for calculating the torques for prescribed joints, and the accelerations for unprescribed joints where 'prescribed' means joint acceleration is prescribed and 'unprescribed'  (or 'torque-specified') means joint torque is known or given. A joint can be either 'prescribed' or 'unprescribed'. Hybrid dynamics is a generalization of the forward and inverse dynamics, i.e., they can be regarded as the extreme cases of hybrid dynamics when all of the joints are 'unprescribed' and when all of the joints are 'prescribed' respectively. 
- Analytical derivatives of the forward/inverse/hybrid dynamics with respect to an arbitrary system parameter. 

