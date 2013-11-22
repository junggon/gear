//====================================================================================================
//
//      Title :   liegroup.h
//
//      Note  :   This is a modified version of the original work by Jinwook Kim (see below).
//                Modified and maintained by Junggon Kim (junggon@gmail.com)
//
//====================================================================================================

//////////////////////////////////////////////////////////////////////////////////
//
//		title		:	LieGroup.h 
//						
//		version		:	v0.97 
//		author		:	Jinwook Kim (zinook@plaza.snu.ac.kr)
//		last update	:	2001.9.6
//
//		Note		:	
//						v0.95 library title changed : robotics.* -> liegroup.*
//						v0.95 is independent of RMatrix
//						v0.95 supports SO3 class
//						v0.95 Inertia class uses smaller number of member variables
//						v0.95 supports friend functions InvAd, InvdAd
//						v0.95 supports /, % operators in SE3 class
//						v0.97 supports articulated inertia class
//						v0.97 supports dad(V, J) - fast calc. of dad(V, J * V)
//
//////////////////////////////////////////////////////////////////////////////////

#ifndef _LIEGROUP_
#define _LIEGROUP_

#include <iostream>
#include <math.h>
#include "greal.h"

class se3;
class dse3;
class SO3;
class SE3;
class Inertia;
class AInertia;

class Vec3
{
private:
	gReal _v[3];
public:
	// constructors
	Vec3() {}
	Vec3(gReal d) { _v[0] = _v[1] = _v[2] = d; }
	Vec3(const gReal v[]) { _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; }
	Vec3(gReal v0, gReal v1, gReal v2) { _v[0] = v0; _v[1] = v1; _v[2] = v2; }
	// operators
	const Vec3 &operator + (void) const { return *this; }						// unary plus 
	Vec3 operator - (void) const { return Vec3(-_v[0],-_v[1],-_v[2]); }		// unary minus 
	gReal &operator [] (int i) { return _v[i]; }
	gReal operator [] (int i) const { return _v[i]; }
	Vec3 &operator = (const Vec3 &v) { _v[0] = v._v[0]; _v[1] = v._v[1]; _v[2] = v._v[2]; return *this; }
	Vec3 &operator = (const gReal v[]) { _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; return *this; }
	Vec3 &operator += (const Vec3 &v) { _v[0] += v._v[0]; _v[1] += v._v[1]; _v[2] += v._v[2];	return *this; }
	Vec3 &operator -= (const Vec3 &v)	{ _v[0] -= v._v[0]; _v[1] -= v._v[1]; _v[2] -= v._v[2]; return *this; }
	Vec3 &operator *= (gReal d) { _v[0] *= d; _v[1] *= d; _v[2] *= d; return *this; }
	Vec3 operator * (gReal d) const { return Vec3(d * _v[0], d * _v[1], d * _v[2]); }
	Vec3 operator + (const Vec3 &v) const { return Vec3(_v[0] + v._v[0], _v[1] + v._v[1], _v[2] + v._v[2]); }
	Vec3 operator - (const Vec3 &v) const { return Vec3(_v[0] - v._v[0], _v[1] - v._v[1], _v[2] - v._v[2]); }
	// methods
	void SetZero(void) { _v[0] = _v[1] = _v[2] = 0.0; }
	gReal Normalize() { gReal mag = Norm(*this); _v[0] /= mag; _v[1] /= mag; _v[2] /= mag; return mag; }
	gReal *GetArray(void) { return _v; }
	const gReal *GetArray(void) const { return _v; }
	gReal &GetAt(int i) { return _v[i]; }
	
	// friend classes
	friend class se3;
	friend class dse3;
	friend class SO3;
	friend class SE3;

	// friend functions
	friend std::ostream &operator << (std::ostream &os, const Vec3 &v);	// std::ostream standard output
	friend Vec3 operator * (gReal d, const Vec3 &v) { return Vec3(d * v._v[0], d * v._v[1], d * v._v[2]); }
	friend gReal Norm(const Vec3 &v) { return sqrt(v._v[0] * v._v[0] + v._v[1] * v._v[1] + v._v[2] * v._v[2]); }
	friend gReal SquareSum(const Vec3 &v) { return v._v[0] * v._v[0] + v._v[1] * v._v[1] + v._v[2] * v._v[2]; }
	friend Vec3 Cross(const Vec3 &p, const Vec3 &q) { return Vec3(p._v[1] * q._v[2] - p._v[2] * q._v[1], p._v[2] * q._v[0] - p._v[0]*q._v[2], p._v[0] * q._v[1] - p._v[1]*q._v[0]); }
	friend gReal Inner(const Vec3 &p, const Vec3 &q) { return (p._v[0] * q._v[0] + p._v[1] * q._v[1] + p._v[2] * q._v[2]); }
	friend se3 Ad(const Vec3 &p, const se3 &s);	
	friend SO3 Exp(const Vec3 &w);
	friend SO3 EulerZYX(const Vec3 &x);
	friend SO3 EulerZYZ(const Vec3 &x);
	friend SO3 EulerXYZ(const Vec3 &x);
	friend SO3 EulerZXY(const Vec3 &x);
};

class se3
{
private:
	gReal _w[6];	// upper three : angular velocity,  lower three : linear velocity
public:
	// constructors
	se3() {}
	se3(gReal k) { _w[0] = _w[1] = _w[2] = _w[3] = _w[4] = _w[5] = k; }
	se3(gReal w0, gReal w1, gReal w2, gReal w3, gReal w4, gReal w5) { _w[0] = w0;	_w[1] = w1;	_w[2] = w2; _w[3] = w3;	_w[4] = w4;	_w[5] = w5; }
	se3(const se3 &s) { _w[0] = s._w[0]; _w[1] = s._w[1]; _w[2] = s._w[2]; _w[3] = s._w[3]; _w[4] = s._w[4]; _w[5] = s._w[5]; }
	se3(const gReal *s) { _w[0] = s[0]; _w[1] = s[1]; _w[2] = s[2]; _w[3] = s[3]; _w[4] = s[4]; _w[5] = s[5]; }
	se3(const Vec3 &w, const Vec3 &v) { _w[0] = w._v[0]; _w[1] = w._v[1]; _w[2] = w._v[2]; _w[3] = v._v[0]; _w[4] = v._v[1]; _w[5] = v._v[2]; }
	// operators
	const se3 &operator + (void) const { return *this; }
	se3 operator - (void) const { return se3(-_w[0], -_w[1], -_w[2], -_w[3], -_w[4], -_w[5]); }
	se3 &operator = (const se3 &s) { _w[0] = s._w[0]; _w[1] = s._w[1]; _w[2] = s._w[2]; _w[3] = s._w[3]; _w[4] = s._w[4]; _w[5] = s._w[5]; return *this; }
	se3 &operator = (gReal d) { _w[0] = d; _w[1] = d; _w[2] = d; _w[3] = d; _w[4] = d; _w[5] = d; return *this; }
	se3 &operator += (const se3 &s) { _w[0] += s._w[0]; _w[1] += s._w[1]; _w[2] += s._w[2]; _w[3] += s._w[3]; _w[4] += s._w[4]; _w[5] += s._w[5]; return *this; }
	se3 &operator += (const gReal *s) { _w[0] += s[0]; _w[1] += s[1]; _w[2] += s[2]; _w[3] += s[3]; _w[4] += s[4]; _w[5] += s[5]; return *this; }
	se3 &operator -= (const se3 &s) { _w[0] -= s._w[0]; _w[1] -= s._w[1]; _w[2] -= s._w[2]; _w[3] -= s._w[3]; _w[4] -= s._w[4]; _w[5] -= s._w[5];	return *this; }
	se3 &operator *= (gReal d) { _w[0] *= d; _w[1] *= d; _w[2] *= d; _w[3] *= d; _w[4] *= d; _w[5] *= d; return *this; }
	se3 &operator /= (gReal d) { d = (gReal)1.0 / d; _w[0] *= d; _w[1] *= d; _w[2] *= d; _w[3] *= d; _w[4] *= d; _w[5] *= d; return *this; }
	se3 operator + (const se3 &s) const { return se3(_w[0] + s._w[0], _w[1] + s._w[1], _w[2] + s._w[2], _w[3] + s._w[3], _w[4] + s._w[4], _w[5] + s._w[5]); }
	se3 operator - (const se3 &s) const { return se3(_w[0] - s._w[0], _w[1] - s._w[1], _w[2] - s._w[2], _w[3] - s._w[3], _w[4] - s._w[4], _w[5] - s._w[5]); }
	se3 operator * (gReal d) const { return se3(d * _w[0], d * _w[1], d * _w[2], d * _w[3], d * _w[4], d * _w[5]); }
	se3 operator / (gReal d) const { d = (gReal)1.0 / d; return se3(d * _w[0], d * _w[1], d * _w[2], d * _w[3], d * _w[4], d * _w[5]); }
	gReal &operator [] (int i) { return _w[i]; }
	gReal operator [] (int i) const { return _w[i]; }
	// methods
	gReal *GetArray(void) { return _w; }
	const gReal *GetArray(void) const { return _w; }
	void SetZero(void) { _w[0] = _w[1] = _w[2] = _w[3] = _w[4] = _w[5] = (gReal)0.0; }
	gReal InnerProductWith(const gReal *s) { return (_w[0] * s[0] + _w[1] * s[1] + _w[2] * s[2] + _w[3] * s[3] + _w[4] * s[4] + _w[5] * s[5]); }
	Vec3 GetW() { return Vec3(_w[0],_w[1],_w[2]); }
	Vec3 GetV() { return Vec3(_w[3],_w[4],_w[5]); }
	Vec3 GetW() const { return Vec3(_w[0],_w[1],_w[2]); }
	Vec3 GetV() const { return Vec3(_w[3],_w[4],_w[5]); }
	void set_Ad(const SE3 &T, const se3 &s);		// *this = Ad(T,s)
	void set_Ad(const SE3 &T, const gReal *s);
	void set_Ad(const gReal *T, const gReal *s);
	void add_Ad(const SE3 &T, const se3 &s);		// *this = Ad(T,s)
	void add_Ad(const SE3 &T, const gReal *s);
	void add_Ad(const gReal *T, const gReal *s);
	void set_ad(const se3 &s1, const se3 &s2);		// *this = ad(s1,s2)
	void set_ad(const se3 &s1, const gReal *s2);
	void set_ad(const gReal *s1, const gReal *s2);
	void add_ad(const se3 &s1, const se3 &s2);		// *this = ad(s1,s2)
	void add_ad(const se3 &s1, const gReal *s2);
	void add_ad(const gReal *s1, const gReal *s2);

	// friend classes
	friend class SE3;
	friend class Inertia;
	friend class AInertia;

	// friend functions
	friend std::ostream &operator << (std::ostream &os, const se3 &s);	// std::ostream standard output
	friend se3 operator * (gReal d, const se3 &s) { return se3(d * s._w[0], d * s._w[1], d * s._w[2], d * s._w[3], d * s._w[4], d * s._w[5]); }
	friend gReal operator * (const dse3 &t, const se3 &s);
	friend gReal operator * (const se3 &s, const dse3 &t);
	friend SE3 Exp(const se3 &s);
	friend SE3 Exp(const se3 &s, gReal theta);
	friend se3 Log(const SE3 &T);
	friend se3 Ad(const SE3 &T, const se3 &s);
//	friend se3 InvAd(const SE3 &T, const se3 &s);
	friend se3 Ad(const Vec3 &p, const se3 &s);
	friend se3 ad(const se3 &s1, const se3 &s2);
	friend dse3 dad(const se3 &s, const dse3 &t);
	friend dse3 dad(const se3 &V, const Inertia &J);
	friend Vec3 GetW(const se3 &s) { return Vec3(s._w[0], s._w[1], s._w[2]); }
	friend Vec3 GetV(const se3 &s) { return Vec3(s._w[3], s._w[4], s._w[5]); }
	friend se3 LieBracket(const se3 &s1, const se3 &s2) { return ad(s1, s2); }
	friend se3 InvSkew(const SE3 &T); // invskew(T - I)
	friend gReal SquareSum(const se3 &s) { return s._w[0] * s._w[0] + s._w[1] * s._w[1] + s._w[2] * s._w[2] + s._w[3] * s._w[3] + s._w[4] * s._w[4] + s._w[5] * s._w[5]; }
	friend void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re = J*s
	friend void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s);				// re = J*s
	friend void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re += J*s
	friend void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s);				// re += J*s
};

class dse3
{
private:
	gReal _m[6];
public:
	// constructors
	dse3() {}
	dse3(gReal k) { _m[0] = _m[1] = _m[2] = _m[3] = _m[4] = _m[5] = k; }
	dse3(gReal m0, gReal m1, gReal m2, gReal m3, gReal m4, gReal m5) { _m[0] = m0; _m[1] = m1; _m[2] = m2; _m[3] = m3; _m[4] = m4; _m[5] = m5; }
	dse3(const dse3 &t) { _m[0] = t._m[0]; _m[1] = t._m[1]; _m[2] = t._m[2];	_m[3] = t._m[3]; _m[4] = t._m[4]; _m[5] = t._m[5]; }
	dse3(const gReal *t) { _m[0] = t[0]; _m[1] = t[1]; _m[2] = t[2]; _m[3] = t[3]; _m[4] = t[4]; _m[5] = t[5]; }
	dse3(const Vec3 &m, const Vec3 &f) { _m[0] = m._v[0]; _m[1] = m._v[1]; _m[2] = m._v[2];	_m[3] = f._v[0]; _m[4] = f._v[1]; _m[5] = f._v[2]; }
	// operators
	const dse3 &operator + (void) const { return *this; }
	dse3 operator - (void) const { return dse3(-_m[0], -_m[1], -_m[2], -_m[3], -_m[4], -_m[5]); }
	dse3 &operator = (const dse3 &t) { _m[0] = t._m[0]; _m[1] = t._m[1]; _m[2] = t._m[2]; _m[3] = t._m[3]; _m[4] = t._m[4]; _m[5] = t._m[5]; return *this; }
	dse3 &operator = (gReal d) { _m[0] = d; _m[1] = d; _m[2] = d; _m[3] = d; _m[4] = d; _m[5] = d; return *this; }
	dse3 &operator += (const dse3 &t) { _m[0] += t._m[0]; _m[1] += t._m[1]; _m[2] += t._m[2]; _m[3] += t._m[3]; _m[4] += t._m[4]; _m[5] += t._m[5]; return *this; }
	dse3 &operator += (const gReal *t) { _m[0] += t[0]; _m[1] += t[1]; _m[2] += t[2]; _m[3] += t[3]; _m[4] += t[4]; _m[5] += t[5]; return *this; }
	dse3 &operator -= (const dse3 &t) { _m[0] -= t._m[0]; _m[1] -= t._m[1]; _m[2] -= t._m[2]; _m[3] -= t._m[3]; _m[4] -= t._m[4]; _m[5] -= t._m[5]; return *this; }
	dse3 &operator *= (gReal d) { _m[0] *= d; _m[1] *= d; _m[2] *= d; _m[3] *= d; _m[4] *= d; _m[5] *= d; return *this; }
	dse3 &operator /= (gReal d) { d = (gReal)1.0 / d; _m[0] *= d; _m[1] *= d; _m[2] *= d; _m[3] *= d; _m[4] *= d; _m[5] *= d; return *this; }
	dse3 operator + (const dse3 &t) const { return dse3(_m[0] + t._m[0], _m[1] + t._m[1], _m[2] + t._m[2], _m[3] + t._m[3], _m[4] + t._m[4], _m[5] + t._m[5]); }	
	dse3 operator - (const dse3 &t) const { return dse3(_m[0] - t._m[0], _m[1] - t._m[1], _m[2] - t._m[2], _m[3] - t._m[3], _m[4] - t._m[4], _m[5] - t._m[5]); }	
	dse3 operator * (gReal d) const { return dse3(d * _m[0], d * _m[1], d * _m[2], d * _m[3], d * _m[4], d * _m[5]); }
	gReal &operator [] (int i) { return _m[i]; }
	gReal operator [] (int i) const { return _m[i]; }
	// methods
	void SetZero(void) { _m[0] = _m[1] = _m[2] = _m[3] = _m[4] = _m[5] = (gReal)0.0; }
	gReal InnerProductWith(const gReal *s) { return (_m[0] * s[0] + _m[1] * s[1] + _m[2] * s[2] + _m[3] * s[3] + _m[4] * s[4] + _m[5] * s[5]); }
	Vec3 GetM() { return Vec3(_m[0],_m[1],_m[2]); }
	Vec3 GetF() { return Vec3(_m[3],_m[4],_m[5]); }
	gReal *GetArray(void) { return _m; }
	const gReal *GetArray(void) const { return _m; }

	// friend class
	friend class Inertia;
	friend class AInertia;

	// friend functions
	friend std::ostream &operator << (std::ostream &os, const dse3 &t);	// std::ostream standard output
	friend dse3 operator * (gReal d, const dse3 &t) { return dse3(d * t._m[0], d * t._m[1], d * t._m[2], d * t._m[3], d * t._m[4], d * t._m[5]); }
	friend gReal operator * (const dse3 &t, const se3 &s) { return (t._m[0] * s._w[0] + t._m[1] * s._w[1] + t._m[2] * s._w[2] + t._m[3] * s._w[3] + t._m[4] * s._w[4] + t._m[5] * s._w[5]); }
	friend gReal operator * (const se3 &s, const dse3 &t) { return (t * s); }
	friend dse3 dAd(const SE3 &T, const dse3 &t);
//	friend dse3 InvdAd(const SE3 &T, const dse3 &t);
	friend dse3 dad(const se3 &s, const dse3 &t);	
	friend dse3 dad(const se3 &V, const Inertia &J);
	friend AInertia KroneckerProduct(const dse3 &x, const dse3 &y);
	friend gReal SquareSum(const dse3 &t) { return t._m[0] * t._m[0] + t._m[1] * t._m[1] + t._m[2] * t._m[2] + t._m[3] * t._m[3] + t._m[4] * t._m[4] + t._m[5] * t._m[5]; }
	friend Vec3 GetM(const dse3 &t) { return Vec3(t._m[0], t._m[1], t._m[2]); }
	friend Vec3 GetF(const dse3 &t) { return Vec3(t._m[3], t._m[4], t._m[5]); }
	friend void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re = J*s
	friend void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s);			// re = J*s
	friend void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re += J*s
	friend void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s);			// re += J*s
};

class SO3
{
private:
	gReal _R[9];
public:
	// constructors
	SO3() { _R[0] = _R[4] = _R[8] = 1.0; _R[1] = _R[2] = _R[3] = _R[5] = _R[6] = _R[7] = 0.0; }
	SO3(const SO3 &R) { _R[0] = R._R[0]; _R[3] = R._R[3]; _R[6] = R._R[6]; _R[1] = R._R[1]; _R[4] = R._R[4]; _R[7] = R._R[7]; _R[2] = R._R[2]; _R[5] = R._R[5]; _R[8] = R._R[8]; }
	SO3(const gReal R[]) { _R[0] = R[0]; _R[3] = R[3]; _R[6] = R[6]; _R[1] = R[1]; _R[4] = R[4]; _R[7] = R[7]; _R[2] = R[2]; _R[5] = R[5]; _R[8] = R[8]; }
	SO3(gReal R0, gReal R1, gReal R2, gReal R3, gReal R4, gReal R5, gReal R6, gReal R7, gReal R8) { _R[0] = R0; _R[1] = R1; _R[2] = R2; _R[3] = R3; _R[4] = R4; _R[5] = R5; _R[6] = R6; _R[7] = R7; _R[8] = R8; }
	SO3(const Vec3 &ex, const Vec3 &ey, const Vec3 &ez) { _R[0] = ex[0]; _R[1] = ex[1]; _R[2] = ex[2]; _R[3] = ey[0]; _R[4] = ey[1]; _R[5] = ey[2]; _R[6] = ez[0]; _R[7] = ez[1]; _R[8] = ez[2]; }
	// operators
	gReal operator () (int i, int j) const { return _R[i+3*j]; }
	gReal &operator () (int i, int j) { return _R[i+3*j]; }
	gReal operator [] (int i) const { return _R[i]; }
	gReal &operator [] (int i) { return _R[i]; }
	SO3 &operator = (const SO3 &R)  { _R[0] = R._R[0]; _R[3] = R._R[3]; _R[6] = R._R[6]; _R[1] = R._R[1]; _R[4] = R._R[4]; _R[7] = R._R[7]; _R[2] = R._R[2]; _R[5] = R._R[5]; _R[8] = R._R[8]; return *this; }
	SO3 &operator *= (const SO3 &R);
	SO3 operator * (const SO3 &R) const;
	Vec3 operator * (const Vec3 &p) const { return Vec3(_R[0] * p._v[0] + _R[3] * p._v[1] + _R[6] * p._v[2], _R[1] * p._v[0] + _R[4] * p._v[1] + _R[7] * p._v[2], _R[2] * p._v[0] + _R[5] * p._v[1] + _R[8] * p._v[2]); }
	SO3 operator ~ (void) const { return SO3(_R[0], _R[3], _R[6], _R[1], _R[4], _R[7], _R[2], _R[5], _R[8]); }
	// methods
	void SetIdentity(void) { _R[0] = _R[4] = _R[8] = 1.0; _R[1] = _R[2] = _R[3] = _R[5] = _R[6] = _R[7] = 0.0; }
	gReal *GetArray(void) { return _R; }
	const gReal *GetArray(void) const { return _R; }
	
	// friend class
	friend class SE3;

	// friend functions
	friend std::ostream &operator << (std::ostream &os, const SO3 &R);	// std::ostream standard output
	friend SO3 Inv(const SO3 &R) { return SO3(R._R[0], R._R[3], R._R[6], R._R[1], R._R[4], R._R[7], R._R[2], R._R[5], R._R[8]); }
	friend SO3 Exp(gReal w0, gReal w1, gReal w2);
	friend SO3 Exp(const Vec3 &w);
	friend Vec3 Log(const SO3 &R);
	friend SO3 RotX(gReal theta);
	friend SO3 RotY(gReal theta);
	friend SO3 RotZ(gReal theta);
	friend SO3 EulerZYX(const Vec3 &x);		// singularity : x[1] = -+ 0.5*PI
	friend Vec3 iEulerZYX(const SO3 &R);
	friend SO3 EulerZYZ(const Vec3 &x);		// singularity : x[1] = 0, PI
	friend Vec3 iEulerZYZ(const SO3 &R);
	friend SO3 EulerZXY(const Vec3 &x);
	friend Vec3 iEulerZXY(const SO3 &R);
	friend SO3 EulerXYZ(const Vec3 &x);
	friend Vec3 iEulerXYZ(const SO3 &R);
	friend SO3 Quat(gReal *quat);						// quaternion(quat[4]) --> SO3
	friend void iQuat(const SO3 &R, gReal *quat);		// SO3 --> quaternion(quat[4])
	friend bool isSO3(gReal *R, gReal eps = 1E-6);	// is R[9] is a rotation matrix?
	friend SO3 GetRotationWithZAxis(Vec3 axis);			// return a rotation matrix whose z-axis is aligned to axis
};

class SE3
{
private:
	gReal _T[16];	// column order 4 X 4 homogeneous transformation matrix
public:
	// constructors
	SE3() { _T[0] = _T[5] = _T[10] = _T[15] = 1.0; _T[1] = _T[2] = _T[3] = _T[4] = _T[6] = _T[7] = _T[8] = _T[9] = _T[11] = _T[12] = _T[13] = _T[14] = 0.0; }
	SE3(const SE3 &T) { _T[0] = T._T[0]; _T[1] = T._T[1]; _T[2] = T._T[2]; _T[4] = T._T[4]; _T[5] = T._T[5]; _T[6] = T._T[6]; _T[8] = T._T[8]; _T[9] = T._T[9]; _T[10] = T._T[10]; _T[12] = T._T[12]; _T[13] = T._T[13]; _T[14] = T._T[14]; _T[3] = _T[7] = _T[11] = 0.0; _T[15] = 1.0; }
	SE3(gReal T0, gReal T1, gReal T2, gReal T4, gReal T5, gReal T6, gReal T8, gReal T9, gReal T10, gReal T12, gReal T13, gReal T14) { _T[0] = T0; _T[1] = T1; _T[2] = T2; _T[4] = T4; _T[5] = T5; _T[6] = T6; _T[8] = T8; _T[9] = T9; _T[10] = T10; _T[12] = T12; _T[13] = T13; _T[14] = T14; _T[3] = _T[7] = _T[11] = 0.0; _T[15] = 1.0; }
	SE3(const SO3 &R, const Vec3 &p) { _T[0] = R._R[0]; _T[4] = R._R[3]; _T[8] = R._R[6]; _T[1] = R._R[1]; _T[5] = R._R[4]; _T[9] = R._R[7]; _T[2] = R._R[2]; _T[6] = R._R[5]; _T[10] = R._R[8]; _T[12] = p._v[0]; _T[13] = p._v[1]; _T[14] = p._v[2]; _T[3] = _T[7] = _T[11] = 0.0; _T[15] = 1.0; }
	SE3(const SO3 &R) { _T[0] = R._R[0]; _T[4] = R._R[3]; _T[8] = R._R[6]; _T[1] = R._R[1]; _T[5] = R._R[4]; _T[9] = R._R[7]; _T[2] = R._R[2]; _T[6] = R._R[5]; _T[10] = R._R[8]; _T[3] = _T[7] = _T[11] = _T[12] = _T[13] = _T[14] = 0.0; _T[15] = 1.0; }
	SE3(const Vec3 &p) { _T[0] = _T[5] = _T[10] = _T[15] = 1.0; _T[1] = _T[2] = _T[3] = _T[4] = _T[6] = _T[7] = _T[8] = _T[9] = _T[11] = 0.0; _T[12] = p._v[0]; _T[13] = p._v[1]; _T[14] = p._v[2]; }
	SE3(const gReal T[]) { _T[0] = T[0]; _T[1] = T[1]; _T[2] = T[2]; _T[4] = T[4]; _T[5] = T[5]; _T[6] = T[6]; _T[8] = T[8]; _T[9] = T[9]; _T[10] = T[10]; _T[12] = T[12]; _T[13] = T[13]; _T[14] = T[14]; _T[3] = _T[7] = _T[11] = 0.0; _T[15] = 1.0; }
	// operators
	gReal operator () (int i, int j) const { return _T[i+4*j]; }
	gReal &operator () (int i, int j) { return _T[i+4*j]; }
	gReal operator [] (int i) const { return _T[i]; }
	gReal &operator [] (int i) { return _T[i]; }
	SE3 &operator = (const SE3 &T) { _T[0] = T._T[0]; _T[1] = T._T[1]; _T[2] = T._T[2]; _T[4] = T._T[4]; _T[5] = T._T[5]; _T[6] = T._T[6]; _T[8] = T._T[8]; _T[9] = T._T[9]; _T[10] = T._T[10]; _T[12] = T._T[12]; _T[13] = T._T[13]; _T[14] = T._T[14]; return *this; }
	SE3 operator * (const SE3 &T) const;
	SE3 operator / (const SE3 &T) const;	// *this * Inv(T), note that A / B * C != A / ( B * C )
	SE3 operator % (const SE3 &T) const;	// Inv(*this) * T, note that A * B % C != A * ( B % C )
	Vec3 operator * (const Vec3 &p) const;
	SE3 &operator *= (const SE3 &T);
	SE3 &operator /= (const SE3 &T);
	SE3 &operator %= (const SE3 &T);
	// methods
	void SetIdentity(void) { _T[0] = _T[5] = _T[10] = _T[15] = 1.0; _T[1] = _T[2] = _T[3] = _T[4] = _T[6] = _T[7] = _T[8] = _T[9] = _T[11] = _T[12] = _T[13] = _T[14] = 0.0; }
	void SetInvOf(const SE3 &T); // *this = Inv(T)
	SE3 &SetRotation(const SO3 &R);
	SE3 &SetPosition(const Vec3 &Pos);
	SE3 &Translate(const Vec3 &Pos);
	SE3 &Rotate(const SO3 &R);
	Vec3 GetPosition(void) const { return Vec3(_T[12], _T[13], _T[14]); }
	SO3 GetRotation(void) const { return SO3(_T[0], _T[1], _T[2], _T[4], _T[5], _T[6], _T[8], _T[9], _T[10]); }
	gReal *GetArray(void) { return _T; }
	const gReal *GetArray(void) const { return _T; }
	
	// friend functions
	friend std::ostream &operator << (std::ostream &os, const SE3 &T);	// std::ostream standard output
	friend SE3 Inv(const SE3 &T);
	friend SE3 Exp(const se3 &s);
	friend SE3 Exp(const se3 &s, gReal theta);
	friend se3 Log(const SE3 &T);
	friend se3 Ad(const SE3 &T, const se3& s);											// return Ad(T,s)
	friend dse3 dAd(const SE3 &T, const dse3 &t);										// return dAd(T,t)
	friend void matSet_Ad(gReal *re, const SE3 &T, const gReal *s);					// re = Ad(T,s)
	friend void matSet_Ad(gReal *re, const SE3 &T, const gReal *s, int num);			// re = [Ad(T,s_0), ..., Ad(T,s_(num-1))] where s_i = (s[6*i+0],...,s[6*i+5])
	friend void matSet_Ad_minus(gReal *re, const SE3 &T, const gReal *s);				// re = -Ad(T,s)
	friend void matSet_Ad_minus(gReal *re, const SE3 &T, const gReal *s, int num);	// re = -[dAd(T,s_0), ..., dAd(T,s_(num-1))] where s_i = (s[6*i+0],...,s[6*i+5])
	friend void matSet_dAd(gReal *re, const SE3 &T, const gReal *t);					// re = dAd(T,t)
	friend void matSet_dAd(gReal *re, const SE3 &T, const gReal *t, int num);			// re = [dAd(T,t_0), ..., dAd(T,t_(num-1))] where t_i = (t[6*i+0],...,t[6*i+5])
	friend void matSet_dAd_minus(gReal *re, const SE3 &T, const gReal *t);			// re = -dAd(T,t)
	friend void matSet_dAd_minus(gReal *re, const SE3 &T, const gReal *t, int num);	// re = -[Ad(T,t_0), ..., Ad(T,t_(num-1))] where t_i = (t[6*i+0],...,t[6*i+5])
//	friend se3 InvAd(const SE3 &T, const se3& s);
//	friend dse3 InvdAd(const SE3 &T, const dse3 &t);
	friend se3 InvSkew(const SE3 &T)  { return se3((gReal)0.5 * (T._T[6] - T._T[9]), (gReal)0.5 * (T._T[8] - T._T[2]), (gReal)0.5 * (T._T[1] - T._T[4]), T._T[12], T._T[13], T._T[14]); }	// invskew(T - I)
	
	// friend class
	friend class se3;
	friend class Inertia;
	friend class AInertia;
};

class Inertia
{
public:							// Inertia = [I, [r]; -[r], m*eye(3)]
	gReal _I[6], _r[3], _m;	// _I[0] = Ixx, _I[1] = Iyy, _I[2] = Izz, _I[3] = Ixy, _I[4] = Ixz, _I[5] = Iyz
public:
	// constructors
	Inertia() { _I[0] = _I[1] = _I[2] = _I[3] = _I[4] = _I[5] = _r[0] = _r[1] = _r[2] = _m = (gReal)0.0; }
	Inertia(gReal m) { _m = _I[0] = _I[1] = _I[2] = m; _I[3] = _I[4] = _I[5] = _r[0] = _r[1] = _r[2] = (gReal)0.0; }
	Inertia(gReal mass, gReal Ixx, gReal Iyy, gReal Izz);
	Inertia(const Inertia &J) { _I[0] = J._I[0]; _I[1] = J._I[1]; _I[2] = J._I[2]; _I[3] = J._I[3]; _I[4] = J._I[4]; _I[5] = J._I[5]; _r[0] = J._r[0]; _r[1] = J._r[1]; _r[2] = J._r[2]; _m = J._m; }
	// operator
	Inertia &operator = (const Inertia &J) { _I[0] = J._I[0]; _I[1] = J._I[1]; _I[2] = J._I[2]; _I[3] = J._I[3]; _I[4] = J._I[4]; _I[5] = J._I[5]; _r[0] = J._r[0]; _r[1] = J._r[1]; _r[2] = J._r[2]; _m = J._m; return *this; }
	Inertia &operator *= (gReal d) { _I[0] *= d; _I[1] *= d; _I[2] *= d; _I[3] *= d; _I[4] *= d; _I[5] *= d; _r[0] *= d; _r[1] *= d; _r[2] *= d; _m *= d; return *this; }
	Inertia operator * (gReal d) const { Inertia J(*this); J *= d; return J; }
	dse3 operator * (const se3 &acceleration) const;
	dse3 operator * (const gReal *s) const;
	// methods
	void SetZero(void) { _m = _I[0] = _I[1] = _I[2] = _I[3] = _I[4] = _I[5] = _r[0] = _r[1] = _r[2] = (gReal)0.0; }
	void SetMass(gReal mass) { _m = mass; }
	gReal GetMass(void) { return _m; }
	Vec3 GetOffDiag() { return Vec3(_r); }
	void SetInertia(gReal Ixx, gReal Iyy, gReal Izz, gReal Ixy, gReal Ixz, gReal Iyz) { _I[0] = Ixx; _I[1] = Iyy; _I[2] = Izz; _I[3] = Ixy; _I[4] = Ixz; _I[5] = Iyz; }
	void SetOffDiag(const gReal r[]) { _r[0] = r[0]; _r[1] = r[1]; _r[2] = r[2]; }
	Inertia Transform(const SE3 &T) const;
	void ToArray(gReal J[]) const; // J = [I, [r]; -[r], m*eye(3)]
	void InvToArray(gReal invJ[]) const; // invJ = Inv([I, [r]; -[r], m*eye(3)])
	
	// freind class
	friend class AInertia;

	// friend function
	friend std::ostream &operator << (std::ostream &os, const Inertia &iner);	// std::ostream standard output
	friend AInertia operator + (const Inertia &A, const AInertia &B);
	friend AInertia operator - (const Inertia &A, const AInertia &B);
	friend dse3 dad(const se3 &V, const Inertia &J);
};

// acticulated inertia class
class AInertia
{
private:						
	gReal _J[21];
	// _J[] = symmetric matrix elements
	// [_J[ 0], _J[ 1], _J[ 2], _J[ 3], _J[ 4], _J[ 5] ]
	// [        _J[ 6], _J[ 7], _J[ 8], _J[ 9], _J[10] ]	
	// [                _J[11], _J[12], _J[13], _J[14] ]
	// [                        _J[15], _J[16], _J[17] ]
	// [                                _J[18], _J[19] ]
	// [                                        _J[20] ]

public:
	// constructors
	AInertia() { _J[0] = _J[1] = _J[2] = _J[3] = _J[4] = _J[5] = _J[6] = _J[7] = _J[8] = _J[9] = _J[10] = _J[11] = _J[12] = _J[13] = _J[14] = _J[15] = _J[16] = _J[17] = _J[18] = _J[19] = _J[20] = 0.0; }
	AInertia(gReal d) { _J[0] = _J[1] = _J[2] = _J[3] = _J[4] = _J[5] = _J[6] = _J[7] = _J[8] = _J[9] = _J[10] = _J[11] = _J[12] = _J[13] = _J[14] = _J[15] = _J[16] = _J[17] = _J[18] = _J[19] = _J[20] = d; }
	AInertia(gReal a0, gReal a1, gReal a2, gReal a3, gReal a4, gReal a5, gReal a6, gReal a7, gReal a8, gReal a9, gReal a10, gReal a11, gReal a12, gReal a13, gReal a14, gReal a15, gReal a16, gReal a17, gReal a18, gReal a19, gReal a20);
	AInertia(const AInertia &J); 
	AInertia(const Inertia &J);
	AInertia(const gReal *M);	// M[36] is a single array representing column-wise 6 x 6 symmetric matrix.
	// operator
	const AInertia &operator + (void) const { return *this; }
	AInertia operator - (void) const;
	dse3 operator * (const se3 &a) const;
	AInertia operator + (const AInertia &J) const;
	AInertia operator + (const Inertia &J) const;
	AInertia operator - (const AInertia &J) const;
	AInertia operator - (const Inertia &J) const;
	AInertia &operator = (const AInertia &J);
	AInertia &operator = (const Inertia &J);
	AInertia &operator += (const AInertia &J);
	AInertia &operator += (const Inertia &J);
	AInertia &operator -= (const AInertia &J);
	AInertia &operator -= (const Inertia &J);
	AInertia &operator += (const gReal *M); // M[36] is a single array representing column-wise 6 x 6 symmetric matrix.
	AInertia &operator -= (const gReal *M); 
	// methods
	void SetZero(void) { _J[0] = _J[1] = _J[2] = _J[3] = _J[4] = _J[5] = _J[6] = _J[7] = _J[8] = _J[9] = _J[10] = _J[11] = _J[12] = _J[13] = _J[14] = _J[15] = _J[16] = _J[17] = _J[18] = _J[19] = _J[20] = 0.0; }
	AInertia Transform_old(const SE3 &T) const;
	AInertia Transform(const SE3 &T) const;
	void AddTransform(const AInertia &J, const SE3 &T);
	void SubtractKroneckerProduct(const dse3 &x, const dse3 &y);		// (*this) -= KroneckerProduct(x,y)
	void SubtractKroneckerProduct(const gReal *x, const gReal *y);
	void SubstractAlphaXYt(gReal alpha, const gReal *x, const gReal *y); // (*this) -= alpha*x*~y where x and y are 6-dimensional vectors and alpha is a scalar
	void SubstractAlphaSSt(gReal alpha, const gReal *s);		// (*this) -= alpha*s*~s where s = 6-dimensional vector and alpha is a scalar 
	AInertia Transform_ad(const se3 &s) const;					// return (*this) * ad(s) + ~((*this) * ad(s))
	void AddTransform_ad(const AInertia &J, const se3 &s);		// (*this) += J * ad(s) + ~(J * ad(s))
	void SubstractTransform_ad(const AInertia &J, const se3 &s); // (*this) -= J * ad(s) + ~(J * ad(s))
	void ToArray(gReal I[]) const;

	// friend class
	friend class Inertia;

	// friend function	
	friend std::ostream &operator << (std::ostream &os, const AInertia &J);
	friend AInertia KroneckerProduct(const dse3 &x, const dse3 &y);								// return 0.5*(x*~y + y*~x)
	friend void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re = J*s
	friend void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s);			// re = J*s
	friend void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s);				// re = J*s
	friend void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s);			// re = J*s
	friend void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s, int num);	// re = [J*s_0, ..., J*s_(num-1)] where s_i = (s[6*i+0],...,s[6*i+5])
	friend void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re += J*s
	friend void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s);			// re += J*s
	friend void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s);				// re += J*s
	friend void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s);			// re += J*s
	friend void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s, int num);	// re += [J*s_0, ..., J*s_(num-1)] where s_i = (s[6*i+0],...,s[6*i+5])
};

se3 ad(const gReal *s1, const gReal *s2);
dse3 dAd(const SE3 &T, const dse3 &t);										// return dAd(T,t)

void set_Mult_Inertia_se3(gReal *re, const Inertia &I, const gReal *s);			// re = I*s
void set_Mult_Inertia_se3(gReal *re, const Inertia &I, const gReal *s, int num);	// re = [I*s_0, ..., I*s_(num-1)] where s_i = (s[6*i+0],...,s[6*i+5])

void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re = J*s
void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s);			// re = J*s
void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s);				// re = J*s
void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s);			// re = J*s
void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s, int num);	// re = [J*s_0, ..., J*s_(num-1)] where s_i = (s[6*i+0],...,s[6*i+5])
void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s);				// re += J*s
void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s);			// re += J*s
void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s);				// re += J*s
void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s);			// re += J*s
void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s, int num);	// re += [J*s_0, ..., J*s_(num-1)] where s_i = (s[6*i+0],...,s[6*i+5])

void matSet(gReal *re, const gReal *a, int n);			// re = a where size(a) = n
void matAdd(gReal *re, const gReal *a, int n);			// re += a where size(a) = n
void matSubtract(gReal *re, const gReal *a, int n);		// re -= a where size(a) = n
void matMult(gReal *re, gReal alpha, int n);				// re *= alpha where alpha is a scalar

void matSet(gReal *re, const gReal *a, int a_row, int a_col);			// re = A where A = [a] = (a_row x a_col)
void matAdd(gReal *re, const gReal *a, int a_row, int a_col);			// re += A where A = [a] = (a_row x a_col)
void matSubtract(gReal *re, const gReal *a, int a_row, int a_col);	// re -= A where A = [a] = (a_row x a_col)
void matMult(gReal *re, gReal alpha, int a_row, int a_col);			// re *= alpha where alpha is a scalar

void matSet_zero(gReal *re, int n);						// re = zero (re[0:n-1] = 0.0)

void matSet_eye(gReal *re, int r);						// re = eye(r,r) (r by r identity matrix)

void matSet_transpose(gReal *re, const gReal *a, int a_row, int a_col);	
															// re = ~A (transpose of A) where A = [a] = (a_row x a_col)

void matSet_multAB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re = A*B where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_row * b_col and a_col = b_row

void matSet_multAtB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re = A'*B where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_col * b_col and a_row = b_row

void matSet_multABt(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re = A*B' where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_row * b_row and a_col = b_col

void matAdd_multAB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re += A*B where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_row * b_col and a_col = b_row

void matAdd_multAtB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re += A'*B where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_col * b_col and a_row = b_row

void matAdd_multABt(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re += A*B' where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_row * b_row and a_col = b_col

void matSubtract_multAB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re -= A*B where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_row * b_col and a_col = b_row

void matSubtract_multAtB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re -= A'*B where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_col * b_col and a_row = b_row

void matSubtract_multABt(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col); 
															// re -= A*B' where A = [a] = (a_row x a_col) and B = [b] = (b_row x b_col).
															// requirements: size(re) = a_row * b_row and a_col = b_col

bool matSet_inv22(gReal *re, const gReal *a, gReal eps=0);		// re = inverse of A = [a] = (2 x 2)
bool matSet_inv22sym(gReal *re, const gReal *a, gReal eps=0);	// re = inverse of symmetric A = [a] = (2 x 2)
bool matSet_inv33(gReal *re, const gReal *a, gReal eps=0);		// re = inverse of A = [a] = (3 x 3)
bool matSet_inv33sym(gReal *re, const gReal *a, gReal eps=0);	// re = inverse of symmetric A = [a] = (3 x 3)
bool matSet_inv44(gReal *re, const gReal *a, gReal eps=0);		// re = inverse of A = [a] = (4 x 4)
bool matSet_inv44sym(gReal *re, const gReal *a, gReal eps=0);	// re = inverse of symmetric A = [a] = (4 x 4)
bool matSet_invNN(gReal *re, const gReal *a, int n);			// re = inverse of A = [a] = (n x n)
bool matSet_invNNfast(gReal *re, gReal *a, int n);				// re = inverse of A = [a] = (n x n) where n <= _MAX_SIZE_INVNNFAST
																//      Caution! A will be corrupted during the computation.

void matSet_ad(gReal *re, const gReal *s1, const gReal *s2);	// re = ad(s1,s2)
void matAdd_ad(gReal *re, const gReal *s1, const gReal *s2);	// re += ad(s1,s2)

void matSet_Ad(gReal *re, const SE3 &T, const gReal *s);					// re = Ad(T,s)
void matSet_Ad(gReal *re, const SE3 &T, const gReal *s, int num);			// re = [Ad(T,s_0), ..., Ad(T,s_(num-1))] where s_i = (s[6*i+0],...,s[6*i+5])
void matSet_Ad_minus(gReal *re, const SE3 &T, const gReal *s);			// re = -Ad(T,s)
void matSet_Ad_minus(gReal *re, const SE3 &T, const gReal *s, int num);	// re = -[dAd(T,s_0), ..., dAd(T,s_(num-1))] where s_i = (s[6*i+0],...,s[6*i+5])
void matSet_dAd(gReal *re, const SE3 &T, const gReal *t);					// re = dAd(T,t)
void matSet_dAd(gReal *re, const SE3 &T, const gReal *t, int num);		// re = [dAd(T,t_0), ..., dAd(T,t_(num-1))] where t_i = (t[6*i+0],...,t[6*i+5])
void matSet_dAd_minus(gReal *re, const SE3 &T, const gReal *t);			// re = -dAd(T,t)
void matSet_dAd_minus(gReal *re, const SE3 &T, const gReal *t, int num);	// re = -[Ad(T,t_0), ..., Ad(T,t_(num-1))] where t_i = (t[6*i+0],...,t[6*i+5])

// subfunctions
int __idamax(int n, gReal *dx);
void __dgefa(gReal *x, int lda, int n, int *jpvt, int &info);
void __dgesl(gReal *x, int lda, int n, int *jpvt, gReal *b, int job);

#include "liegroup.inl"

#endif