//====================================================================================================
//
//      Title :   liegroup.inl
//
//      Note  :   This is a modified version of the original work by Jinwook Kim (see liegroup.h).
//                Modified and maintained by Junggon Kim (junggon@gmail.com)
//
//====================================================================================================

#define _PI_LIE		3.14159265358979
#define _SMALL_LIE	1E-9
#define _TINY_LIE	1E-12
#define _MAX_SIZE_INVNNFAST 6 // for matSet_invNNfast()

inline std::ostream &operator << (std::ostream &os, const Vec3 &v)
{
	os << "[ ";
	for ( int i = 0; i < 3; i++ )
	{	
		if ( v._v[i] >= 0.0 ) os << " ";
		os << v._v[i] << " ";
	}
	os << "];" << std::endl;
    return os;
}

inline std::ostream &operator << (std::ostream &os, const se3 &s)
{
	os << "[ ";
	for ( int i = 0; i < 6; i++ )
	{	
		if ( s._w[i] >= 0.0 ) os << " ";
		os << s._w[i] << " ";
	}
	os << "];" << std::endl;
    return os;
}

inline std::ostream &operator << (std::ostream &os, const dse3 &t)
{
	os << "[ ";
	for ( int i = 0; i < 6; i++ )
	{	
		if ( t._m[i] >= 0.0 ) os << " ";
		os << t._m[i] << " ";
	}
	os << "];" << std::endl;
    return os;
}

inline std::ostream &operator << (std::ostream &os, const SO3 &R)
{
	os << "[" << std::endl;
	for ( int i = 0; i < 3; i++ ) 
	{
		os << "  ";
		for ( int j = 0; j < 3; j++ )
		{
			if ( R._R[i+j*3] >= 0.0 ) os << " ";
			os << R._R[i+j*3];
			if ( j == 2 ) os << " ;" << std::endl;
			else os << "  ";
		}
	}
	os << "];" << std::endl;
	return os;
}

inline std::ostream &operator << (std::ostream &os, const SE3 &T)
{
	os << "[" << std::endl;
	for ( int i = 0; i < 4; i++ )
	{
		os << "  ";
		for ( int j = 0; j < 4; j++ )
		{
			if ( T._T[i+j*4] >= 0.0 ) os << " ";
			os << T._T[i+j*4];
			if ( j == 3 ) os << " ;" << std::endl;
			else os << "  ";
		}
	}
	os << "];" << std::endl;
	return os;
}

inline std::ostream &operator << (std::ostream &os, const Inertia &iner)
{
	os << "I : " << iner._I[0] << " " << iner._I[3] << " " << iner._I[4] << std::endl;
	os << "    " << iner._I[3] << " " << iner._I[1] << " " << iner._I[5] << std::endl;
	os << "    " << iner._I[4] << " " << iner._I[5] << " " << iner._I[2] << std::endl;
	os << "r : " << iner._r[0] << " " << iner._r[1] << " " << iner._r[2] << std::endl;
	os << "m : " << iner._m << std::endl;
    return os;
}

inline std::ostream &operator << (std::ostream &os, const AInertia &J)
{
	os << J._J[0] << " " << J._J[1] << " " << J._J[2] << " " << J._J[3] << " " << J._J[4] << " " << J._J[5] << std::endl;
	os << J._J[1] << " " << J._J[6] << " " << J._J[7] << " " << J._J[8] << " " << J._J[9] << " " << J._J[10] << std::endl;
	os << J._J[2] << " " << J._J[7] << " " << J._J[11] << " " << J._J[12] << " " << J._J[13] << " " << J._J[14] << std::endl;
	os << J._J[3] << " " << J._J[8] << " " << J._J[12] << " " << J._J[15] << " " << J._J[16] << " " << J._J[17] << std::endl;
	os << J._J[4] << " " << J._J[9] << " " << J._J[13] << " " << J._J[16] << " " << J._J[18] << " " << J._J[19] << std::endl;
	os << J._J[5] << " " << J._J[10] << " " << J._J[14] << " " << J._J[17] << " " << J._J[19] << " " << J._J[20] << std::endl;
	return os;
}

inline void matSet_Ad(gReal *re, const SE3 &T, const gReal *s)
{
	if ( s == NULL ) { matSet_zero(re, 6); return; }

	gReal tmp[3] = { T._T[0] * s[0] + T._T[4] * s[1] + T._T[8] * s[2], 
					T._T[1] * s[0] + T._T[5] * s[1] + T._T[9] * s[2], 
					T._T[2] * s[0] + T._T[6] * s[1] + T._T[10] * s[2] };
	re[0] = tmp[0];
	re[1] = tmp[1];
	re[2] = tmp[2];
	re[3] = T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * s[3] + T._T[4] * s[4] + T._T[8] * s[5];
	re[4] = T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * s[3] + T._T[5] * s[4] + T._T[9] * s[5];
	re[5] = T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * s[3] + T._T[6] * s[4] + T._T[10] * s[5];
}

inline void matSet_Ad_minus(gReal *re, const SE3 &T, const gReal *s)
{
	if ( s == NULL ) { matSet_zero(re, 6); return; }

	gReal tmp[3] = { T._T[0] * s[0] + T._T[4] * s[1] + T._T[8] * s[2], 
					T._T[1] * s[0] + T._T[5] * s[1] + T._T[9] * s[2], 
					T._T[2] * s[0] + T._T[6] * s[1] + T._T[10] * s[2] };
	re[0] = tmp[0]; 
	re[1] = tmp[1]; 
	re[2] = tmp[2]; 
	re[3] = T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * s[3] + T._T[4] * s[4] + T._T[8] * s[5]; 
	re[4] = T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * s[3] + T._T[5] * s[4] + T._T[9] * s[5]; 
	re[5] = T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * s[3] + T._T[6] * s[4] + T._T[10] * s[5]; 

	re[0] *= -1.0; re[1] *= -1.0; re[2] *= -1.0; re[3] *= -1.0; re[4] *= -1.0; re[5] *= -1.0;
}

inline void matSet_Ad(gReal *re, const SE3 &T, const gReal *s, int num)
{
	if ( s == NULL ) { matSet_zero(re, 6*num); return; }

	for (int i=0; i<num; i++) {
		matSet_Ad(&re[6*i], T, &s[6*i]);
	}
}

inline void matSet_Ad_minus(gReal *re, const SE3 &T, const gReal *s, int num)
{
	if ( s == NULL ) { matSet_zero(re, 6*num); return; }

	for (int i=0; i<num; i++) {
		matSet_Ad_minus(&re[6*i], T, &s[6*i]);
	}
}

inline void matSet_dAd(gReal *re, const SE3 &T, const gReal *t)
{
	if ( t == NULL ) { matSet_zero(re, 6); return; }

	gReal tmp[3] = { t[0] - T._T[13] * t[5] + T._T[14] * t[4], 
					t[1] - T._T[14] * t[3] + T._T[12] * t[5], 
					t[2] - T._T[12] * t[4] + T._T[13] * t[3] };

	re[0] = T._T[0] * tmp[0] + T._T[1] * tmp[1] + T._T[2] * tmp[2];
	re[1] = T._T[4] * tmp[0] + T._T[5] * tmp[1] + T._T[6] * tmp[2];
	re[2] = T._T[8] * tmp[0] + T._T[9] * tmp[1] + T._T[10] * tmp[2];
	re[3] = T._T[0] * t[3] + T._T[1] * t[4] + T._T[2] * t[5];
	re[4] = T._T[4] * t[3] + T._T[5] * t[4] + T._T[6] * t[5];
	re[5] = T._T[8] * t[3] + T._T[9] * t[4] + T._T[10] * t[5];
}

inline void matSet_dAd_minus(gReal *re, const SE3 &T, const gReal *t)
{
	if ( t == NULL ) { matSet_zero(re, 6); return; }

	gReal tmp[3] = { t[0] - T._T[13] * t[5] + T._T[14] * t[4], 
					t[1] - T._T[14] * t[3] + T._T[12] * t[5], 
					t[2] - T._T[12] * t[4] + T._T[13] * t[3] };

	re[0] = T._T[0] * tmp[0] + T._T[1] * tmp[1] + T._T[2] * tmp[2];
	re[1] = T._T[4] * tmp[0] + T._T[5] * tmp[1] + T._T[6] * tmp[2];
	re[2] = T._T[8] * tmp[0] + T._T[9] * tmp[1] + T._T[10] * tmp[2];
	re[3] = T._T[0] * t[3] + T._T[1] * t[4] + T._T[2] * t[5];
	re[4] = T._T[4] * t[3] + T._T[5] * t[4] + T._T[6] * t[5];
	re[5] = T._T[8] * t[3] + T._T[9] * t[4] + T._T[10] * t[5];

	re[0] *= -1.0; re[1] *= -1.0; re[2] *= -1.0; re[3] *= -1.0; re[4] *= -1.0; re[5] *= -1.0;
}

inline void matSet_dAd(gReal *re, const SE3 &T, const gReal *t, int num)
{
	if ( t == NULL ) { matSet_zero(re, 6*num); return; }

	for (int i=0; i<num; i++) {
		matSet_dAd(&re[6*i], T, &t[6*i]);
	}
}

inline void matSet_dAd_minus(gReal *re, const SE3 &T, const gReal *t, int num)
{
	if ( t == NULL ) { matSet_zero(re, 6*num); return; }

	for (int i=0; i<num; i++) {
		matSet_dAd_minus(&re[6*i], T, &t[6*i]);
	}
}

// *this = T * s * Inv(T)
inline void se3::set_Ad(const SE3 &T, const se3 &s)
{
	gReal tmp[3] = { T._T[0] * s._w[0] + T._T[4] * s._w[1] + T._T[8] * s._w[2], 
					T._T[1] * s._w[0] + T._T[5] * s._w[1] + T._T[9] * s._w[2], 
					T._T[2] * s._w[0] + T._T[6] * s._w[1] + T._T[10] * s._w[2] };
	_w[0] = tmp[0];
	_w[1] = tmp[1];
	_w[2] = tmp[2];
	_w[3] = T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * s._w[3] + T._T[4] * s._w[4] + T._T[8] * s._w[5];
	_w[4] = T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * s._w[3] + T._T[5] * s._w[4] + T._T[9] * s._w[5];
	_w[5] = T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * s._w[3] + T._T[6] * s._w[4] + T._T[10] * s._w[5];
}

inline void se3::set_Ad(const SE3 &T, const gReal *s)
{
	if ( s == NULL ) { SetZero(); return; }

	gReal tmp[3] = { T._T[0] * s[0] + T._T[4] * s[1] + T._T[8] * s[2], 
					T._T[1] * s[0] + T._T[5] * s[1] + T._T[9] * s[2], 
					T._T[2] * s[0] + T._T[6] * s[1] + T._T[10] * s[2] };
	_w[0] = tmp[0];
	_w[1] = tmp[1];
	_w[2] = tmp[2];
	_w[3] = T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * s[3] + T._T[4] * s[4] + T._T[8] * s[5];
	_w[4] = T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * s[3] + T._T[5] * s[4] + T._T[9] * s[5];
	_w[5] = T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * s[3] + T._T[6] * s[4] + T._T[10] * s[5];
}

inline void se3::set_Ad(const gReal *T, const gReal *s)
{
	if ( T == NULL || s == NULL ) { SetZero(); return; }

	gReal tmp[3] = { T[0] * s[0] + T[4] * s[1] + T[8] * s[2], 
					T[1] * s[0] + T[5] * s[1] + T[9] * s[2], 
					T[2] * s[0] + T[6] * s[1] + T[10] * s[2] };
	_w[0] = tmp[0];
	_w[1] = tmp[1];
	_w[2] = tmp[2];
	_w[3] = T[13] * tmp[2] - T[14] * tmp[1] + T[0] * s[3] + T[4] * s[4] + T[8] * s[5];
	_w[4] = T[14] * tmp[0] - T[12] * tmp[2] + T[1] * s[3] + T[5] * s[4] + T[9] * s[5];
	_w[5] = T[12] * tmp[1] - T[13] * tmp[0] + T[2] * s[3] + T[6] * s[4] + T[10] * s[5];
}

inline void se3::set_ad(const se3 &s1, const se3 &s2)
{
	_w[0] = s1._w[1] * s2._w[2] - s1._w[2] * s2._w[1];
	_w[1] = s1._w[2] * s2._w[0] - s1._w[0] * s2._w[2];
	_w[2] = s1._w[0] * s2._w[1] - s1._w[1] * s2._w[0];
	_w[3] = s1._w[1] * s2._w[5] - s1._w[2] * s2._w[4] - s2._w[1] * s1._w[5] + s2._w[2] * s1._w[4];
	_w[4] = s1._w[2] * s2._w[3] - s1._w[0] * s2._w[5] - s2._w[2] * s1._w[3] + s2._w[0] * s1._w[5];
	_w[5] = s1._w[0] * s2._w[4] - s1._w[1] * s2._w[3] - s2._w[0] * s1._w[4] + s2._w[1] * s1._w[3];
}

inline void se3::set_ad(const se3 &s1, const gReal *s2)
{
	if ( s2 == NULL ) { SetZero(); return; }

	_w[0] = s1._w[1] * s2[2] - s1._w[2] * s2[1];
	_w[1] = s1._w[2] * s2[0] - s1._w[0] * s2[2];
	_w[2] = s1._w[0] * s2[1] - s1._w[1] * s2[0];
	_w[3] = s1._w[1] * s2[5] - s1._w[2] * s2[4] - s2[1] * s1._w[5] + s2[2] * s1._w[4];
	_w[4] = s1._w[2] * s2[3] - s1._w[0] * s2[5] - s2[2] * s1._w[3] + s2[0] * s1._w[5];
	_w[5] = s1._w[0] * s2[4] - s1._w[1] * s2[3] - s2[0] * s1._w[4] + s2[1] * s1._w[3];
}

inline void se3::set_ad(const gReal *s1, const gReal *s2)
{
	if ( s1 == NULL || s2 == NULL ) { SetZero(); return; }

	_w[0] = s1[1] * s2[2] - s1[2] * s2[1];
	_w[1] = s1[2] * s2[0] - s1[0] * s2[2];
	_w[2] = s1[0] * s2[1] - s1[1] * s2[0];
	_w[3] = s1[1] * s2[5] - s1[2] * s2[4] - s2[1] * s1[5] + s2[2] * s1[4];
	_w[4] = s1[2] * s2[3] - s1[0] * s2[5] - s2[2] * s1[3] + s2[0] * s1[5];
	_w[5] = s1[0] * s2[4] - s1[1] * s2[3] - s2[0] * s1[4] + s2[1] * s1[3];
}

inline void se3::add_Ad(const SE3 &T, const se3 &s)
{
	gReal tmp[3] = { T._T[0] * s._w[0] + T._T[4] * s._w[1] + T._T[8] * s._w[2], 
					T._T[1] * s._w[0] + T._T[5] * s._w[1] + T._T[9] * s._w[2], 
					T._T[2] * s._w[0] + T._T[6] * s._w[1] + T._T[10] * s._w[2] };
	_w[0] += tmp[0];
	_w[1] += tmp[1];
	_w[2] += tmp[2];
	_w[3] += T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * s._w[3] + T._T[4] * s._w[4] + T._T[8] * s._w[5];
	_w[4] += T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * s._w[3] + T._T[5] * s._w[4] + T._T[9] * s._w[5];
	_w[5] += T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * s._w[3] + T._T[6] * s._w[4] + T._T[10] * s._w[5];
}

inline void se3::add_Ad(const SE3 &T, const gReal *s)
{
	if ( s == NULL ) { return; }

	gReal tmp[3] = { T._T[0] * s[0] + T._T[4] * s[1] + T._T[8] * s[2], 
					T._T[1] * s[0] + T._T[5] * s[1] + T._T[9] * s[2], 
					T._T[2] * s[0] + T._T[6] * s[1] + T._T[10] * s[2] };
	_w[0] += tmp[0];
	_w[1] += tmp[1];
	_w[2] += tmp[2];
	_w[3] += T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * s[3] + T._T[4] * s[4] + T._T[8] * s[5];
	_w[4] += T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * s[3] + T._T[5] * s[4] + T._T[9] * s[5];
	_w[5] += T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * s[3] + T._T[6] * s[4] + T._T[10] * s[5];
}

inline void se3::add_Ad(const gReal *T, const gReal *s)
{
	if ( T == NULL || s == NULL ) { return; }

	gReal tmp[3] = { T[0] * s[0] + T[4] * s[1] + T[8] * s[2], 
					T[1] * s[0] + T[5] * s[1] + T[9] * s[2], 
					T[2] * s[0] + T[6] * s[1] + T[10] * s[2] };
	_w[0] += tmp[0];
	_w[1] += tmp[1];
	_w[2] += tmp[2];
	_w[3] += T[13] * tmp[2] - T[14] * tmp[1] + T[0] * s[3] + T[4] * s[4] + T[8] * s[5];
	_w[4] += T[14] * tmp[0] - T[12] * tmp[2] + T[1] * s[3] + T[5] * s[4] + T[9] * s[5];
	_w[5] += T[12] * tmp[1] - T[13] * tmp[0] + T[2] * s[3] + T[6] * s[4] + T[10] * s[5];
}

inline void se3::add_ad(const se3 &s1, const se3 &s2)
{
	_w[0] += s1._w[1] * s2._w[2] - s1._w[2] * s2._w[1];
	_w[1] += s1._w[2] * s2._w[0] - s1._w[0] * s2._w[2];
	_w[2] += s1._w[0] * s2._w[1] - s1._w[1] * s2._w[0];
	_w[3] += s1._w[1] * s2._w[5] - s1._w[2] * s2._w[4] - s2._w[1] * s1._w[5] + s2._w[2] * s1._w[4];
	_w[4] += s1._w[2] * s2._w[3] - s1._w[0] * s2._w[5] - s2._w[2] * s1._w[3] + s2._w[0] * s1._w[5];
	_w[5] += s1._w[0] * s2._w[4] - s1._w[1] * s2._w[3] - s2._w[0] * s1._w[4] + s2._w[1] * s1._w[3];
}

inline void se3::add_ad(const se3 &s1, const gReal *s2)
{
	if ( s2 == NULL ) { return; }

	_w[0] += s1._w[1] * s2[2] - s1._w[2] * s2[1];
	_w[1] += s1._w[2] * s2[0] - s1._w[0] * s2[2];
	_w[2] += s1._w[0] * s2[1] - s1._w[1] * s2[0];
	_w[3] += s1._w[1] * s2[5] - s1._w[2] * s2[4] - s2[1] * s1._w[5] + s2[2] * s1._w[4];
	_w[4] += s1._w[2] * s2[3] - s1._w[0] * s2[5] - s2[2] * s1._w[3] + s2[0] * s1._w[5];
	_w[5] += s1._w[0] * s2[4] - s1._w[1] * s2[3] - s2[0] * s1._w[4] + s2[1] * s1._w[3];
}

inline void se3::add_ad(const gReal *s1, const gReal *s2)
{
	if ( s1 == NULL || s2 == NULL ) { return; }

	_w[0] += s1[1] * s2[2] - s1[2] * s2[1];
	_w[1] += s1[2] * s2[0] - s1[0] * s2[2];
	_w[2] += s1[0] * s2[1] - s1[1] * s2[0];
	_w[3] += s1[1] * s2[5] - s1[2] * s2[4] - s2[1] * s1[5] + s2[2] * s1[4];
	_w[4] += s1[2] * s2[3] - s1[0] * s2[5] - s2[2] * s1[3] + s2[0] * s1[5];
	_w[5] += s1[0] * s2[4] - s1[1] * s2[3] - s2[0] * s1[4] + s2[1] * s1[3];
}


inline SE3 Exp(const se3 &S)
{
	gReal theta = sqrt(S._w[0] * S._w[0] + S._w[1] * S._w[1] + S._w[2] * S._w[2]), itheta = (gReal)1.0 / theta;
	gReal s[6] = { itheta * S._w[0], itheta * S._w[1], itheta * S._w[2], itheta * S._w[3], itheta * S._w[4], itheta * S._w[5] };
	gReal st = sin(theta), ct = cos(theta), vt = (gReal)1.0 - ct, ut = (theta - st) * (s[0] * s[3] + s[1] * s[4] + s[2] * s[5]), t0 = s[2] * st, t1 = s[1] * st, t2 = s[0] * st;

	if ( fabs(theta) < _SMALL_LIE ) return SE3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, S._w[3], S._w[4], S._w[5]);
	
	return SE3( s[0] * s[0] * vt + ct,
				s[0] * s[1] * vt + t0,
				s[0] * s[2] * vt - t1,
				s[0] * s[1] * vt - t0,
				s[1] * s[1] * vt + ct,
				s[1] * s[2] * vt + t2,
				s[0] * s[2] * vt + t1,
				s[1] * s[2] * vt - t2,
				s[2] * s[2] * vt + ct,
				st * s[3] + ut * s[0] + vt * (s[1] * s[5] - s[2] * s[4]),
				st * s[4] + ut * s[1] + vt * (s[2] * s[3] - s[0] * s[5]),
				st * s[5] + ut * s[2] + vt * (s[0] * s[4] - s[1] * s[3]) );
}

inline SE3 Exp(const se3 &s, gReal theta)
{
	gReal mag = s._w[0] * s._w[0] + s._w[1] * s._w[1] + s._w[2] * s._w[2];

	if ( fabs(theta) < _SMALL_LIE || mag < _TINY_LIE ) return SE3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, theta * s._w[3], theta * s._w[4], theta * s._w[5]);
	else if ( fabs(mag - 1.0) > _TINY_LIE ) return Exp(theta * s);

	gReal st = sin(theta), ct = cos(theta), vt = (gReal)1.0 - ct, ut = (theta - st) * (s._w[0] * s._w[3] + s._w[1] * s._w[4] + s._w[2] * s._w[5]), t0 = s._w[2] * st, t1 = s._w[1] * st, t2 = s._w[0] * st;

	return SE3( s._w[0] * s._w[0] * vt + ct,
				s._w[0] * s._w[1] * vt + t0,
				s._w[0] * s._w[2] * vt - t1,
				s._w[0] * s._w[1] * vt - t0,
				s._w[1] * s._w[1] * vt + ct,
				s._w[1] * s._w[2] * vt + t2,
				s._w[0] * s._w[2] * vt + t1,
				s._w[1] * s._w[2] * vt - t2,
				s._w[2] * s._w[2] * vt + ct,
				st * s._w[3] + ut * s._w[0] + vt * (s._w[1] * s._w[5] - s._w[2] * s._w[4]),
				st * s._w[4] + ut * s._w[1] + vt * (s._w[2] * s._w[3] - s._w[0] * s._w[5]),
				st * s._w[5] + ut * s._w[2] + vt * (s._w[0] * s._w[4] - s._w[1] * s._w[3]) );
}

inline se3 Log(const SE3 &T)
{
	gReal d = (gReal)0.5 * (T._T[0] + T._T[5] + T._T[10] - (gReal)1.0);
	if ( d > 1.0 ) { d = (gReal)1.0; }
	if ( d < -1.0 ) { d = (gReal)-1.0; }
	gReal theta = acos(d);
	if ( fabs(theta) < _SMALL_LIE ) return se3(0.0, 0.0, 0.0, T._T[12], T._T[13], T._T[14]);
	gReal cof = theta / ((gReal)2.0 * sin(theta));
	gReal x = cof * (T._T[6] - T._T[9]), y = cof * (T._T[8] - T._T[2]), z = cof * (T._T[1] - T._T[4]);
	theta = sqrt(x * x + y * y + z * z);
	cof = ((gReal)2.0 * sin(theta) - theta * ((gReal)1.0 + cos(theta))) / ((gReal)2.0 * theta * theta * sin(theta));
	return se3( x, y, z,
				((gReal)1.0 - cof * (y * y + z * z)) * T._T[12] + ((gReal)0.5 * z + cof * x * y) * T._T[13] + (cof * x * z - (gReal)0.5 * y) * T._T[14],
				((gReal)1.0 - cof * (x * x + z * z)) * T._T[13] + ((gReal)0.5 * x + cof * z * y) * T._T[14] + (cof * x * y - (gReal)0.5 * z) * T._T[12],
				((gReal)1.0 - cof * (y * y + x * x)) * T._T[14] + ((gReal)0.5 * y + cof * x * z) * T._T[12] + (cof * y * z - (gReal)0.5 * x) * T._T[13]);
}

// re = T * s * Inv(T)
inline se3 Ad(const SE3 &T, const se3 &s)
{
	gReal tmp[3] = { T._T[0] * s._w[0] + T._T[4] * s._w[1] + T._T[8] * s._w[2], 
					T._T[1] * s._w[0] + T._T[5] * s._w[1] + T._T[9] * s._w[2], 
					T._T[2] * s._w[0] + T._T[6] * s._w[1] + T._T[10] * s._w[2] };

	return se3(	tmp[0], tmp[1], tmp[2],
				T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * s._w[3] + T._T[4] * s._w[4] + T._T[8] * s._w[5],
				T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * s._w[3] + T._T[5] * s._w[4] + T._T[9] * s._w[5],
				T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * s._w[3] + T._T[6] * s._w[4] + T._T[10] * s._w[5]);
}

//// re = Inv(T) * s * T
//inline se3 InvAd(const SE3 &T, const se3 &s)
//{
//	gReal tmp[3] = { s._w[3] + s._w[1] * T._T[14] - s._w[2] * T._T[13], 
//					s._w[4] + s._w[2] * T._T[12] - s._w[0] * T._T[14], 
//					s._w[5] + s._w[0] * T._T[13] - s._w[1] * T._T[12] };
//
//	return se3(	T._T[0] * s._w[0] + T._T[1] * s._w[1] + T._T[2] * s._w[2],
//				T._T[4] * s._w[0] + T._T[5] * s._w[1] + T._T[6] * s._w[2],
//				T._T[8] * s._w[0] + T._T[9] * s._w[1] + T._T[10] * s._w[2],
//				T._T[0] * tmp[0] + T._T[1] * tmp[1] + T._T[2] * tmp[2],
//				T._T[4] * tmp[0] + T._T[5] * tmp[1] + T._T[6] * tmp[2],
//				T._T[8] * tmp[0] + T._T[9] * tmp[1] + T._T[10] * tmp[2]	);
//}

inline se3 Ad(const Vec3 &p, const se3 &s)
{
	return se3(	s._w[0], 
				s._w[1], 
				s._w[2], 
				p._v[1] * s._w[2] - p._v[2] * s._w[1] + s._w[3], 
				p._v[2] * s._w[0] - p._v[0] * s._w[2] + s._w[4], 
				p._v[0] * s._w[1] - p._v[1] * s._w[0] + s._w[5]	);
}

inline se3 ad(const se3 &s1, const se3 &s2)
{
	return se3(	s1._w[1] * s2._w[2] - s1._w[2] * s2._w[1],
				s1._w[2] * s2._w[0] - s1._w[0] * s2._w[2],
				s1._w[0] * s2._w[1] - s1._w[1] * s2._w[0],
				s1._w[1] * s2._w[5] - s1._w[2] * s2._w[4] - s2._w[1] * s1._w[5] + s2._w[2] * s1._w[4],
				s1._w[2] * s2._w[3] - s1._w[0] * s2._w[5] - s2._w[2] * s1._w[3] + s2._w[0] * s1._w[5],
				s1._w[0] * s2._w[4] - s1._w[1] * s2._w[3] - s2._w[0] * s1._w[4] + s2._w[1] * s1._w[3] );
}

inline se3 ad(const gReal *s1, const gReal *s2)
{
	if ( s1 == NULL || s2 == NULL ) { return se3(0.0); }

	return se3(	s1[1] * s2[2] - s1[2] * s2[1],
				s1[2] * s2[0] - s1[0] * s2[2],
				s1[0] * s2[1] - s1[1] * s2[0],
				s1[1] * s2[5] - s1[2] * s2[4] - s2[1] * s1[5] + s2[2] * s1[4],
				s1[2] * s2[3] - s1[0] * s2[5] - s2[2] * s1[3] + s2[0] * s1[5],
				s1[0] * s2[4] - s1[1] * s2[3] - s2[0] * s1[4] + s2[1] * s1[3] );
}

inline dse3 dAd(const SE3 &T, const dse3 &t)
{
	gReal tmp[3] = { t._m[0] - T._T[13] * t._m[5] + T._T[14] * t._m[4], 
					t._m[1] - T._T[14] * t._m[3] + T._T[12] * t._m[5], 
					t._m[2] - T._T[12] * t._m[4] + T._T[13] * t._m[3] };

	return dse3(T._T[0] * tmp[0] + T._T[1] * tmp[1] + T._T[2] * tmp[2],
				T._T[4] * tmp[0] + T._T[5] * tmp[1] + T._T[6] * tmp[2],
				T._T[8] * tmp[0] + T._T[9] * tmp[1] + T._T[10] * tmp[2],
				T._T[0] * t._m[3] + T._T[1] * t._m[4] + T._T[2] * t._m[5],
				T._T[4] * t._m[3] + T._T[5] * t._m[4] + T._T[6] * t._m[5],
				T._T[8] * t._m[3] + T._T[9] * t._m[4] + T._T[10] * t._m[5]);
}

//inline dse3 InvdAd(const SE3 &T, const dse3 &t)
//{
//	gReal tmp[3] = { T._T[0] * t._m[3] + T._T[4] * t._m[4] + T._T[8] * t._m[5], 
//					T._T[1] * t._m[3] + T._T[5] * t._m[4] + T._T[9] * t._m[5], 
//					T._T[2] * t._m[3] + T._T[6] * t._m[4] + T._T[10] * t._m[5] };
//
//	return dse3(T._T[13] * tmp[2] - T._T[14] * tmp[1] + T._T[0] * t._m[0] + T._T[4] * t._m[1] + T._T[8] * t._m[2],
//				T._T[14] * tmp[0] - T._T[12] * tmp[2] + T._T[1] * t._m[0] + T._T[5] * t._m[1] + T._T[9] * t._m[2],
//				T._T[12] * tmp[1] - T._T[13] * tmp[0] + T._T[2] * t._m[0] + T._T[6] * t._m[1] + T._T[10] * t._m[2],
//				tmp[0], tmp[1], tmp[2]);
//}

inline dse3 dad(const se3 &s, const dse3 &t)
{
	return dse3(t._m[1] * s._w[2] - t._m[2] * s._w[1] + t._m[4] * s._w[5] - t._m[5] * s._w[4],
				t._m[2] * s._w[0] - t._m[0] * s._w[2] + t._m[5] * s._w[3] - t._m[3] * s._w[5],
				t._m[0] * s._w[1] - t._m[1] * s._w[0] + t._m[3] * s._w[4] - t._m[4] * s._w[3],
				t._m[4] * s._w[2] - t._m[5] * s._w[1],
				t._m[5] * s._w[0] - t._m[3] * s._w[2],
				t._m[3] * s._w[1] - t._m[4] * s._w[0]	);
}

// slightly efficient computation of dad(V, J * V)
inline dse3 dad(const se3 &V, const Inertia &J)
{
	gReal ww = V._w[0] * V._w[0] + V._w[1] * V._w[1] + V._w[2] * V._w[2];
	gReal wr = V._w[0] * J._r[0] + V._w[1] * J._r[1] + V._w[2] * J._r[2];
	gReal vr = V._w[3] * J._r[0] + V._w[4] * J._r[1] + V._w[5] * J._r[2];

	return dse3((J._I[3] * V._w[0] + J._I[1] * V._w[1] + J._I[5] * V._w[2]) * V._w[2] - (J._I[4] * V._w[0] + J._I[5] * V._w[1] + J._I[2] * V._w[2]) * V._w[1] - vr * V._w[0] + wr * V._w[3],
				(J._I[4] * V._w[0] + J._I[5] * V._w[1] + J._I[2] * V._w[2]) * V._w[0] - (J._I[0] * V._w[0] + J._I[3] * V._w[1] + J._I[4] * V._w[2]) * V._w[2] - vr * V._w[1] + wr * V._w[4],
				(J._I[0] * V._w[0] + J._I[3] * V._w[1] + J._I[4] * V._w[2]) * V._w[1] - (J._I[3] * V._w[0] + J._I[1] * V._w[1] + J._I[5] * V._w[2]) * V._w[0] - vr * V._w[2] + wr * V._w[5],
				ww * J._r[0] - wr * V._w[0] + J._m * (V._w[4] * V._w[2] - V._w[5] * V._w[1]),
				ww * J._r[1] - wr * V._w[1] + J._m * (V._w[5] * V._w[0] - V._w[3] * V._w[2]),
				ww * J._r[2] - wr * V._w[2] + J._m * (V._w[3] * V._w[1] - V._w[4] * V._w[0]));
}

inline SO3 SO3::operator * (const SO3 &R) const
{
	return SO3(	_R[0] * R._R[0] + _R[3] * R._R[1] + _R[6] * R._R[2], _R[1] * R._R[0] + _R[4] * R._R[1] + _R[7] * R._R[2], _R[2] * R._R[0] + _R[5] * R._R[1] + _R[8] * R._R[2],
				_R[0] * R._R[3] + _R[3] * R._R[4] + _R[6] * R._R[5], _R[1] * R._R[3] + _R[4] * R._R[4] + _R[7] * R._R[5], _R[2] * R._R[3] + _R[5] * R._R[4] + _R[8] * R._R[5],
				_R[0] * R._R[6] + _R[3] * R._R[7] + _R[6] * R._R[8], _R[1] * R._R[6] + _R[4] * R._R[7] + _R[7] * R._R[8], _R[2] * R._R[6] + _R[5] * R._R[7] + _R[8] * R._R[8]);
}

inline SO3 &SO3::operator *= (const SO3 &R)
{
	*this = *this * R;
	return *this;
}

inline SO3 Exp(gReal w0, gReal w1, gReal w2)
{
	gReal theta = sqrt(w0 * w0 + w1 * w1 + w2 * w2);
	if ( fabs(theta) < _TINY_LIE ) return SO3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
	
	w0 /= theta; w1 /= theta; w2 /= theta;
	gReal st = sin(theta), ct = cos(theta), vt = (gReal)1.0 - ct, t0 = w2 * st, t1 = w1 * st, t2 = w0 * st;

	return SO3( w0 * w0 * vt + ct, w0 * w1 * vt + t0, w0 * w2 * vt - t1,
				w0 * w1 * vt - t0, w1 * w1 * vt + ct, w1 * w2 * vt + t2,
				w0 * w2 * vt + t1, w1 * w2 * vt - t2, w2 * w2 * vt + ct );
}

inline SO3 Exp(const Vec3 &w)
{
	return Exp(w._v[0], w._v[1], w._v[2]);
}

inline Vec3 Log(const SO3 &R)
{
	gReal d = (gReal)0.5 * (R._R[0] + R._R[4] + R._R[8] - (gReal)1.0); 
	if ( d > 1.0 ) { d = (gReal)1.0; }
	if ( d < -1.0 ) { d = (gReal)-1.0; }
	gReal theta = acos(d), cof = theta / ((gReal)2.0 * sin(theta));

	if ( fabs(theta) < _SMALL_LIE ) return Vec3(0.0, 0.0, 0.0);

	return Vec3(cof * (R._R[5] - R._R[7]), cof * (R._R[6] - R._R[2]), cof * (R._R[1] - R._R[3]));
}

inline SO3 RotX(gReal theta)
{
	gReal c = cos(theta), s = sin(theta);
	return SO3(1, 0, 0, 0, c, s, 0, -s, c);
}

inline SO3 RotY(gReal theta)
{
	gReal c = cos(theta), s = sin(theta);
	return SO3(c, 0, -s, 0, 1, 0, s, 0, c);
}

inline SO3 RotZ(gReal theta)
{
	gReal c = cos(theta), s = sin(theta);
	return SO3(c, s, 0, -s, c, 0, 0, 0, 1);
}

inline SO3 EulerZYX(const Vec3 &x)
{
	gReal c0 = cos(x._v[0]), s0 = sin(x._v[0]), c1 = cos(x._v[1]), s1 = sin(x._v[1]), c2 = cos(x._v[2]), s2 = sin(x._v[2]);
	return SO3(c0 * c1, s0 * c1, -s1, c0 * s1 * s2 - s0 * c2, s0 * s1 * s2 + c0 * c2, c1 * s2, c0 * s1 * c2 + s0 * s2, s0 * s1 * c2 - c0 * s2, c1 * c2);
}

inline Vec3 iEulerZYX(const SO3 &R)
{
	return Vec3(atan2(R._R[1], R._R[0]), atan2(-R._R[2], sqrt(R._R[0] * R._R[0] + R._R[1] * R._R[1])), atan2(R._R[5], R._R[8]));
}

inline SO3 EulerZYZ(const Vec3 &x)
{
	gReal ca = cos(x._v[0]), sa = sin(x._v[0]), cb = cos(x._v[1]), sb = sin(x._v[1]), cg = cos(x._v[2]), sg = sin(x._v[2]);
	return SO3(ca * cb * cg - sa * sg, sa * cb * cg + ca * sg, -sb*cg,	-ca * cb * sg - sa * cg, ca * cg - sa * cb * sg, sb * sg, ca * sb, sa * sb, cb);
}
	
inline Vec3 iEulerZYZ(const SO3 &R)
{
	return Vec3(atan2(R._R[7], R._R[6]), atan2(sqrt(R._R[2] * R._R[2] + R._R[5] * R._R[5]), R._R[8]), atan2(R._R[5], -R._R[2]));
}

inline SO3 EulerZXY(const Vec3 &x)
{
	gReal c0 = cos(x._v[0]), s0 = sin(x._v[0]), c1 = cos(x._v[1]), s1 = sin(x._v[1]), c2 = cos(x._v[2]), s2 = sin(x._v[2]);
	return SO3(c0*c2 - s0*s1*s2, c2*s0 + c0*s1*s2, -c1*s2, -c1*s0, c0*c1, s1, c0*s2 + c2*s0*s1, s0*s2 - c0*c2*s1, c1*c2);
}

inline Vec3 iEulerZXY(const SO3 &R)
{
	return Vec3(atan2(-R._R[3], R._R[4]), atan2(R._R[5], sqrt(R._R[3]*R._R[3]+R._R[4]*R._R[4])), atan2(-R._R[2], R._R[8]));
}

inline SO3 EulerXYZ(const Vec3 &x)
{
	gReal c0 = cos(x._v[0]), s0 = sin(x._v[0]), c1 = cos(x._v[1]), s1 = sin(x._v[1]), c2 = cos(x._v[2]), s2 = sin(x._v[2]);
	return SO3(c1*c2, c0*s2 + c2*s0*s1, s0*s2 - c0*c2*s1, -c1*s2, c0*c2 - s0*s1*s2, c2*s0 + c0*s1*s2, s1, -c1*s0, c0*c1);
}

inline Vec3 iEulerXYZ(const SO3 &R)
{
	return Vec3(atan2(-R._R[7], R._R[8]), atan2(R._R[6], sqrt(R._R[7]*R._R[7]+R._R[8]*R._R[8])), atan2(-R._R[3], R._R[0]));
}

inline SO3 Quat(gReal *quat)
{
	// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm

	SO3 re;
	gReal sqw = quat[0]*quat[0];
	gReal sqx = quat[1]*quat[1];
	gReal sqy = quat[2]*quat[2];
	gReal sqz = quat[3]*quat[3];

	// invs (inverse square length) is only required if quaternion is not already normalized
	gReal invs = 1 / (sqx + sqy + sqz + sqw);
	re._R[0] = ( sqx - sqy - sqz + sqw)*invs ; // since sqw + sqx + sqy + sqz =1/invs*invs
	re._R[4] = (-sqx + sqy - sqz + sqw)*invs ;
	re._R[8] = (-sqx - sqy + sqz + sqw)*invs ;

	gReal tmp1 = quat[1]*quat[2];
	gReal tmp2 = quat[3]*quat[0];
	re._R[1] = (gReal)2.0 * (tmp1 + tmp2)*invs ;
	re._R[3] = (gReal)2.0 * (tmp1 - tmp2)*invs ;

	tmp1 = quat[1]*quat[3];
	tmp2 = quat[2]*quat[0];
	re._R[2] = (gReal)2.0 * (tmp1 - tmp2)*invs ;
	re._R[6] = (gReal)2.0 * (tmp1 + tmp2)*invs ;
	tmp1 = quat[2]*quat[3];
	tmp2 = quat[1]*quat[0];
	re._R[5] = (gReal)2.0 * (tmp1 + tmp2)*invs ;
	re._R[7] = (gReal)2.0 * (tmp1 - tmp2)*invs ; 

	return re;
}

inline void iQuat(const SO3 &R, gReal *quat)
{
	//// convert rotation matrix (R) to unit quaternion (quat)
	//quat[0] = (gReal)0.5 * sqrt((gReal)1.+R._R[0]+R._R[4]+R._R[8]);
	//quat[1] = (R._R[5]-R._R[7])/(4*quat[0]);
	//quat[2] = (R._R[6]-R._R[2])/(4*quat[0]);
	//quat[3] = (R._R[1]-R._R[3])/(4*quat[0]);

	// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
	gReal trace = R._R[0] + R._R[4] + R._R[8]; // I removed + 1.0f; see discussion with Ethan
	if( trace > 0 ) {// I changed M_EPSILON to 0
	gReal s = 0.5 / sqrt(trace+ 1.0);
	quat[0] = 0.25 / s;
	quat[1] = ( R._R[5] - R._R[7] ) * s;
	quat[2] = ( R._R[6] - R._R[2] ) * s;
	quat[3] = ( R._R[1] - R._R[3] ) * s;
	} else {
		if ( R._R[0] > R._R[4] && R._R[0] > R._R[8] ) {
			gReal s = 2.0 * sqrt( 1.0 + R._R[0] - R._R[4] - R._R[8]);
			quat[0] = (R._R[5] - R._R[7] ) / s;
			quat[1] = 0.25 * s;
			quat[2] = (R._R[3] + R._R[1] ) / s;
			quat[3] = (R._R[6] + R._R[2] ) / s;
		} else if (R._R[4] > R._R[8]) {
			gReal s = 2.0 * sqrt( 1.0 + R._R[4] - R._R[0] - R._R[8]);
			quat[0] = (R._R[6] - R._R[2] ) / s;
			quat[1] = (R._R[3] + R._R[1] ) / s;
			quat[2] = 0.25 * s;
			quat[3] = (R._R[7] + R._R[5] ) / s;
		} else {
			gReal s = 2.0 * sqrt( 1.0 + R._R[8] - R._R[0] - R._R[4] );
			quat[0] = (R._R[1] - R._R[3] ) / s;
			quat[1] = (R._R[6] + R._R[2] ) / s;
			quat[2] = (R._R[7] + R._R[5] ) / s;
			quat[3] = 0.25 * s;
		}
	}
}

inline bool isSO3(gReal *R, gReal eps)
{
	if (fabs(R[0]*R[1] + R[3]*R[4] + R[6]*R[7]) > eps) return false;
	if (fabs(R[0]*R[2] + R[3]*R[5] + R[6]*R[8]) > eps) return false;
	if (fabs(R[1]*R[2] + R[4]*R[5] + R[7]*R[8]) > eps) return false;
	if (fabs(R[0]*R[0] + R[3]*R[3] + R[6]*R[6] - 1) > eps) return false;
	if (fabs(R[1]*R[1] + R[4]*R[4] + R[7]*R[7] - 1) > eps) return false;
	if (fabs(R[2]*R[2] + R[5]*R[5] + R[8]*R[8] - 1) > eps) return false;
	if (fabs(R[0]*R[4]*R[8] + R[3]*R[7]*R[2] + R[6]*R[1]*R[5] - R[0]*R[7]*R[5] - R[3]*R[1]*R[8] - R[6]*R[4]*R[2] -1) > eps) return false;
	return true;
}

inline SO3 GetRotationWithZAxis(Vec3 axis)
{
	Vec3 z(0,0,1);
	axis.Normalize();

	if ( Norm(axis-z) < _TINY_LIE )		// if axis = z
		return SO3(); // identity
	else if ( Norm(axis+z) < _TINY_LIE )	// else if axis = -z
		return SO3(1,0,0,0,-1,0,0,0,-1);
	else
	{
		Vec3 m = Cross(z,axis);
		gReal theta = asin(Norm(m));
		if ( Inner(z,axis) < 0 ) theta = (gReal)_PI_LIE - theta;

		m.Normalize();
		return Exp(m*theta);
	}
}

inline SE3 &SE3::operator *= (const SE3 &T)
{
	gReal x0, x1, x2;
	
	_T[12] += _T[0] * T._T[12] + _T[4] * T._T[13] + _T[8] * T._T[14];
	_T[13] += _T[1] * T._T[12] + _T[5] * T._T[13] + _T[9] * T._T[14];
	_T[14] += _T[2] * T._T[12] + _T[6] * T._T[13] + _T[10] * T._T[14];
	
	x0 = _T[0] * T._T[0] + _T[4] * T._T[1] + _T[8] * T._T[2];
	x1 = _T[0] * T._T[4] + _T[4] * T._T[5] + _T[8] * T._T[6];
	x2 = _T[0] * T._T[8] + _T[4] * T._T[9] + _T[8] * T._T[10];
	_T[0] = x0; _T[4] = x1; _T[8] = x2;
	x0 = _T[1] * T._T[0] + _T[5] * T._T[1] + _T[9] * T._T[2];
	x1 = _T[1] * T._T[4] + _T[5] * T._T[5] + _T[9] * T._T[6];
	x2 = _T[1] * T._T[8] + _T[5] * T._T[9] + _T[9] * T._T[10];
	_T[1] = x0; _T[5] =x1; _T[9] = x2;
	x0 = _T[2] * T._T[0] + _T[6] * T._T[1] + _T[10] * T._T[2];
	x1 = _T[2] * T._T[4] + _T[6] * T._T[5] + _T[10] * T._T[6];
	x2 = _T[2] * T._T[8] + _T[6] * T._T[9] + _T[10] * T._T[10];
	_T[2] = x0; _T[6] = x1; _T[10] = x2;
	
	return  *this;
}

inline SE3 &SE3::operator /= (const SE3 &T)
{
	gReal tmp[9] = {	_T[0] * T._T[0] + _T[4] * T._T[4] + _T[8] * T._T[8],
					_T[1] * T._T[0] + _T[5] * T._T[4] + _T[9] * T._T[8],
					_T[2] * T._T[0] + _T[6] * T._T[4] + _T[10] * T._T[8],
					_T[0] * T._T[1] + _T[4] * T._T[5] + _T[8] * T._T[9],
					_T[1] * T._T[1] + _T[5] * T._T[5] + _T[9] * T._T[9],
					_T[2] * T._T[1] + _T[6] * T._T[5] + _T[10] * T._T[9],
					_T[0] * T._T[2] + _T[4] * T._T[6] + _T[8] * T._T[10],
					_T[1] * T._T[2] + _T[5] * T._T[6] + _T[9] * T._T[10],
					_T[2] * T._T[2] + _T[6] * T._T[6] + _T[10] * T._T[10] };
		
	_T[0] = tmp[0]; _T[1] = tmp[1]; _T[2] = tmp[2];
	_T[4] = tmp[3]; _T[5] = tmp[4]; _T[6] = tmp[5];
	_T[8] = tmp[6]; _T[9] = tmp[7]; _T[10] = tmp[8], 
	_T[12] -= tmp[0] * T._T[12] + tmp[3] * T._T[13] + tmp[6] * T._T[14];
	_T[13] -= tmp[1] * T._T[12] + tmp[4] * T._T[13] + tmp[7] * T._T[14];
	_T[14] -= tmp[2] * T._T[12] + tmp[5] * T._T[13] + tmp[8] * T._T[14];

	return *this;
}

inline SE3 &SE3::operator %= (const SE3 &T)
{
	gReal tmp[12] = { _T[0], _T[1], _T[2], _T[4], _T[5], _T[6], _T[8], _T[9], _T[10], T._T[12] - _T[12], T._T[13] - _T[13], T._T[14] - _T[14] };
	
	_T[0] = tmp[0] * T._T[0] + tmp[1] * T._T[1] + tmp[2] * T._T[2];
	_T[1] = tmp[3] * T._T[0] + tmp[4] * T._T[1] + tmp[5] * T._T[2];
	_T[2] = tmp[6] * T._T[0] + tmp[7] * T._T[1] + tmp[8] * T._T[2];
	_T[4] = tmp[0] * T._T[4] + tmp[1] * T._T[5] + tmp[2] * T._T[6];
	_T[5] = tmp[3] * T._T[4] + tmp[4] * T._T[5] + tmp[5] * T._T[6];
	_T[6] = tmp[6] * T._T[4] + tmp[7] * T._T[5] + tmp[8] * T._T[6];
	_T[8] = tmp[0] * T._T[8] + tmp[1] * T._T[9] + tmp[2] * T._T[10];
	_T[9] = tmp[3] * T._T[8] + tmp[4] * T._T[9] + tmp[5] * T._T[10];
	_T[10] = tmp[6] * T._T[8] + tmp[7] * T._T[9] + tmp[8] * T._T[10];
	_T[12] = tmp[0] * tmp[9] + tmp[1] * tmp[10] + tmp[2] * tmp[11];
	_T[13] = tmp[3] * tmp[9] + tmp[4] * tmp[10] + tmp[5] * tmp[11];
	_T[14] = tmp[6] * tmp[9] + tmp[7] * tmp[10] + tmp[8] * tmp[11];

	return *this;
}

inline SE3 SE3::operator * (const SE3 &T) const
{
	return SE3(	_T[0] * T._T[0] + _T[4] * T._T[1] + _T[8] * T._T[2],
				_T[1] * T._T[0] + _T[5] * T._T[1] + _T[9] * T._T[2],
				_T[2] * T._T[0] + _T[6] * T._T[1] + _T[10] * T._T[2],
				_T[0] * T._T[4] + _T[4] * T._T[5] + _T[8] * T._T[6],
				_T[1] * T._T[4] + _T[5] * T._T[5] + _T[9] * T._T[6],
				_T[2] * T._T[4] + _T[6] * T._T[5] + _T[10] * T._T[6],
				_T[0] * T._T[8] + _T[4] * T._T[9] + _T[8] * T._T[10],
				_T[1] * T._T[8] + _T[5] * T._T[9] + _T[9] * T._T[10],
				_T[2] * T._T[8] + _T[6] * T._T[9] + _T[10] * T._T[10],
				_T[12] + _T[0] * T._T[12] + _T[4] * T._T[13] + _T[8] * T._T[14],
				_T[13] + _T[1] * T._T[12] + _T[5] * T._T[13] + _T[9] * T._T[14],
				_T[14] + _T[2] * T._T[12] + _T[6] * T._T[13] + _T[10] * T._T[14] );
}

inline Vec3 SE3::operator * (const Vec3 &p) const
{
	return Vec3(_T[12] + _T[0] * p._v[0] + _T[4] * p._v[1] + _T[8] * p._v[2],
				_T[13] + _T[1] * p._v[0] + _T[5] * p._v[1] + _T[9] * p._v[2],
				_T[14] + _T[2] * p._v[0] + _T[6] * p._v[1] + _T[10] * p._v[2] );
}

inline SE3 SE3::operator / (const SE3 &T) const
{
	gReal tmp[9] = {	_T[0] * T._T[0] + _T[4] * T._T[4] + _T[8] * T._T[8],
					_T[1] * T._T[0] + _T[5] * T._T[4] + _T[9] * T._T[8],
					_T[2] * T._T[0] + _T[6] * T._T[4] + _T[10] * T._T[8],
					_T[0] * T._T[1] + _T[4] * T._T[5] + _T[8] * T._T[9],
					_T[1] * T._T[1] + _T[5] * T._T[5] + _T[9] * T._T[9],
					_T[2] * T._T[1] + _T[6] * T._T[5] + _T[10] * T._T[9],
					_T[0] * T._T[2] + _T[4] * T._T[6] + _T[8] * T._T[10],
					_T[1] * T._T[2] + _T[5] * T._T[6] + _T[9] * T._T[10],
					_T[2] * T._T[2] + _T[6] * T._T[6] + _T[10] * T._T[10] };
		
	return SE3(	tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8], 
				_T[12] - tmp[0] * T._T[12] - tmp[3] * T._T[13] - tmp[6] * T._T[14],
				_T[13] - tmp[1] * T._T[12] - tmp[4] * T._T[13] - tmp[7] * T._T[14],
				_T[14] - tmp[2] * T._T[12] - tmp[5] * T._T[13] - tmp[8] * T._T[14] );
}

inline SE3 SE3::operator % (const SE3 &T) const
{
	return SE3(	_T[0] * T._T[0] + _T[1] * T._T[1] + _T[2] * T._T[2],
				_T[4] * T._T[0] + _T[5] * T._T[1] + _T[6] * T._T[2],
				_T[8] * T._T[0] + _T[9] * T._T[1] + _T[10] * T._T[2],
				_T[0] * T._T[4] + _T[1] * T._T[5] + _T[2] * T._T[6],
				_T[4] * T._T[4] + _T[5] * T._T[5] + _T[6] * T._T[6],
				_T[8] * T._T[4] + _T[9] * T._T[5] + _T[10] * T._T[6],
				_T[0] * T._T[8] + _T[1] * T._T[9] + _T[2] * T._T[10],
				_T[4] * T._T[8] + _T[5] * T._T[9] + _T[6] * T._T[10],
				_T[8] * T._T[8] + _T[9] * T._T[9] + _T[10] * T._T[10],
				_T[0] * (T._T[12] - _T[12]) + _T[1] * (T._T[13] - _T[13]) + _T[2] * (T._T[14] - _T[14]),
				_T[4] * (T._T[12] - _T[12]) + _T[5] * (T._T[13] - _T[13]) + _T[6] * (T._T[14] - _T[14]),
				_T[8] * (T._T[12] - _T[12]) + _T[9] * (T._T[13] - _T[13]) + _T[10] * (T._T[14] - _T[14]) );
}

inline void SE3::SetInvOf(const SE3 &T)
{
	_T[0] = T._T[0];	_T[4] = T._T[1];	_T[8] = T._T[2];
	_T[1] = T._T[4];	_T[5] = T._T[5];	_T[9] = T._T[6];
	_T[2] = T._T[8];	_T[6] = T._T[9];	_T[10] = T._T[10];

	_T[12] = -T._T[0] * T._T[12] - T._T[1] * T._T[13] - T._T[2] * T._T[14];
	_T[13] = -T._T[4] * T._T[12] - T._T[5] * T._T[13] - T._T[6] * T._T[14];
	_T[14] = -T._T[8] * T._T[12] - T._T[9] * T._T[13] - T._T[10] * T._T[14];
}

inline SE3 &SE3::SetRotation(const SO3 &R)
{
	_T[0] = R._R[0]; _T[4] = R._R[3]; _T[8] = R._R[6];
	_T[1] = R._R[1]; _T[5] = R._R[4]; _T[9] = R._R[7];
	_T[2] = R._R[2]; _T[6] = R._R[5]; _T[10] = R._R[8];
	return  *this;
}

inline SE3 &SE3::SetPosition(const Vec3 &Pos)
{
	_T[12] = Pos._v[0]; _T[13] = Pos._v[1]; _T[14] = Pos._v[2];
	return  *this;
}

inline SE3 &SE3::Translate(const Vec3 &Pos)
{
	_T[12] += Pos._v[0]; _T[13] += Pos._v[1]; _T[14] += Pos._v[2];
	return  *this;
}

inline SE3 &SE3::Rotate(const SO3 &R)
{
	gReal r1, r2, r3;
	r1 = _T[0] * R._R[0] + _T[4] * R._R[1] + _T[8] * R._R[2];
	r2 = _T[0] * R._R[3] + _T[4] * R._R[4] + _T[8] * R._R[5];
	r3 = _T[0] * R._R[6] + _T[4] * R._R[7] + _T[8] * R._R[8];
	_T[0] = r1; _T[4] = r2; _T[8] = r3;
	r1 = _T[1] * R._R[0] + _T[5] * R._R[1] + _T[9] * R._R[2];
	r2 = _T[1] * R._R[3] + _T[5] * R._R[4] + _T[9] * R._R[5];
	r3 = _T[1] * R._R[6] + _T[5] * R._R[7] + _T[9] * R._R[8];
	_T[1] = r1; _T[5] = r2; _T[9] = r3;
	r1 = _T[2] * R._R[0] + _T[6] * R._R[1] + _T[10] * R._R[2];
	r2 = _T[2] * R._R[3] + _T[6] * R._R[4] + _T[10] * R._R[5];
	r3 = _T[2] * R._R[6] + _T[6] * R._R[7] + _T[10] * R._R[8];
	_T[2] = r1; _T[6] = r2; _T[10] = r3;
	return  *this;
}

inline SE3 Inv(const SE3 &T)
{
	return SE3(	T._T[0], T._T[4], T._T[8], T._T[1], T._T[5], T._T[9], T._T[2], T._T[6], T._T[10],
				-T._T[0] * T._T[12] - T._T[1] * T._T[13] - T._T[2] * T._T[14],
				-T._T[4] * T._T[12] - T._T[5] * T._T[13] - T._T[6] * T._T[14],
				-T._T[8] * T._T[12] - T._T[9] * T._T[13] - T._T[10] * T._T[14]);
}


inline Inertia::Inertia(gReal mass, gReal Ixx, gReal Iyy, gReal Izz)
{
	_m = mass;
	_I[0] = Ixx; _I[1] = Iyy; _I[2] = Izz;
	_I[3] = _I[4] = _I[5] = _r[0] = _r[1] = _r[2] = 0.0;
}
	
inline Inertia Inertia::Transform(const SE3 &T) const
{
	// operation count: multiplication = 101, addition = 48, subtraction = 18
	Inertia J(_m);
	
	J._I[0] = _I[0] + _m * T._T[14] * T._T[14] + _m * T._T[13] * T._T[13] - (gReal)2.0 * _r[2] * T._T[14] - (gReal)2.0 * _r[1] * T._T[13];
	J._I[3] = _I[3] + T._T[13] * _r[0] + T._T[12] * _r[1] - _m * T._T[13] * T._T[12];
	J._I[4] = _I[4] + T._T[14] * _r[0] + T._T[12] * _r[2] - _m * T._T[14] * T._T[12];
	J._I[1] = _I[1] + _m * T._T[14] * T._T[14] + _m * T._T[12] * T._T[12] - (gReal)2.0 * _r[2] * T._T[14] - (gReal)2.0 * _r[0] * T._T[12];
	J._I[5] = _I[5] + T._T[14] * _r[1] + T._T[13] * _r[2] - _m * T._T[14] * T._T[13];
	J._I[2] = _I[2] + _m * T._T[13] * T._T[13] + _m * T._T[12] * T._T[12] - (gReal)2.0 * _r[1] * T._T[13] - (gReal)2.0 * _r[0] * T._T[12];

	// _tmp = Transpose of Rotation Part of T * J._I
	gReal tmp[9] = {	T._T[0] * J._I[0] + T._T[1] * J._I[3] + T._T[2] * J._I[4],
					T._T[4] * J._I[0] + T._T[5] * J._I[3] + T._T[6] * J._I[4],
					T._T[8] * J._I[0] + T._T[9] * J._I[3] + T._T[10] * J._I[4],
					T._T[0] * J._I[3] + T._T[1] * J._I[1] + T._T[2] * J._I[5],
					T._T[4] * J._I[3] + T._T[5] * J._I[1] + T._T[6] * J._I[5],
					T._T[8] * J._I[3] + T._T[9] * J._I[1] + T._T[10] * J._I[5],
					T._T[0] * J._I[4] + T._T[1] * J._I[5] + T._T[2] * J._I[2],
					T._T[4] * J._I[4] + T._T[5] * J._I[5] + T._T[6] * J._I[2],
					T._T[8] * J._I[4] + T._T[9] * J._I[5] + T._T[10] * J._I[2] };
	// J._I = tmp * Rotation Part of T
	J._I[0] = tmp[0] * T._T[0] + tmp[3] * T._T[1] + tmp[6] * T._T[2];
	J._I[3] = tmp[1] * T._T[0] + tmp[4] * T._T[1] + tmp[7] * T._T[2];
	J._I[4] = tmp[2] * T._T[0] + tmp[5] * T._T[1] + tmp[8] * T._T[2];
	J._I[1] = tmp[1] * T._T[4] + tmp[4] * T._T[5] + tmp[7] * T._T[6];
	J._I[5] = tmp[2] * T._T[4] + tmp[5] * T._T[5] + tmp[8] * T._T[6];
	J._I[2] = tmp[2] * T._T[8] + tmp[5] * T._T[9] + tmp[8] * T._T[10];

	J._r[0] = T._T[0] * (_r[0] - _m * T._T[12]) + T._T[1] * (_r[1] - _m * T._T[13]) + T._T[2] * (_r[2] - _m * T._T[14]);
	J._r[1] = T._T[4] * (_r[0] - _m * T._T[12]) + T._T[5] * (_r[1] - _m * T._T[13]) + T._T[6] * (_r[2] - _m * T._T[14]);
	J._r[2] = T._T[8] * (_r[0] - _m * T._T[12]) + T._T[9] * (_r[1] - _m * T._T[13]) + T._T[10] * (_r[2] - _m * T._T[14]);
	
	return J;
}

inline dse3 Inertia::operator * (const se3 &s) const
{
	return dse3(_I[0] * s._w[0] + _I[3] * s._w[1] + _I[4] * s._w[2] + _r[1] * s._w[5] - _r[2] * s._w[4],
				_I[3] * s._w[0] + _I[1] * s._w[1] + _I[5] * s._w[2] + _r[2] * s._w[3] - _r[0] * s._w[5],
				_I[4] * s._w[0] + _I[5] * s._w[1] + _I[2] * s._w[2] + _r[0] * s._w[4] - _r[1] * s._w[3],
				s._w[1] * _r[2] - s._w[2] * _r[1] + _m * s._w[3],
				s._w[2] * _r[0] - s._w[0] * _r[2] + _m * s._w[4],
				s._w[0] * _r[1] - s._w[1] * _r[0] + _m * s._w[5]);
}

inline dse3 Inertia::operator * (const gReal *s) const
{
	if ( s == NULL ) { return dse3(0.0); }

	return dse3(_I[0] * s[0] + _I[3] * s[1] + _I[4] * s[2] + _r[1] * s[5] - _r[2] * s[4],
				_I[3] * s[0] + _I[1] * s[1] + _I[5] * s[2] + _r[2] * s[3] - _r[0] * s[5],
				_I[4] * s[0] + _I[5] * s[1] + _I[2] * s[2] + _r[0] * s[4] - _r[1] * s[3],
				s[1] * _r[2] - s[2] * _r[1] + _m * s[3],
				s[2] * _r[0] - s[0] * _r[2] + _m * s[4],
				s[0] * _r[1] - s[1] * _r[0] + _m * s[5]);
}

inline void Inertia::ToArray(gReal J[]) const
{
	J[0]  =  _I[0];		J[6]  =  _I[3];		J[12] =  _I[4];		J[18] =    0.0;		J[24] = -_r[2];		J[30] =  _r[1];
	J[1]  =  _I[3];		J[7]  =  _I[1];		J[13] =  _I[5];		J[19] =  _r[2];		J[25] =    0.0;		J[31] = -_r[0];
	J[2]  =  _I[4];		J[8]  =  _I[5];		J[14] =  _I[2];		J[20] = -_r[1];		J[26] =  _r[0];		J[32] =    0.0;
	J[3]  =    0.0;		J[9]  =  _r[2];		J[15] = -_r[1];		J[21] =     _m;		J[27] =    0.0;		J[33] =    0.0;
	J[4]  = -_r[2];		J[10] =    0.0;		J[16] =  _r[0];		J[22] =    0.0;		J[28] =     _m;		J[34] =    0.0;
	J[5]  =  _r[1];		J[11] = -_r[0];		J[17] =    0.0;		J[23] =    0.0;		J[29] =    0.0;		J[35] =     _m;
}

inline void Inertia::InvToArray(gReal invJ[]) const
{
	// J = [I, [r]; -[r], m*eye(3)]
	// invJ = Inv([I, [r]; -[r], m*eye(3)]) = [ K, -(1/m)K[r]; (1/m)[r]K, (1/m)eye(3)-(1/m^2)[r]K[r] ] where K = Inv(I+(1/m)[r]^2)

	gReal I2[9], K[9], B[9], C[9]; // I2 = I + (1/m)[r]^2, K = Inv(I2), B = (1/m)[r]K, C = (1/m)eye(3)-(1/m^2)[r]K[r]
	gReal I[9] = {_I[0], _I[3], _I[4], _I[3], _I[1], _I[5], _I[4], _I[5], _I[2]};
	gReal eye[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	gReal r[9] = {0, _r[2], -_r[1], -_r[2], 0, _r[0], _r[1], -_r[0], 0};
	gReal rr_m[9] = {-_r[2]*_r[2]-_r[1]*_r[1], _r[0]*_r[1], _r[0]*_r[2], _r[0]*_r[1], -_r[2]*_r[2]-_r[0]*_r[0], _r[1]*_r[2], _r[0]*_r[2], _r[1]*_r[2], -_r[1]*_r[1]-_r[0]*_r[0]};
	matMult(rr_m, (1./_m), 3, 3); // rr_m = (1/m)[r]^2
	matSet(I2, I, 3, 3); // I2 = I
	matAdd(I2, rr_m, 3, 3); // I2 += rr_m
	matSet_inv33sym(K, I2); // K = Inv(I2)
	matSet_multAB(B, r, K, 3, 3, 3, 3); // B = [r]K
	matMult(B, (1./_m), 3, 3); // B *= (1/m)
	matSet_multAB(C, B, r, 3, 3, 3, 3); // C = B[r] = (1/m)[r]K[r]
	matSubtract(C, eye, 3, 3); // C -= eye --> C = (1/m)[r]K[r] - eye
	matMult(C, -(1./_m), 3, 3); // C *= -(1/m) --> C = (1/m)eye - (1/m^2)[r]K[r]

	// invJ = [K, B^T; 
	//         B, C  ]
	invJ[0] = K[0]; invJ[6]  = K[3]; invJ[12] = K[6]; invJ[18] = B[0]; invJ[24] = B[1]; invJ[30] = B[2];
	invJ[1] = K[1]; invJ[7]  = K[4]; invJ[13] = K[7]; invJ[19] = B[3]; invJ[25] = B[4]; invJ[31] = B[5];
	invJ[2] = K[2]; invJ[8]  = K[5]; invJ[14] = K[8]; invJ[20] = B[6]; invJ[26] = B[7]; invJ[32] = B[8];
	invJ[3] = B[0]; invJ[9]  = B[3]; invJ[15] = B[6]; invJ[21] = C[0]; invJ[27] = C[3]; invJ[33] = C[6];
	invJ[4] = B[1]; invJ[10] = B[4]; invJ[16] = B[7]; invJ[22] = C[1]; invJ[28] = C[4]; invJ[34] = C[7];
	invJ[5] = B[2]; invJ[11] = B[5]; invJ[17] = B[8]; invJ[23] = C[2]; invJ[29] = C[5]; invJ[35] = C[8];
}

inline void set_Mult_Inertia_se3(gReal *re, const Inertia &I, const gReal *s)
{
	re[0] = I._I[0] * s[0] + I._I[3] * s[1] + I._I[4] * s[2] + I._r[1] * s[5] - I._r[2] * s[4];
	re[1] = I._I[3] * s[0] + I._I[1] * s[1] + I._I[5] * s[2] + I._r[2] * s[3] - I._r[0] * s[5];
	re[2] = I._I[4] * s[0] + I._I[5] * s[1] + I._I[2] * s[2] + I._r[0] * s[4] - I._r[1] * s[3];
	re[3] = s[1] * I._r[2] - s[2] * I._r[1] + I._m * s[3];
	re[4] = s[2] * I._r[0] - s[0] * I._r[2] + I._m * s[4];
	re[5] = s[0] * I._r[1] - s[1] * I._r[0] + I._m * s[5];
}

inline void set_Mult_Inertia_se3(gReal *re, const Inertia &I, const gReal *s, int num)
{
	for (int i=0; i<num; i++) {
		set_Mult_Inertia_se3(&re[6*i], I, &s[6*i]);
	}
}

inline AInertia::AInertia(gReal a0, gReal a1, gReal a2, gReal a3, gReal a4, gReal a5, gReal a6, gReal a7, gReal a8, gReal a9, gReal a10, gReal a11, gReal a12, gReal a13, gReal a14, gReal a15, gReal a16, gReal a17, gReal a18, gReal a19, gReal a20)
{
	_J[0] = a0;		_J[1] = a1;		_J[2] = a2; 	_J[3] = a3;		_J[4] = a4;		_J[5] = a5;
					_J[6] = a6;		_J[7] = a7;		_J[8] = a8;		_J[9] = a9;		_J[10] = a10;	
									_J[11] = a11;	_J[12] = a12;	_J[13] = a13;	_J[14] = a14;
													_J[15] = a15;	_J[16] = a16;	_J[17] = a17;
																	_J[18] = a18;	_J[19] = a19;	
																					_J[20] = a20;
}

inline AInertia::AInertia(const AInertia &J)
{
	_J[0] = J._J[0]; _J[1] = J._J[1]; _J[2] = J._J[2];   _J[3] = J._J[3];   _J[4] = J._J[4];   _J[5] = J._J[5]; 
					 _J[6] = J._J[6]; _J[7] = J._J[7];   _J[8] = J._J[8];   _J[9] = J._J[9];   _J[10] = J._J[10]; 
					 				  _J[11] = J._J[11]; _J[12] = J._J[12]; _J[13] = J._J[13]; _J[14] = J._J[14];
														 _J[15] = J._J[15]; _J[16] = J._J[16]; _J[17] = J._J[17]; 
														 					_J[18] = J._J[18]; _J[19] = J._J[19]; 
														 									   _J[20] = J._J[20];
}

inline AInertia::AInertia(const Inertia &J)
{
	_J[0] = J._I[0];	_J[1] = J._I[3];	_J[2] = J._I[4];	_J[3] = 0.0;		_J[4] = -J._r[2];	_J[5] = J._r[1];	
						_J[6] = J._I[1];	_J[7] = J._I[5];	_J[8] = J._r[2];	_J[9] = 0.0;		_J[10] = -J._r[0];
											_J[11] = J._I[2];	_J[12] = -J._r[1];	_J[13] = J._r[0];	_J[14] = 0.0;
																_J[15] = J._m;		_J[16] = 0.0;		_J[17] = 0.0;		
																					_J[18] = J._m;		_J[19] = 0.0;
																										_J[20] = J._m;
}

inline AInertia::AInertia(const gReal *M)
{
	_J[0] = M[0]; _J[1] = M[ 6]; _J[2] = M[12];  _J[3] = M[18];  _J[4] = M[24];  _J[5] = M[30];
	              _J[6] = M[ 7]; _J[7] = M[13];  _J[8] = M[19];  _J[9] = M[25];  _J[10] = M[31];
	                             _J[11] = M[14]; _J[12] = M[20]; _J[13] = M[26]; _J[14] = M[32];
	                                             _J[15] = M[21]; _J[16] = M[27]; _J[17] = M[33];
	                                                             _J[18] = M[28]; _J[19] = M[34];
	                                                                             _J[20] = M[35];
}

inline AInertia AInertia::operator - (void) const
{
	return AInertia(-_J[0], -_J[1], -_J[2], -_J[3], -_J[4], -_J[5], -_J[6], -_J[7], -_J[8], -_J[9], -_J[10], -_J[11], -_J[12], -_J[13], -_J[14], -_J[15], -_J[16], -_J[17], -_J[18], -_J[19], -_J[20]);
}

inline dse3 AInertia::operator * (const se3 &a) const
{
	return dse3(_J[0] * a._w[0] + _J[1] * a._w[1] + _J[2] * a._w[2] + _J[3] * a._w[3] + _J[4] * a._w[4] + _J[5] * a._w[5],
				_J[1] * a._w[0] + _J[6] * a._w[1] + _J[7] * a._w[2] + _J[8] * a._w[3] + _J[9] * a._w[4] + _J[10] * a._w[5],
				_J[2] * a._w[0] + _J[7] * a._w[1] + _J[11] * a._w[2] + _J[12] * a._w[3] + _J[13] * a._w[4] + _J[14] * a._w[5],
				_J[3] * a._w[0] + _J[8] * a._w[1] + _J[12] * a._w[2] + _J[15] * a._w[3] + _J[16] * a._w[4] + _J[17] * a._w[5],
				_J[4] * a._w[0] + _J[9] * a._w[1] + _J[13] * a._w[2] + _J[16] * a._w[3] + _J[18] * a._w[4] + _J[19] * a._w[5],
				_J[5] * a._w[0] + _J[10] * a._w[1] + _J[14] * a._w[2] + _J[17] * a._w[3] + _J[19] * a._w[4] + _J[20] * a._w[5]);
}

inline AInertia AInertia::operator + (const AInertia &J) const
{
	return AInertia(_J[0] + J._J[0], _J[1] + J._J[1], _J[2] + J._J[2], _J[3] + J._J[3], _J[4] + J._J[4], _J[5] + J._J[5], _J[6] + J._J[6], _J[7] + J._J[7], _J[8] + J._J[8], _J[9] + J._J[9], _J[10] + J._J[10], _J[11] + J._J[11], _J[12] + J._J[12], _J[13] + J._J[13], _J[14] + J._J[14], _J[15] + J._J[15], _J[16] + J._J[16], _J[17] + J._J[17], _J[18] + J._J[18], _J[19] + J._J[19], _J[20] + J._J[20]);
}

inline AInertia AInertia::operator + (const Inertia &I) const
{
	return AInertia(	_J[0] + I._I[0],	_J[1] + I._I[3],	_J[2] + I._I[4],	_J[3],				_J[4] - I._r[2],	_J[5] + I._r[1],	
											_J[6] + I._I[1],	_J[7] + I._I[5],	_J[8] + I._r[2],	_J[9],				_J[10] - I._r[0],
																_J[11] + I._I[2],	_J[12] - I._r[1],	_J[13] + I._r[0],	_J[14],
																					_J[15] + I._m,		_J[16],				_J[17],		
																										_J[18] + I._m,		_J[19],
																															_J[20] + I._m	);
}

inline AInertia AInertia::operator - (const AInertia &J) const
{
	return AInertia(_J[0] - J._J[0], _J[1] - J._J[1], _J[2] - J._J[2], _J[3] - J._J[3], _J[4] - J._J[4], _J[5] - J._J[5], _J[6] - J._J[6], _J[7] - J._J[7], _J[8] - J._J[8], _J[9] - J._J[9], _J[10] - J._J[10], _J[11] - J._J[11], _J[12] - J._J[12], _J[13] - J._J[13], _J[14] - J._J[14], _J[15] - J._J[15], _J[16] - J._J[16], _J[17] - J._J[17], _J[18] - J._J[18], _J[19] - J._J[19], _J[20] - J._J[20]);
}

inline AInertia AInertia::operator - (const Inertia &J) const
{
	return AInertia(	_J[0] - J._I[0],	_J[1] - J._I[3],	_J[2] - J._I[4],	_J[3],				_J[4] + J._r[2],	_J[5] - J._r[1],	
											_J[6] - J._I[1],	_J[7] - J._I[5],	_J[8] - J._r[2],	_J[9],				_J[10] + J._r[0],
																_J[11] - J._I[2],	_J[12] + J._r[1],	_J[13] - J._r[0],	_J[14],
																					_J[15] - J._m,		_J[16],				_J[17],		
																										_J[18] - J._m,		_J[19],
																															_J[20] - J._m	);
}

inline AInertia &AInertia::operator = (const AInertia &J)
{
	_J[0] = J._J[0]; _J[1] = J._J[1]; _J[2] = J._J[2];   _J[3] = J._J[3];   _J[4] = J._J[4];   _J[5] = J._J[5]; 
					 _J[6] = J._J[6]; _J[7] = J._J[7];   _J[8] = J._J[8];   _J[9] = J._J[9];   _J[10] = J._J[10]; 
					 				  _J[11] = J._J[11]; _J[12] = J._J[12]; _J[13] = J._J[13]; _J[14] = J._J[14];
														 _J[15] = J._J[15]; _J[16] = J._J[16]; _J[17] = J._J[17]; 
														 					_J[18] = J._J[18]; _J[19] = J._J[19]; 
														 									   _J[20] = J._J[20];
	return *this;
}

inline AInertia &AInertia::operator = (const Inertia &J)
{
	_J[0] = J._I[0];	_J[1] = J._I[3];	_J[2] = J._I[4];	_J[3] = 0.0;		_J[4] = -J._r[2];	_J[5] = J._r[1];	
						_J[6] = J._I[1];	_J[7] = J._I[5];	_J[8] = J._r[2];	_J[9] = 0.0;		_J[10] = -J._r[0];
											_J[11] = J._I[2];	_J[12] = -J._r[1];	_J[13] = J._r[0];	_J[14] = 0.0;
																_J[15] = J._m;		_J[16] = 0.0;		_J[17] = 0.0;		
																					_J[18] = J._m;		_J[19] = 0.0;
																										_J[20] = J._m;
	return *this;
}

inline AInertia &AInertia::operator += (const gReal *M)
{
	if ( M == NULL ) { return *this; }

	_J[0] += M[0]; 	_J[1] += M[ 6];	 _J[2] += M[12]; 	_J[3] += M[18]; 	_J[4] += M[24]; 	_J[5] += M[30];
	               	_J[6] += M[ 7];	 _J[7] += M[13]; 	_J[8] += M[19]; 	_J[9] += M[25]; 	_J[10] += M[31];
	                               	_J[11] += M[14]; 	_J[12] += M[20]; 	_J[13] += M[26]; 	_J[14] += M[32];
	                                               		_J[15] += M[21]; 	_J[16] += M[27]; 	_J[17] += M[33];
	                                                               			_J[18] += M[28]; 	_J[19] += M[34];
	                                                                               				_J[20] += M[35];
	return *this;
}

inline AInertia &AInertia::operator -= (const gReal *M)
{
	if ( M == NULL ) { return *this; }

	_J[0] -= M[0]; 	_J[1] -= M[ 6];	 _J[2] -= M[12]; 	_J[3] -= M[18]; 	_J[4] -= M[24]; 	_J[5] -= M[30];
	               	_J[6] -= M[ 7];	 _J[7] -= M[13]; 	_J[8] -= M[19]; 	_J[9] -= M[25]; 	_J[10] -= M[31];
	                               	_J[11] -= M[14]; 	_J[12] -= M[20]; 	_J[13] -= M[26]; 	_J[14] -= M[32];
	                                               		_J[15] -= M[21]; 	_J[16] -= M[27]; 	_J[17] -= M[33];
	                                                               			_J[18] -= M[28]; 	_J[19] -= M[34];
	                                                                               				_J[20] -= M[35];
	return *this;
}

inline AInertia &AInertia::operator += (const AInertia &J)
{
	_J[0] += J._J[0]; 	_J[1] += J._J[1]; 	_J[2] += J._J[2];   _J[3] += J._J[3];   _J[4] += J._J[4];   _J[5] += J._J[5]; 
					 	_J[6] += J._J[6]; 	_J[7] += J._J[7];   _J[8] += J._J[8];   _J[9] += J._J[9];   _J[10] += J._J[10]; 
					 				  		_J[11] += J._J[11]; _J[12] += J._J[12]; _J[13] += J._J[13]; _J[14] += J._J[14];
														 		_J[15] += J._J[15]; _J[16] += J._J[16]; _J[17] += J._J[17]; 
														 							_J[18] += J._J[18]; _J[19] += J._J[19]; 
														 									   			_J[20] += J._J[20];
	return *this;
}

inline AInertia &AInertia::operator -= (const AInertia &J)
{
	_J[0] -= J._J[0]; 	_J[1] -= J._J[1]; 	_J[2] -= J._J[2];   _J[3] -= J._J[3];   _J[4] -= J._J[4];   _J[5] -= J._J[5]; 
					 	_J[6] -= J._J[6]; 	_J[7] -= J._J[7];   _J[8] -= J._J[8];   _J[9] -= J._J[9];   _J[10] -= J._J[10]; 
					 				  		_J[11] -= J._J[11]; _J[12] -= J._J[12]; _J[13] -= J._J[13]; _J[14] -= J._J[14];
														 		_J[15] -= J._J[15]; _J[16] -= J._J[16]; _J[17] -= J._J[17]; 
														 							_J[18] -= J._J[18]; _J[19] -= J._J[19]; 
														 									   			_J[20] -= J._J[20];
	return *this;
}

inline AInertia &AInertia::operator += (const Inertia &J)
{
	_J[0] += J._I[0];	_J[1] += J._I[3];	_J[2] += J._I[4];						_J[4] -= J._r[2];	_J[5] += J._r[1];	
						_J[6] += J._I[1];	_J[7] += J._I[5];	_J[8] += J._r[2];						_J[10] -= J._r[0];
											_J[11] += J._I[2];	_J[12] -= J._r[1];	_J[13] += J._r[0];	
																_J[15] += J._m;		
																					_J[18] += J._m;	
																										_J[20] += J._m;
	return *this;
}

inline AInertia &AInertia::operator -= (const Inertia &J)
{
	_J[0] -= J._I[0];	_J[1] -= J._I[3];	_J[2] -= J._I[4];						_J[4] += J._r[2];	_J[5] -= J._r[1];	
						_J[6] -= J._I[1];	_J[7] -= J._I[5];	_J[8] -= J._r[2];						_J[10] += J._r[0];
											_J[11] -= J._I[2];	_J[12] += J._r[1];	_J[13] -= J._r[0];	
																_J[15] -= J._m;		
																					_J[18] -= J._m;	
																										_J[20] -= J._m;
	return *this;
}

inline void AInertia::ToArray(gReal I[]) const
{
	I[0] = _J[0]; I[ 6] = _J[1]; I[12] = _J[2]; I[18] = _J[3]; I[24] = _J[4]; I[30] = _J[5];
	I[1] = _J[1]; I[ 7] = _J[6]; I[13] = _J[7]; I[19] = _J[8]; I[25] = _J[9]; I[31] = _J[10];
	I[2] = _J[2]; I[ 8] = _J[7]; I[14] = _J[11]; I[20] = _J[12]; I[26] = _J[13]; I[32] = _J[14];
	I[3] = _J[3]; I[ 9] = _J[8]; I[15] = _J[12]; I[21] = _J[15]; I[27] = _J[16]; I[33] = _J[17];
	I[4] = _J[4]; I[10] = _J[9]; I[16] = _J[13]; I[22] = _J[16]; I[28] = _J[18]; I[34] = _J[19];
	I[5] = _J[5]; I[11] = _J[10]; I[17] = _J[14]; I[23] = _J[17]; I[29] = _J[19]; I[35] = _J[20];
}

inline AInertia AInertia::Transform(const SE3 &T) const
{
	// operation count: multiplication = 186, addition = 117, subtract = 21

	gReal d0 = _J[3] + T._T[14] * _J[16] - T._T[13] * _J[17];
	gReal d1 = _J[8] - T._T[14] * _J[15] + T._T[12] * _J[17];
	gReal d2 = _J[12] + T._T[13] * _J[15] - T._T[12] * _J[16];
	gReal d3 = _J[4] + T._T[14] * _J[18] - T._T[13] * _J[19];
	gReal d4 = _J[9] - T._T[14] * _J[16] + T._T[12] * _J[19];
	gReal d5 = _J[13] + T._T[13] * _J[16] - T._T[12] * _J[18];
	gReal d6 = _J[5] + T._T[14] * _J[19] - T._T[13] * _J[20];
	gReal d7 = _J[10] - T._T[14] * _J[17] + T._T[12] * _J[20];
	gReal d8 = _J[14] + T._T[13] * _J[17] - T._T[12] * _J[19];
	gReal e0 = _J[0] + T._T[14] * _J[4] - T._T[13] * _J[5] + d3 * T._T[14] - d6 * T._T[13];
	gReal e3 = _J[1] + T._T[14] * _J[9] - T._T[13] * _J[10] - d0 * T._T[14] + d6 * T._T[12];
	gReal e4 = _J[6] - T._T[14] * _J[8] + T._T[12] * _J[10] - d1 * T._T[14] + d7 * T._T[12];
	gReal e6 = _J[2] + T._T[14] * _J[13] - T._T[13] * _J[14] + d0 * T._T[13] - d3 * T._T[12];
	gReal e7 = _J[7] - T._T[14] * _J[12] + T._T[12] * _J[14] + d1 * T._T[13] - d4 * T._T[12];
	gReal e8 = _J[11] + T._T[13] * _J[12] - T._T[12] * _J[13] + d2 * T._T[13] - d5 * T._T[12];
	gReal f0 = T._T[0] * e0 + T._T[1] * e3 + T._T[2] * e6;
	gReal f1 = T._T[0] * e3 + T._T[1] * e4 + T._T[2] * e7;
	gReal f2 = T._T[0] * e6 + T._T[1] * e7 + T._T[2] * e8;
	gReal f3 = T._T[0] * d0 + T._T[1] * d1 + T._T[2] * d2;
	gReal f4 = T._T[0] * d3 + T._T[1] * d4 + T._T[2] * d5;
	gReal f5 = T._T[0] * d6 + T._T[1] * d7 + T._T[2] * d8;
	gReal f6 = T._T[4] * e0 + T._T[5] * e3 + T._T[6] * e6;
	gReal f7 = T._T[4] * e3 + T._T[5] * e4 + T._T[6] * e7;
	gReal f8 = T._T[4] * e6 + T._T[5] * e7 + T._T[6] * e8;
	gReal g0 = T._T[4] * d0 + T._T[5] * d1 + T._T[6] * d2;
	gReal g1 = T._T[4] * d3 + T._T[5] * d4 + T._T[6] * d5;
	gReal g2 = T._T[4] * d6 + T._T[5] * d7 + T._T[6] * d8;
	gReal g3 = T._T[8] * d0 + T._T[9] * d1 + T._T[10] * d2;
	gReal g4 = T._T[8] * d3 + T._T[9] * d4 + T._T[10] * d5;
	gReal g5 = T._T[8] * d6 + T._T[9] * d7 + T._T[10] * d8;
	gReal h0 = T._T[0] * _J[15] + T._T[1] * _J[16] + T._T[2] * _J[17];
	gReal h1 = T._T[0] * _J[16] + T._T[1] * _J[18] + T._T[2] * _J[19];
	gReal h2 = T._T[0] * _J[17] + T._T[1] * _J[19] + T._T[2] * _J[20];
	gReal h3 = T._T[4] * _J[15] + T._T[5] * _J[16] + T._T[6] * _J[17];
	gReal h4 = T._T[4] * _J[16] + T._T[5] * _J[18] + T._T[6] * _J[19];
	gReal h5 = T._T[4] * _J[17] + T._T[5] * _J[19] + T._T[6] * _J[20];

	return AInertia(f0 * T._T[0] + f1 * T._T[1] + f2 * T._T[2],
					f0 * T._T[4] + f1 * T._T[5] + f2 * T._T[6],
					f0 * T._T[8] + f1 * T._T[9] + f2 * T._T[10],
					f3 * T._T[0] + f4 * T._T[1] + f5 * T._T[2],
					f3 * T._T[4] + f4 * T._T[5] + f5 * T._T[6],
					f3 * T._T[8] + f4 * T._T[9] + f5 * T._T[10],
					f6 * T._T[4] + f7 * T._T[5] + f8 * T._T[6],
					f6 * T._T[8] + f7 * T._T[9] + f8 * T._T[10],
					g0 * T._T[0] + g1 * T._T[1] + g2 * T._T[2],
					g0 * T._T[4] + g1 * T._T[5] + g2 * T._T[6],
					g0 * T._T[8] + g1 * T._T[9] + g2 * T._T[10],
					(T._T[8] * e0 + T._T[9] * e3 + T._T[10] * e6) * T._T[8] + (T._T[8] * e3 + T._T[9] * e4 + T._T[10] * e7) * T._T[9] + (T._T[8] * e6 + T._T[9] * e7 + T._T[10] * e8) * T._T[10],
					g3 * T._T[0] + g4 * T._T[1] + g5 * T._T[2],
					g3 * T._T[4] + g4 * T._T[5] + g5 * T._T[6],
					g3 * T._T[8] + g4 * T._T[9] + g5 * T._T[10],
					h0 * T._T[0] + h1 * T._T[1] + h2 * T._T[2],
					h0 * T._T[4] + h1 * T._T[5] + h2 * T._T[6],
					h0 * T._T[8] + h1 * T._T[9] + h2 * T._T[10],
					h3 * T._T[4] + h4 * T._T[5] + h5 * T._T[6],
					h3 * T._T[8] + h4 * T._T[9] + h5 * T._T[10],
					(T._T[8] * _J[15] + T._T[9] * _J[16] + T._T[10] * _J[17]) * T._T[8] + (T._T[8] * _J[16] + T._T[9] * _J[18] + T._T[10] * _J[19]) * T._T[9] + (T._T[8] * _J[17] + T._T[9] * _J[19] + T._T[10] * _J[20]) * T._T[10]);
}

inline void AInertia::AddTransform(const AInertia &J, const SE3 &T)
{	
	gReal d0 = J._J[3] + T._T[14] * J._J[16] - T._T[13] * J._J[17];
	gReal d1 = J._J[8] - T._T[14] * J._J[15] + T._T[12] * J._J[17];
	gReal d2 = J._J[12] + T._T[13] * J._J[15] - T._T[12] * J._J[16];
	gReal d3 = J._J[4] + T._T[14] * J._J[18] - T._T[13] * J._J[19];
	gReal d4 = J._J[9] - T._T[14] * J._J[16] + T._T[12] * J._J[19];
	gReal d5 = J._J[13] + T._T[13] * J._J[16] - T._T[12] * J._J[18];
	gReal d6 = J._J[5] + T._T[14] * J._J[19] - T._T[13] * J._J[20];
	gReal d7 = J._J[10] - T._T[14] * J._J[17] + T._T[12] * J._J[20];
	gReal d8 = J._J[14] + T._T[13] * J._J[17] - T._T[12] * J._J[19];
	gReal e0 = J._J[0] + T._T[14] * J._J[4] - T._T[13] * J._J[5] + d3 * T._T[14] - d6 * T._T[13];
	gReal e3 = J._J[1] + T._T[14] * J._J[9] - T._T[13] * J._J[10] - d0 * T._T[14] + d6 * T._T[12];
	gReal e4 = J._J[6] - T._T[14] * J._J[8] + T._T[12] * J._J[10] - d1 * T._T[14] + d7 * T._T[12];
	gReal e6 = J._J[2] + T._T[14] * J._J[13] - T._T[13] * J._J[14] + d0 * T._T[13] - d3 * T._T[12];
	gReal e7 = J._J[7] - T._T[14] * J._J[12] + T._T[12] * J._J[14] + d1 * T._T[13] - d4 * T._T[12];
	gReal e8 = J._J[11] + T._T[13] * J._J[12] - T._T[12] * J._J[13] + d2 * T._T[13] - d5 * T._T[12];
	gReal f0 = T._T[0] * e0 + T._T[1] * e3 + T._T[2] * e6;
	gReal f1 = T._T[0] * e3 + T._T[1] * e4 + T._T[2] * e7;
	gReal f2 = T._T[0] * e6 + T._T[1] * e7 + T._T[2] * e8;
	gReal f3 = T._T[0] * d0 + T._T[1] * d1 + T._T[2] * d2;
	gReal f4 = T._T[0] * d3 + T._T[1] * d4 + T._T[2] * d5;
	gReal f5 = T._T[0] * d6 + T._T[1] * d7 + T._T[2] * d8;
	gReal f6 = T._T[4] * e0 + T._T[5] * e3 + T._T[6] * e6;
	gReal f7 = T._T[4] * e3 + T._T[5] * e4 + T._T[6] * e7;
	gReal f8 = T._T[4] * e6 + T._T[5] * e7 + T._T[6] * e8;
	gReal g0 = T._T[4] * d0 + T._T[5] * d1 + T._T[6] * d2;
	gReal g1 = T._T[4] * d3 + T._T[5] * d4 + T._T[6] * d5;
	gReal g2 = T._T[4] * d6 + T._T[5] * d7 + T._T[6] * d8;
	gReal g3 = T._T[8] * d0 + T._T[9] * d1 + T._T[10] * d2;
	gReal g4 = T._T[8] * d3 + T._T[9] * d4 + T._T[10] * d5;
	gReal g5 = T._T[8] * d6 + T._T[9] * d7 + T._T[10] * d8;
	gReal h0 = T._T[0] * J._J[15] + T._T[1] * J._J[16] + T._T[2] * J._J[17];
	gReal h1 = T._T[0] * J._J[16] + T._T[1] * J._J[18] + T._T[2] * J._J[19];
	gReal h2 = T._T[0] * J._J[17] + T._T[1] * J._J[19] + T._T[2] * J._J[20];
	gReal h3 = T._T[4] * J._J[15] + T._T[5] * J._J[16] + T._T[6] * J._J[17];
	gReal h4 = T._T[4] * J._J[16] + T._T[5] * J._J[18] + T._T[6] * J._J[19];
	gReal h5 = T._T[4] * J._J[17] + T._T[5] * J._J[19] + T._T[6] * J._J[20];

	_J[0] += f0 * T._T[0] + f1 * T._T[1] + f2 * T._T[2];
	_J[1] += f0 * T._T[4] + f1 * T._T[5] + f2 * T._T[6];
	_J[2] += f0 * T._T[8] + f1 * T._T[9] + f2 * T._T[10];
	_J[3] += f3 * T._T[0] + f4 * T._T[1] + f5 * T._T[2];
	_J[4] += f3 * T._T[4] + f4 * T._T[5] + f5 * T._T[6];
	_J[5] += f3 * T._T[8] + f4 * T._T[9] + f5 * T._T[10];
	_J[6] += f6 * T._T[4] + f7 * T._T[5] + f8 * T._T[6];
	_J[7] += f6 * T._T[8] + f7 * T._T[9] + f8 * T._T[10];
	_J[8] += g0 * T._T[0] + g1 * T._T[1] + g2 * T._T[2];
	_J[9] += g0 * T._T[4] + g1 * T._T[5] + g2 * T._T[6];
	_J[10] += g0 * T._T[8] + g1 * T._T[9] + g2 * T._T[10];
	_J[11] += (T._T[8] * e0 + T._T[9] * e3 + T._T[10] * e6) * T._T[8] + (T._T[8] * e3 + T._T[9] * e4 + T._T[10] * e7) * T._T[9] + (T._T[8] * e6 + T._T[9] * e7 + T._T[10] * e8) * T._T[10];
	_J[12] += g3 * T._T[0] + g4 * T._T[1] + g5 * T._T[2];
	_J[13] += g3 * T._T[4] + g4 * T._T[5] + g5 * T._T[6];
	_J[14] += g3 * T._T[8] + g4 * T._T[9] + g5 * T._T[10];
	_J[15] += h0 * T._T[0] + h1 * T._T[1] + h2 * T._T[2];
	_J[16] += h0 * T._T[4] + h1 * T._T[5] + h2 * T._T[6];
	_J[17] += h0 * T._T[8] + h1 * T._T[9] + h2 * T._T[10];
	_J[18] += h3 * T._T[4] + h4 * T._T[5] + h5 * T._T[6];
	_J[19] += h3 * T._T[8] + h4 * T._T[9] + h5 * T._T[10];
	_J[20] += (T._T[8] * J._J[15] + T._T[9] * J._J[16] + T._T[10] * J._J[17]) * T._T[8] + (T._T[8] * J._J[16] + T._T[9] * J._J[18] + T._T[10] * J._J[19]) * T._T[9] + (T._T[8] * J._J[17] + T._T[9] * J._J[19] + T._T[10] * J._J[20]) * T._T[10];
}

inline AInertia AInertia::Transform_old(const SE3 &T) const
{
	// operation count: multiplication = 294, addition = 189, subtract = 21
	AInertia re;

	gReal _R[] = {	T._T[0], T._T[1], T._T[2], T._T[4], T._T[5], T._T[6], T._T[8], T._T[9], T._T[10] };

	gReal _D[] = {	_J[3] + T._T[14] * _J[16] - T._T[13] * _J[17], _J[8] - T._T[14] * _J[15] + T._T[12] * _J[17], _J[12] + T._T[13] * _J[15] - T._T[12] * _J[16],
					_J[4] + T._T[14] * _J[18] - T._T[13] * _J[19], _J[9] - T._T[14] * _J[16] + T._T[12] * _J[19], _J[13] + T._T[13] * _J[16] - T._T[12] * _J[18],
					_J[5] + T._T[14] * _J[19] - T._T[13] * _J[20], _J[10] - T._T[14] * _J[17] + T._T[12] * _J[20], _J[14] + T._T[13] * _J[17] - T._T[12] * _J[19] };

	gReal _E[] = {	_J[0] + T._T[14] * _J[4] - T._T[13] * _J[5] + _D[3] * T._T[14] - _D[6] * T._T[13], 0.0, 0.0,
					_J[1] + T._T[14] * _J[9] - T._T[13] * _J[10] - _D[0] * T._T[14] + _D[6] * T._T[12], _J[6] - T._T[14] * _J[8] + T._T[12] * _J[10] - _D[1] * T._T[14] + _D[7] * T._T[12], 0.0,
					_J[2] + T._T[14] * _J[13] - T._T[13] * _J[14] + _D[0] * T._T[13] - _D[3] * T._T[12], _J[7] - T._T[14] * _J[12] + T._T[12] * _J[14] + _D[1] * T._T[13] - _D[4] * T._T[12], _J[11] + T._T[13] * _J[12] - T._T[12] * _J[13] + _D[2] * T._T[13] - _D[5] * T._T[12]	};

	re._J[0] = (_R[0] * _E[0] + _R[1] * _E[3] + _R[2] * _E[6]) * _R[0] + (_R[0] * _E[3] + _R[1] * _E[4] + _R[2] * _E[7]) * _R[1] + (_R[0] * _E[6] + _R[1] * _E[7] + _R[2] * _E[8]) * _R[2];
	re._J[1] = (_R[0] * _E[0] + _R[1] * _E[3] + _R[2] * _E[6]) * _R[3] + (_R[0] * _E[3] + _R[1] * _E[4] + _R[2] * _E[7]) * _R[4] + (_R[0] * _E[6] + _R[1] * _E[7] + _R[2] * _E[8]) * _R[5];
	re._J[6] = (_R[3] * _E[0] + _R[4] * _E[3] + _R[5] * _E[6]) * _R[3] + (_R[3] * _E[3] + _R[4] * _E[4] + _R[5] * _E[7]) * _R[4] + (_R[3] * _E[6] + _R[4] * _E[7] + _R[5] * _E[8]) * _R[5];
	re._J[2] = (_R[0] * _E[0] + _R[1] * _E[3] + _R[2] * _E[6]) * _R[6] + (_R[0] * _E[3] + _R[1] * _E[4] + _R[2] * _E[7]) * _R[7] + (_R[0] * _E[6] + _R[1] * _E[7] + _R[2] * _E[8]) * _R[8];
	re._J[7] = (_R[3] * _E[0] + _R[4] * _E[3] + _R[5] * _E[6]) * _R[6] + (_R[3] * _E[3] + _R[4] * _E[4] + _R[5] * _E[7]) * _R[7] + (_R[3] * _E[6] + _R[4] * _E[7] + _R[5] * _E[8]) * _R[8];
	re._J[11] = (_R[6] * _E[0] + _R[7] * _E[3] + _R[8] * _E[6]) * _R[6] + (_R[6] * _E[3] + _R[7] * _E[4] + _R[8] * _E[7]) * _R[7] + (_R[6] * _E[6] + _R[7] * _E[7] + _R[8] * _E[8]) * _R[8];

	re._J[3] = (_R[0] * _D[0] + _R[1] * _D[1] + _R[2] * _D[2]) * _R[0] + (_R[0] * _D[3] + _R[1] * _D[4] + _R[2] * _D[5]) * _R[1] + (_R[0] * _D[6] + _R[1] * _D[7] + _R[2] * _D[8]) * _R[2];
	re._J[8] = (_R[3] * _D[0] + _R[4] * _D[1] + _R[5] * _D[2]) * _R[0] + (_R[3] * _D[3] + _R[4] * _D[4] + _R[5] * _D[5]) * _R[1] + (_R[3] * _D[6] + _R[4] * _D[7] + _R[5] * _D[8]) * _R[2];
	re._J[12] = (_R[6] * _D[0] + _R[7] * _D[1] + _R[8] * _D[2]) * _R[0] + (_R[6] * _D[3] + _R[7] * _D[4] + _R[8] * _D[5]) * _R[1] + (_R[6] * _D[6] + _R[7] * _D[7] + _R[8] * _D[8]) * _R[2];
	re._J[4] = (_R[0] * _D[0] + _R[1] * _D[1] + _R[2] * _D[2]) * _R[3] + (_R[0] * _D[3] + _R[1] * _D[4] + _R[2] * _D[5]) * _R[4] + (_R[0] * _D[6] + _R[1] * _D[7] + _R[2] * _D[8]) * _R[5];
	re._J[9] = (_R[3] * _D[0] + _R[4] * _D[1] + _R[5] * _D[2]) * _R[3] + (_R[3] * _D[3] + _R[4] * _D[4] + _R[5] * _D[5]) * _R[4] + (_R[3] * _D[6] + _R[4] * _D[7] + _R[5] * _D[8]) * _R[5];
	re._J[13] = (_R[6] * _D[0] + _R[7] * _D[1] + _R[8] * _D[2]) * _R[3] + (_R[6] * _D[3] + _R[7] * _D[4] + _R[8] * _D[5]) * _R[4] + (_R[6] * _D[6] + _R[7] * _D[7] + _R[8] * _D[8]) * _R[5];
	re._J[5] = (_R[0] * _D[0] + _R[1] * _D[1] + _R[2] * _D[2]) * _R[6] + (_R[0] * _D[3] + _R[1] * _D[4] + _R[2] * _D[5]) * _R[7] + (_R[0] * _D[6] + _R[1] * _D[7] + _R[2] * _D[8]) * _R[8];
	re._J[10] = (_R[3] * _D[0] + _R[4] * _D[1] + _R[5] * _D[2]) * _R[6] + (_R[3] * _D[3] + _R[4] * _D[4] + _R[5] * _D[5]) * _R[7] + (_R[3] * _D[6] + _R[4] * _D[7] + _R[5] * _D[8]) * _R[8];
	re._J[14] = (_R[6] * _D[0] + _R[7] * _D[1] + _R[8] * _D[2]) * _R[6] + (_R[6] * _D[3] + _R[7] * _D[4] + _R[8] * _D[5]) * _R[7] + (_R[6] * _D[6] + _R[7] * _D[7] + _R[8] * _D[8]) * _R[8];

	re._J[15] = (_R[0] * _J[15] + _R[1] * _J[16] + _R[2] * _J[17]) * _R[0] + (_R[0] * _J[16] + _R[1] * _J[18] + _R[2] * _J[19]) * _R[1] + (_R[0] * _J[17] + _R[1] * _J[19] + _R[2] * _J[20]) * _R[2];
	re._J[16] = (_R[0] * _J[15] + _R[1] * _J[16] + _R[2] * _J[17]) * _R[3] + (_R[0] * _J[16] + _R[1] * _J[18] + _R[2] * _J[19]) * _R[4] + (_R[0] * _J[17] + _R[1] * _J[19] + _R[2] * _J[20]) * _R[5];
	re._J[18] = (_R[3] * _J[15] + _R[4] * _J[16] + _R[5] * _J[17]) * _R[3] + (_R[3] * _J[16] + _R[4] * _J[18] + _R[5] * _J[19]) * _R[4] + (_R[3] * _J[17] + _R[4] * _J[19] + _R[5] * _J[20]) * _R[5];
	re._J[17] = (_R[0] * _J[15] + _R[1] * _J[16] + _R[2] * _J[17]) * _R[6] + (_R[0] * _J[16] + _R[1] * _J[18] + _R[2] * _J[19]) * _R[7] + (_R[0] * _J[17] + _R[1] * _J[19] + _R[2] * _J[20]) * _R[8];
	re._J[19] = (_R[3] * _J[15] + _R[4] * _J[16] + _R[5] * _J[17]) * _R[6] + (_R[3] * _J[16] + _R[4] * _J[18] + _R[5] * _J[19]) * _R[7] + (_R[3] * _J[17] + _R[4] * _J[19] + _R[5] * _J[20]) * _R[8];
	re._J[20] = (_R[6] * _J[15] + _R[7] * _J[16] + _R[8] * _J[17]) * _R[6] + (_R[6] * _J[16] + _R[7] * _J[18] + _R[8] * _J[19]) * _R[7] + (_R[6] * _J[17] + _R[7] * _J[19] + _R[8] * _J[20]) * _R[8];

	return re;
}

inline void AInertia::SubstractAlphaXYt(gReal alpha, const gReal *x, const gReal *y)
{
	if ( x == NULL || y == NULL ) return; 

	gReal ax0 = alpha*x[0];
	gReal ax1 = alpha*x[1];
	gReal ax2 = alpha*x[2];
	gReal ax3 = alpha*x[3];
	gReal ax4 = alpha*x[4];
	gReal ax5 = alpha*x[5];
	_J[0] -= ax0*y[0];	_J[1] -= ax0*y[1];	_J[2] -= ax0*y[2];	 _J[3] -= ax0*y[3];	_J[4] -= ax0*y[4];	_J[5] -= ax0*y[5];
						_J[6] -= ax1*y[1];	_J[7] -= ax1*y[2];	 _J[8] -= ax1*y[3];	_J[9] -= ax1*y[4];	_J[10] -= ax1*y[5];
											_J[11] -= ax2*y[2]; _J[12] -= ax2*y[3];	_J[13] -= ax2*y[4];	_J[14] -= ax2*y[5];
																_J[15] -= ax3*y[3];	_J[16] -= ax3*y[4];	_J[17] -= ax3*y[5];
																					_J[18] -= ax4*y[4];	_J[19] -= ax4*y[5];
																										_J[20] -= ax5*y[5];
}

inline void AInertia::SubstractAlphaSSt(gReal a, const gReal *s)
{
	if ( s == NULL ) return;

	_J[0] -= a*s[0]*s[0];	_J[1] -= a*s[0]*s[1];	_J[2] -= a*s[0]*s[2];	_J[3] -= a*s[0]*s[3];	_J[4] -= a*s[0]*s[4];	_J[5] -= a*s[0]*s[5];
							_J[6] -= a*s[1]*s[1];	_J[7] -= a*s[1]*s[2];	_J[8] -= a*s[1]*s[3];	_J[9] -= a*s[1]*s[4];	_J[10] -= a*s[1]*s[5];
													_J[11] -= a*s[2]*s[2];	_J[12] -= a*s[2]*s[3];	_J[13] -= a*s[2]*s[4];	_J[14] -= a*s[2]*s[5];
																			_J[15] -= a*s[3]*s[3];	_J[16] -= a*s[3]*s[4];	_J[17] -= a*s[3]*s[5];
																									_J[18] -= a*s[4]*s[4];	_J[19] -= a*s[4]*s[5];
																															_J[20] -= a*s[5]*s[5];
}

inline AInertia AInertia::Transform_ad(const se3 &s) const
{
	AInertia re;

	re._J[0] = (gReal)2. * (-_J[2] * s._w[1] + _J[1] * s._w[2] - _J[5] * s._w[4] + _J[4] * s._w[5]);
	re._J[1] = _J[2] * s._w[0] - _J[7] * s._w[1] - _J[0] * s._w[2] + _J[6] * s._w[2] + _J[5] * s._w[3] - _J[10] * s._w[4] - _J[3] * s._w[5] + _J[9] * s._w[5];
	re._J[2] = -_J[1] * s._w[0] + _J[0] * s._w[1] - _J[11] * s._w[1] + _J[7] * s._w[2] - _J[4] * s._w[3] - _J[14] * s._w[4] + _J[3] * s._w[4] + _J[13] * s._w[5];
	re._J[3] = -_J[12] * s._w[1] - _J[5] * s._w[1] + _J[4] * s._w[2] + _J[8] * s._w[2] - _J[17] * s._w[4] + _J[16] * s._w[5];
	re._J[4] = _J[5] * s._w[0] - _J[13] * s._w[1] - _J[3] * s._w[2] + _J[9] * s._w[2] - _J[19] * s._w[4] + _J[18] * s._w[5];
	re._J[5] = -_J[4] * s._w[0] - _J[14] * s._w[1] + _J[3] * s._w[1] + _J[10] * s._w[2] - _J[20] * s._w[4] + _J[19] * s._w[5];

	re._J[6] = (gReal)2. * (_J[7] * s._w[0] - _J[1] * s._w[2] + _J[10] * s._w[3] - _J[8] * s._w[5]);
	re._J[7] = _J[11] * s._w[0] - _J[6] * s._w[0] + _J[1] * s._w[1] - _J[2] * s._w[2] + _J[14] * s._w[3] - _J[9] * s._w[3] + _J[8] * s._w[4] - _J[12] * s._w[5];
	re._J[8] = _J[12] * s._w[0] - _J[10] * s._w[1] - _J[3] * s._w[2] + _J[9] * s._w[2] + _J[17] * s._w[3] - _J[15] * s._w[5];
	re._J[9] = _J[10] * s._w[0] + _J[13] * s._w[0] - _J[4] * s._w[2] - _J[8] * s._w[2] + _J[19] * s._w[3] - _J[16] * s._w[5];
	re._J[10] = _J[14] * s._w[0] - _J[9] * s._w[0] + _J[8] * s._w[1] - _J[5] * s._w[2] + _J[20] * s._w[3] - _J[17] * s._w[5];

	re._J[11] = (gReal)2. * (-_J[7] * s._w[0] + _J[2] * s._w[1] - _J[13] * s._w[3] + _J[12] * s._w[4]);
	re._J[12] = -_J[8] * s._w[0] - _J[14] * s._w[1] + _J[3] * s._w[1] + _J[13] * s._w[2] - _J[16] * s._w[3] + _J[15] * s._w[4];
	re._J[13] = _J[14] * s._w[0] - _J[9] * s._w[0] + _J[4] * s._w[1] - _J[12] * s._w[2] - _J[18] * s._w[3] + _J[16] * s._w[4];
	re._J[14] = -_J[10] * s._w[0] - _J[13] * s._w[0] + _J[12] * s._w[1] + _J[5] * s._w[1] - _J[19] * s._w[3] + _J[17] * s._w[4];

	re._J[15] = (gReal)-2. * _J[17] * s._w[1] + (gReal)2. * _J[16] * s._w[2];
	re._J[16] = _J[17] * s._w[0] - _J[19] * s._w[1] + (-_J[15] + _J[18]) * s._w[2];
	re._J[17] = -_J[16] * s._w[0] + _J[15] * s._w[1] - _J[20] * s._w[1] + _J[19] * s._w[2];

	re._J[18] = (gReal)2. * _J[19] * s._w[0] - (gReal)2. * _J[16] * s._w[2];
	re._J[19] = -_J[18] * s._w[0] + _J[20] * s._w[0] + _J[16] * s._w[1] - _J[17] * s._w[2];

	re._J[20] = (gReal)-2. * _J[19] * s._w[0] + (gReal)2. * _J[17] * s._w[1];

	return re;
}

inline void AInertia::AddTransform_ad(const AInertia &J, const se3 &s)
{
	_J[0] += (gReal)2. * (-J._J[2] * s._w[1] + J._J[1] * s._w[2] - J._J[5] * s._w[4] + J._J[4] * s._w[5]);
	_J[1] += J._J[2] * s._w[0] - J._J[7] * s._w[1] - J._J[0] * s._w[2] + J._J[6] * s._w[2] + J._J[5] * s._w[3] - J._J[10] * s._w[4] - J._J[3] * s._w[5] + J._J[9] * s._w[5];
	_J[2] += -J._J[1] * s._w[0] + J._J[0] * s._w[1] - J._J[11] * s._w[1] + J._J[7] * s._w[2] - J._J[4] * s._w[3] - J._J[14] * s._w[4] + J._J[3] * s._w[4] + J._J[13] * s._w[5];
	_J[3] += -J._J[12] * s._w[1] - J._J[5] * s._w[1] + J._J[4] * s._w[2] + J._J[8] * s._w[2] - J._J[17] * s._w[4] + J._J[16] * s._w[5];
	_J[4] += J._J[5] * s._w[0] - J._J[13] * s._w[1] - J._J[3] * s._w[2] + J._J[9] * s._w[2] - J._J[19] * s._w[4] + J._J[18] * s._w[5];
	_J[5] += -J._J[4] * s._w[0] - J._J[14] * s._w[1] + J._J[3] * s._w[1] + J._J[10] * s._w[2] - J._J[20] * s._w[4] + J._J[19] * s._w[5];

	_J[6] += (gReal)2. * (J._J[7] * s._w[0] - J._J[1] * s._w[2] + J._J[10] * s._w[3] - J._J[8] * s._w[5]);
	_J[7] += J._J[11] * s._w[0] - J._J[6] * s._w[0] + J._J[1] * s._w[1] - J._J[2] * s._w[2] + J._J[14] * s._w[3] - J._J[9] * s._w[3] + J._J[8] * s._w[4] - J._J[12] * s._w[5];
	_J[8] += J._J[12] * s._w[0] - J._J[10] * s._w[1] - J._J[3] * s._w[2] + J._J[9] * s._w[2] + J._J[17] * s._w[3] - J._J[15] * s._w[5];
	_J[9] += J._J[10] * s._w[0] + J._J[13] * s._w[0] - J._J[4] * s._w[2] - J._J[8] * s._w[2] + J._J[19] * s._w[3] - J._J[16] * s._w[5];
	_J[10] += J._J[14] * s._w[0] - J._J[9] * s._w[0] + J._J[8] * s._w[1] - J._J[5] * s._w[2] + J._J[20] * s._w[3] - J._J[17] * s._w[5];

	_J[11] += (gReal)2. * (-J._J[7] * s._w[0] + J._J[2] * s._w[1] - J._J[13] * s._w[3] + J._J[12] * s._w[4]);
	_J[12] += -J._J[8] * s._w[0] - J._J[14] * s._w[1] + J._J[3] * s._w[1] + J._J[13] * s._w[2] - J._J[16] * s._w[3] + J._J[15] * s._w[4];
	_J[13] += J._J[14] * s._w[0] - J._J[9] * s._w[0] + J._J[4] * s._w[1] - J._J[12] * s._w[2] - J._J[18] * s._w[3] + J._J[16] * s._w[4];
	_J[14] += -J._J[10] * s._w[0] - J._J[13] * s._w[0] + J._J[12] * s._w[1] + J._J[5] * s._w[1] - J._J[19] * s._w[3] + J._J[17] * s._w[4];

	_J[15] += (gReal)-2. * J._J[17] * s._w[1] + (gReal)2. * J._J[16] * s._w[2];
	_J[16] += J._J[17] * s._w[0] - J._J[19] * s._w[1] + (-J._J[15] + J._J[18]) * s._w[2];
	_J[17] += -J._J[16] * s._w[0] + J._J[15] * s._w[1] - J._J[20] * s._w[1] + J._J[19] * s._w[2];

	_J[18] += (gReal)2. * J._J[19] * s._w[0] - (gReal)2. * J._J[16] * s._w[2];
	_J[19] += -J._J[18] * s._w[0] + J._J[20] * s._w[0] + J._J[16] * s._w[1] - J._J[17] * s._w[2];

	_J[20] += (gReal)-2. * J._J[19] * s._w[0] + (gReal)2. * J._J[17] * s._w[1];
}

inline void AInertia::SubstractTransform_ad(const AInertia &J, const se3 &s)
{
	_J[0] -= (gReal)2. * (-J._J[2] * s._w[1] + J._J[1] * s._w[2] - J._J[5] * s._w[4] + J._J[4] * s._w[5]);
	_J[1] -= J._J[2] * s._w[0] - J._J[7] * s._w[1] - J._J[0] * s._w[2] + J._J[6] * s._w[2] + J._J[5] * s._w[3] - J._J[10] * s._w[4] - J._J[3] * s._w[5] + J._J[9] * s._w[5];
	_J[2] -= -J._J[1] * s._w[0] + J._J[0] * s._w[1] - J._J[11] * s._w[1] + J._J[7] * s._w[2] - J._J[4] * s._w[3] - J._J[14] * s._w[4] + J._J[3] * s._w[4] + J._J[13] * s._w[5];
	_J[3] -= -J._J[12] * s._w[1] - J._J[5] * s._w[1] + J._J[4] * s._w[2] + J._J[8] * s._w[2] - J._J[17] * s._w[4] + J._J[16] * s._w[5];
	_J[4] -= J._J[5] * s._w[0] - J._J[13] * s._w[1] - J._J[3] * s._w[2] + J._J[9] * s._w[2] - J._J[19] * s._w[4] + J._J[18] * s._w[5];
	_J[5] -= -J._J[4] * s._w[0] - J._J[14] * s._w[1] + J._J[3] * s._w[1] + J._J[10] * s._w[2] - J._J[20] * s._w[4] + J._J[19] * s._w[5];

	_J[6] -= (gReal)2. * (J._J[7] * s._w[0] - J._J[1] * s._w[2] + J._J[10] * s._w[3] - J._J[8] * s._w[5]);
	_J[7] -= J._J[11] * s._w[0] - J._J[6] * s._w[0] + J._J[1] * s._w[1] - J._J[2] * s._w[2] + J._J[14] * s._w[3] - J._J[9] * s._w[3] + J._J[8] * s._w[4] - J._J[12] * s._w[5];
	_J[8] -= J._J[12] * s._w[0] - J._J[10] * s._w[1] - J._J[3] * s._w[2] + J._J[9] * s._w[2] + J._J[17] * s._w[3] - J._J[15] * s._w[5];
	_J[9] -= J._J[10] * s._w[0] + J._J[13] * s._w[0] - J._J[4] * s._w[2] - J._J[8] * s._w[2] + J._J[19] * s._w[3] - J._J[16] * s._w[5];
	_J[10] -= J._J[14] * s._w[0] - J._J[9] * s._w[0] + J._J[8] * s._w[1] - J._J[5] * s._w[2] + J._J[20] * s._w[3] - J._J[17] * s._w[5];

	_J[11] -= (gReal)2. * (-J._J[7] * s._w[0] + J._J[2] * s._w[1] - J._J[13] * s._w[3] + J._J[12] * s._w[4]);
	_J[12] -= -J._J[8] * s._w[0] - J._J[14] * s._w[1] + J._J[3] * s._w[1] + J._J[13] * s._w[2] - J._J[16] * s._w[3] + J._J[15] * s._w[4];
	_J[13] -= J._J[14] * s._w[0] - J._J[9] * s._w[0] + J._J[4] * s._w[1] - J._J[12] * s._w[2] - J._J[18] * s._w[3] + J._J[16] * s._w[4];
	_J[14] -= -J._J[10] * s._w[0] - J._J[13] * s._w[0] + J._J[12] * s._w[1] + J._J[5] * s._w[1] - J._J[19] * s._w[3] + J._J[17] * s._w[4];

	_J[15] -= (gReal)-2. * J._J[17] * s._w[1] + (gReal)2. * J._J[16] * s._w[2];
	_J[16] -= J._J[17] * s._w[0] - J._J[19] * s._w[1] + (-J._J[15] + J._J[18]) * s._w[2];
	_J[17] -= -J._J[16] * s._w[0] + J._J[15] * s._w[1] - J._J[20] * s._w[1] + J._J[19] * s._w[2];

	_J[18] -= (gReal)2. * J._J[19] * s._w[0] - (gReal)2. * J._J[16] * s._w[2];
	_J[19] -= -J._J[18] * s._w[0] + J._J[20] * s._w[0] + J._J[16] * s._w[1] - J._J[17] * s._w[2];

	_J[20] -= (gReal)-2. * J._J[19] * s._w[0] + (gReal)2. * J._J[17] * s._w[1];
}

// 0.5 * ( x * ~y + y * ~x)
inline AInertia KroneckerProduct(const dse3 &x, const dse3 &y)
{
	gReal y_m0 = (gReal)0.5 * y._m[0];
	gReal y_m1 = (gReal)0.5 * y._m[1];
	gReal y_m2 = (gReal)0.5 * y._m[2];
	gReal y_m3 = (gReal)0.5 * y._m[3];
	gReal y_m4 = (gReal)0.5 * y._m[4];
	gReal y_m5 = (gReal)0.5 * y._m[5];

	return AInertia(x._m[0] * y._m[0], 
					x._m[0] * y_m1 + x._m[1] * y_m0, 
					x._m[0] * y_m2 + x._m[2] * y_m0,
					x._m[0] * y_m3 + x._m[3] * y_m0,
					x._m[0] * y_m4 + x._m[4] * y_m0,
					x._m[0] * y_m5 + x._m[5] * y_m0,
					x._m[1] * y._m[1], 
					x._m[1] * y_m2 + x._m[2] * y_m1,
					x._m[1] * y_m3 + x._m[3] * y_m1,
					x._m[1] * y_m4 + x._m[4] * y_m1,
					x._m[1] * y_m5 + x._m[5] * y_m1,
					x._m[2] * y._m[2],
					x._m[2] * y_m3 + x._m[3] * y_m2,
					x._m[2] * y_m4 + x._m[4] * y_m2,
					x._m[2] * y_m5 + x._m[5] * y_m2,
					x._m[3] * y._m[3],
					x._m[3] * y_m4 + x._m[4] * y_m3,
					x._m[3] * y_m5 + x._m[5] * y_m3,
					x._m[4] * y._m[4],
					x._m[4] * y_m5 + x._m[5] * y_m4,
					x._m[5] * y._m[5]);
}

inline void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s)
{
	re._m[0] = J._J[0] * s._w[0] + J._J[1] * s._w[1] + J._J[2] * s._w[2] + J._J[3] * s._w[3] + J._J[4] * s._w[4] + J._J[5] * s._w[5];
	re._m[1] = J._J[1] * s._w[0] + J._J[6] * s._w[1] + J._J[7] * s._w[2] + J._J[8] * s._w[3] + J._J[9] * s._w[4] + J._J[10] * s._w[5];
	re._m[2] = J._J[2] * s._w[0] + J._J[7] * s._w[1] + J._J[11] * s._w[2] + J._J[12] * s._w[3] + J._J[13] * s._w[4] + J._J[14] * s._w[5];
	re._m[3] = J._J[3] * s._w[0] + J._J[8] * s._w[1] + J._J[12] * s._w[2] + J._J[15] * s._w[3] + J._J[16] * s._w[4] + J._J[17] * s._w[5];
	re._m[4] = J._J[4] * s._w[0] + J._J[9] * s._w[1] + J._J[13] * s._w[2] + J._J[16] * s._w[3] + J._J[18] * s._w[4] + J._J[19] * s._w[5];
	re._m[5] = J._J[5] * s._w[0] + J._J[10] * s._w[1] + J._J[14] * s._w[2] + J._J[17] * s._w[3] + J._J[19] * s._w[4] + J._J[20] * s._w[5];
}

inline void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s)
{
	re[0] = J._J[0] * s._w[0] + J._J[1] * s._w[1] + J._J[2] * s._w[2] + J._J[3] * s._w[3] + J._J[4] * s._w[4] + J._J[5] * s._w[5];
	re[1] = J._J[1] * s._w[0] + J._J[6] * s._w[1] + J._J[7] * s._w[2] + J._J[8] * s._w[3] + J._J[9] * s._w[4] + J._J[10] * s._w[5];
	re[2] = J._J[2] * s._w[0] + J._J[7] * s._w[1] + J._J[11] * s._w[2] + J._J[12] * s._w[3] + J._J[13] * s._w[4] + J._J[14] * s._w[5];
	re[3] = J._J[3] * s._w[0] + J._J[8] * s._w[1] + J._J[12] * s._w[2] + J._J[15] * s._w[3] + J._J[16] * s._w[4] + J._J[17] * s._w[5];
	re[4] = J._J[4] * s._w[0] + J._J[9] * s._w[1] + J._J[13] * s._w[2] + J._J[16] * s._w[3] + J._J[18] * s._w[4] + J._J[19] * s._w[5];
	re[5] = J._J[5] * s._w[0] + J._J[10] * s._w[1] + J._J[14] * s._w[2] + J._J[17] * s._w[3] + J._J[19] * s._w[4] + J._J[20] * s._w[5];
}

inline void set_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s)
{
	if ( s == NULL ) { re.SetZero(); return; }

	re._m[0] = J._J[0] * s[0] + J._J[1] * s[1] + J._J[2] * s[2] + J._J[3] * s[3] + J._J[4] * s[4] + J._J[5] * s[5];
	re._m[1] = J._J[1] * s[0] + J._J[6] * s[1] + J._J[7] * s[2] + J._J[8] * s[3] + J._J[9] * s[4] + J._J[10] * s[5];
	re._m[2] = J._J[2] * s[0] + J._J[7] * s[1] + J._J[11] * s[2] + J._J[12] * s[3] + J._J[13] * s[4] + J._J[14] * s[5];
	re._m[3] = J._J[3] * s[0] + J._J[8] * s[1] + J._J[12] * s[2] + J._J[15] * s[3] + J._J[16] * s[4] + J._J[17] * s[5];
	re._m[4] = J._J[4] * s[0] + J._J[9] * s[1] + J._J[13] * s[2] + J._J[16] * s[3] + J._J[18] * s[4] + J._J[19] * s[5];
	re._m[5] = J._J[5] * s[0] + J._J[10] * s[1] + J._J[14] * s[2] + J._J[17] * s[3] + J._J[19] * s[4] + J._J[20] * s[5];
}

inline void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s)
{
	if ( s == NULL ) { matSet_zero(re, 6); return; }

	re[0] = J._J[0] * s[0] + J._J[1] * s[1] + J._J[2] * s[2] + J._J[3] * s[3] + J._J[4] * s[4] + J._J[5] * s[5];
	re[1] = J._J[1] * s[0] + J._J[6] * s[1] + J._J[7] * s[2] + J._J[8] * s[3] + J._J[9] * s[4] + J._J[10] * s[5];
	re[2] = J._J[2] * s[0] + J._J[7] * s[1] + J._J[11] * s[2] + J._J[12] * s[3] + J._J[13] * s[4] + J._J[14] * s[5];
	re[3] = J._J[3] * s[0] + J._J[8] * s[1] + J._J[12] * s[2] + J._J[15] * s[3] + J._J[16] * s[4] + J._J[17] * s[5];
	re[4] = J._J[4] * s[0] + J._J[9] * s[1] + J._J[13] * s[2] + J._J[16] * s[3] + J._J[18] * s[4] + J._J[19] * s[5];
	re[5] = J._J[5] * s[0] + J._J[10] * s[1] + J._J[14] * s[2] + J._J[17] * s[3] + J._J[19] * s[4] + J._J[20] * s[5];
}

inline void set_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s, int num)
{
	if ( s == NULL ) { matSet_zero(re, 6*num); return; }

	for (int i=0; i<num; i++) {
		set_Mult_AInertia_se3(&re[6*i], J, &s[6*i]);
	}
}

inline void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const se3 &s)
{
	re._m[0] += J._J[0] * s._w[0] + J._J[1] * s._w[1] + J._J[2] * s._w[2] + J._J[3] * s._w[3] + J._J[4] * s._w[4] + J._J[5] * s._w[5];
	re._m[1] += J._J[1] * s._w[0] + J._J[6] * s._w[1] + J._J[7] * s._w[2] + J._J[8] * s._w[3] + J._J[9] * s._w[4] + J._J[10] * s._w[5];
	re._m[2] += J._J[2] * s._w[0] + J._J[7] * s._w[1] + J._J[11] * s._w[2] + J._J[12] * s._w[3] + J._J[13] * s._w[4] + J._J[14] * s._w[5];
	re._m[3] += J._J[3] * s._w[0] + J._J[8] * s._w[1] + J._J[12] * s._w[2] + J._J[15] * s._w[3] + J._J[16] * s._w[4] + J._J[17] * s._w[5];
	re._m[4] += J._J[4] * s._w[0] + J._J[9] * s._w[1] + J._J[13] * s._w[2] + J._J[16] * s._w[3] + J._J[18] * s._w[4] + J._J[19] * s._w[5];
	re._m[5] += J._J[5] * s._w[0] + J._J[10] * s._w[1] + J._J[14] * s._w[2] + J._J[17] * s._w[3] + J._J[19] * s._w[4] + J._J[20] * s._w[5];
}

inline void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const se3 &s)
{
	re[0] += J._J[0] * s._w[0] + J._J[1] * s._w[1] + J._J[2] * s._w[2] + J._J[3] * s._w[3] + J._J[4] * s._w[4] + J._J[5] * s._w[5];
	re[1] += J._J[1] * s._w[0] + J._J[6] * s._w[1] + J._J[7] * s._w[2] + J._J[8] * s._w[3] + J._J[9] * s._w[4] + J._J[10] * s._w[5];
	re[2] += J._J[2] * s._w[0] + J._J[7] * s._w[1] + J._J[11] * s._w[2] + J._J[12] * s._w[3] + J._J[13] * s._w[4] + J._J[14] * s._w[5];
	re[3] += J._J[3] * s._w[0] + J._J[8] * s._w[1] + J._J[12] * s._w[2] + J._J[15] * s._w[3] + J._J[16] * s._w[4] + J._J[17] * s._w[5];
	re[4] += J._J[4] * s._w[0] + J._J[9] * s._w[1] + J._J[13] * s._w[2] + J._J[16] * s._w[3] + J._J[18] * s._w[4] + J._J[19] * s._w[5];
	re[5] += J._J[5] * s._w[0] + J._J[10] * s._w[1] + J._J[14] * s._w[2] + J._J[17] * s._w[3] + J._J[19] * s._w[4] + J._J[20] * s._w[5];
}

inline void add_Mult_AInertia_se3(dse3 &re, const AInertia &J, const gReal *s)
{
	if ( s == NULL ) { return; }

	re._m[0] += J._J[0] * s[0] + J._J[1] * s[1] + J._J[2] * s[2] + J._J[3] * s[3] + J._J[4] * s[4] + J._J[5] * s[5];
	re._m[1] += J._J[1] * s[0] + J._J[6] * s[1] + J._J[7] * s[2] + J._J[8] * s[3] + J._J[9] * s[4] + J._J[10] * s[5];
	re._m[2] += J._J[2] * s[0] + J._J[7] * s[1] + J._J[11] * s[2] + J._J[12] * s[3] + J._J[13] * s[4] + J._J[14] * s[5];
	re._m[3] += J._J[3] * s[0] + J._J[8] * s[1] + J._J[12] * s[2] + J._J[15] * s[3] + J._J[16] * s[4] + J._J[17] * s[5];
	re._m[4] += J._J[4] * s[0] + J._J[9] * s[1] + J._J[13] * s[2] + J._J[16] * s[3] + J._J[18] * s[4] + J._J[19] * s[5];
	re._m[5] += J._J[5] * s[0] + J._J[10] * s[1] + J._J[14] * s[2] + J._J[17] * s[3] + J._J[19] * s[4] + J._J[20] * s[5];
}

inline void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s)
{
	if ( s == NULL ) { return; }

	re[0] += J._J[0] * s[0] + J._J[1] * s[1] + J._J[2] * s[2] + J._J[3] * s[3] + J._J[4] * s[4] + J._J[5] * s[5];
	re[1] += J._J[1] * s[0] + J._J[6] * s[1] + J._J[7] * s[2] + J._J[8] * s[3] + J._J[9] * s[4] + J._J[10] * s[5];
	re[2] += J._J[2] * s[0] + J._J[7] * s[1] + J._J[11] * s[2] + J._J[12] * s[3] + J._J[13] * s[4] + J._J[14] * s[5];
	re[3] += J._J[3] * s[0] + J._J[8] * s[1] + J._J[12] * s[2] + J._J[15] * s[3] + J._J[16] * s[4] + J._J[17] * s[5];
	re[4] += J._J[4] * s[0] + J._J[9] * s[1] + J._J[13] * s[2] + J._J[16] * s[3] + J._J[18] * s[4] + J._J[19] * s[5];
	re[5] += J._J[5] * s[0] + J._J[10] * s[1] + J._J[14] * s[2] + J._J[17] * s[3] + J._J[19] * s[4] + J._J[20] * s[5];
}

inline void add_Mult_AInertia_se3(gReal *re, const AInertia &J, const gReal *s, int num)
{
	if ( s == NULL ) { return; }

	for (int i=0; i<num; i++) {
		add_Mult_AInertia_se3(&re[6*i], J, &s[6*i]);
	}
}

inline void matSet(gReal *re, const gReal *a, int n)
{
	for (int i=0; i<n; i++) { *(re++) = *(a++); }
}

inline void matAdd(gReal *re, const gReal *a, int n)
{
	for (int i=0; i<n; i++) { *(re++) += *(a++); }
}

inline void matSubtract(gReal *re, const gReal *a, int n)
{
	for (int i=0; i<n; i++) { *(re++) -= *(a++); }
}

inline void matMult(gReal *re, gReal alpha, int n)
{
	for (int i=0; i<n; i++) { *(re++) *= alpha; }
}

inline void matSet(gReal *re, const gReal *a, int a_row, int a_col)
{
	int n = a_row * a_col;
	for (int i=0; i<n; i++) { *(re++) = *(a++); }
}

inline void matAdd(gReal *re, const gReal *a, int a_row, int a_col)
{
	int n = a_row * a_col;
	for (int i=0; i<n; i++) { *(re++) += *(a++); }
}

inline void matSubtract(gReal *re, const gReal *a, int a_row, int a_col)
{
	int n = a_row * a_col;
	for (int i=0; i<n; i++) { *(re++) -= *(a++); }
}

inline void matMult(gReal *re, gReal alpha, int a_row, int a_col)
{
	int n = a_row * a_col;
	for (int i=0; i<n; i++) { *(re++) *= alpha; }
}

inline void matSet_zero(gReal *re, int n)
{
	for (int i=0; i<n; i++) { *(re++) = 0; }
}

inline void matSet_eye(gReal *re, int r)
{
	matSet_zero(re, r*r);
	for (int i=0; i<r; i++) { re[i*i] = 1; }
}

inline void matSet_transpose(gReal *re, const gReal *a, int a_row, int a_col)
{
	int i = 0, r = a_row, c;
	const gReal *_mt;
	while ( r-- )
	{
		_mt = a + (i++);
		c = a_col;
		while ( c-- )
		{
			*(re++) = *_mt;
			_mt += a_row;
		}
	}
}

inline void matSet_multAB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int i, bc = b_col, k, ar;
	const gReal *tmpa, *tmpb = b;
	gReal sum;
	while ( bc-- )
	{
		ar = a_row;
		i = 0;
		while ( ar-- )
		{
			tmpa = a + (i++);
			sum = (gReal)0.0;
			k = a_col;
			while ( k-- )
			{
				sum += *tmpa * *tmpb;
				tmpa += a_row;
				tmpb++;
			}
			tmpb -= b_row;
			*(re++) = sum;				
		}
		tmpb += b_row;
	}
}

inline void matSet_multAtB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int ac, bc = b_col, ar;
	const gReal *tmpa, *tmpb = b;
	gReal sum;
	while ( bc-- )
	{
		tmpa = a;
		ac = a_col;
		while ( ac-- )
		{
			sum = (gReal)0.0;
			ar = a_row;
			while ( ar-- ) sum += *(tmpa++) * *(tmpb++);
				
			*(re++) = sum;
			tmpb -= b_row;
		}
		tmpb += b_row;
	}
}

inline void matSet_multABt(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int i, j = 0, br = b_row, ar, ac;
	const gReal *tmpa, *tmpb;
	gReal sum;	
	while ( br-- )
	{
		ar = a_row;
		i = 0;
		while ( ar-- )
		{
			tmpa = a + (i++);
			tmpb = b + j;
			sum = (gReal)0.0;
			ac = a_col;
			while ( ac-- )
			{
				sum += *tmpa * *tmpb;
				tmpa += a_row;
				tmpb += b_row;
			}
			*(re++) = sum;
		}
		j++;
	}
}

inline void matAdd_multAB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int i, bc = b_col, k, ar;
	const gReal *tmpa, *tmpb = b;
	gReal sum;
	while ( bc-- )
	{
		ar = a_row;
		i = 0;
		while ( ar-- )
		{
			tmpa = a + (i++);
			sum = (gReal)0.0;
			k = a_col;
			while ( k-- )
			{
				sum += *tmpa * *tmpb;
				tmpa += a_row;
				tmpb++;
			}
			tmpb -= b_row;
			*(re++) += sum;				
		}
		tmpb += b_row;
	}
}

inline void matAdd_multAtB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int ac, bc = b_col, ar;
	const gReal *tmpa, *tmpb = b;
	gReal sum;
	while ( bc-- )
	{
		tmpa = a;
		ac = a_col;
		while ( ac-- )
		{
			sum = (gReal)0.0;
			ar = a_row;
			while ( ar-- ) sum += *(tmpa++) * *(tmpb++);
				
			*(re++) += sum;
			tmpb -= b_row;
		}
		tmpb += b_row;
	}
}

inline void matAdd_multABt(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int i, j = 0, br = b_row, ar, ac;
	const gReal *tmpa, *tmpb;
	gReal sum;	
	while ( br-- )
	{
		ar = a_row;
		i = 0;
		while ( ar-- )
		{
			tmpa = a + (i++);
			tmpb = b + j;
			sum = (gReal)0.0;
			ac = a_col;
			while ( ac-- )
			{
				sum += *tmpa * *tmpb;
				tmpa += a_row;
				tmpb += b_row;
			}
			*(re++) += sum;
		}
		j++;
	}
}

inline void matSubtract_multAB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int i, bc = b_col, k, ar;
	const gReal *tmpa, *tmpb = b;
	gReal sum;
	while ( bc-- )
	{
		ar = a_row;
		i = 0;
		while ( ar-- )
		{
			tmpa = a + (i++);
			sum = (gReal)0.0;
			k = a_col;
			while ( k-- )
			{
				sum += *tmpa * *tmpb;
				tmpa += a_row;
				tmpb++;
			}
			tmpb -= b_row;
			*(re++) -= sum;				
		}
		tmpb += b_row;
	}
}

inline void matSubtract_multAtB(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int ac, bc = b_col, ar;
	const gReal *tmpa, *tmpb = b;
	gReal sum;
	while ( bc-- )
	{
		tmpa = a;
		ac = a_col;
		while ( ac-- )
		{
			sum = (gReal)0.0;
			ar = a_row;
			while ( ar-- ) sum += *(tmpa++) * *(tmpb++);
				
			*(re++) -= sum;
			tmpb -= b_row;
		}
		tmpb += b_row;
	}
}

inline void matSubtract_multABt(gReal *re, const gReal *a, const gReal *b, int a_row, int a_col, int b_row, int b_col)
{
	int i, j = 0, br = b_row, ar, ac;
	const gReal *tmpa, *tmpb;
	gReal sum;	
	while ( br-- )
	{
		ar = a_row;
		i = 0;
		while ( ar-- )
		{
			tmpa = a + (i++);
			tmpb = b + j;
			sum = (gReal)0.0;
			ac = a_col;
			while ( ac-- )
			{
				sum += *tmpa * *tmpb;
				tmpa += a_row;
				tmpb += b_row;
			}
			*(re++) -= sum;
		}
		j++;
	}
}

inline bool matSet_inv22(gReal *re, const gReal *a, gReal eps)
{
	gReal det = a[0]*a[3]-a[2]*a[1];
	if ( fabs(det) < eps ) return false;
	gReal idet = 1/det;
	re[0] = a[3]*idet;  re[2] = -a[2]*idet;
	re[1] = -a[1]*idet; re[3] = a[0]*idet;
	return true;
}

inline bool matSet_inv22sym(gReal *re, const gReal *a, gReal eps)
{
	gReal det = a[0]*a[3]-a[2]*a[1];
	if ( fabs(det) < eps ) return false;
	gReal idet = 1/det;
	re[0] = a[3]*idet;  re[2] = -a[2]*idet;
						re[3] = a[0]*idet;
	re[1] = re[2];
	return true;
}

inline bool matSet_inv33(gReal *re, const gReal *a, gReal eps)
{
	gReal det = a[0]*a[4]*a[8]+a[2]*a[3]*a[7]+a[1]*a[5]*a[6]-(a[2]*a[4]*a[6]+a[1]*a[3]*a[8]+a[0]*a[5]*a[7]);
	if ( fabs(det) < eps ) return false;
	gReal idet = 1/det;
	re[0] = (a[4]*a[8]-a[5]*a[7])*idet;	re[3] = (a[6]*a[5]-a[3]*a[8])*idet;	re[6] = (a[3]*a[7]-a[6]*a[4])*idet;
	re[1] = (a[7]*a[2]-a[1]*a[8])*idet;	re[4] = (a[0]*a[8]-a[6]*a[2])*idet;	re[7] = (a[6]*a[1]-a[0]*a[7])*idet;
	re[2] = (a[1]*a[5]-a[4]*a[2])*idet;	re[5] = (a[3]*a[2]-a[0]*a[5])*idet;	re[8] = (a[0]*a[4]-a[3]*a[1])*idet;
	return true;
}

inline bool matSet_inv33sym(gReal *re, const gReal *a, gReal eps)
{
	gReal det = a[0]*a[4]*a[8]+a[2]*a[3]*a[7]+a[1]*a[5]*a[6]-(a[2]*a[4]*a[6]+a[1]*a[3]*a[8]+a[0]*a[5]*a[7]);
	if ( fabs(det) < eps ) return false;
	gReal idet = 1/det;
	re[0] = (a[4]*a[8]-a[5]*a[7])*idet;	re[3] = (a[6]*a[5]-a[3]*a[8])*idet;	re[6] = (a[3]*a[7]-a[6]*a[4])*idet;
										re[4] = (a[0]*a[8]-a[6]*a[2])*idet;	re[7] = (a[6]*a[1]-a[0]*a[7])*idet;
																			re[8] = (a[0]*a[4]-a[3]*a[1])*idet;
	re[1] = re[3]; re[2] = re[6]; re[5] = re[7];
	return true;
}

inline bool matSet_inv44(gReal *re, const gReal *a, gReal eps)
{
	// copied from a c-code written by Chris Atkeson, 12/11/2012
	gReal det = a[12]*a[9]*a[6]*a[3] - a[8]*a[13]*a[6]*a[3] - a[12]*a[5]*a[10]*a[3] + a[4]*a[13]*a[10]*a[3] + a[8]*a[5]*a[14]*a[3] - a[4]*a[9]*a[14]*a[3] 
	   - a[12]*a[9]*a[2]*a[7] + a[8]*a[13]*a[2]*a[7] + a[12]*a[1]*a[10]*a[7] - a[0]*a[13]*a[10]*a[7] - a[8]*a[1]*a[14]*a[7] + a[0]*a[9]*a[14]*a[7] 
	   + a[12]*a[5]*a[2]*a[11] - a[4]*a[13]*a[2]*a[11] - a[12]*a[1]*a[6]*a[11] + a[0]*a[13]*a[6]*a[11] + a[4]*a[1]*a[14]*a[11] - a[0]*a[5]*a[14]*a[11] 
	   - a[8]*a[5]*a[2]*a[15] + a[4]*a[9]*a[2]*a[15] + a[8]*a[1]*a[6]*a[15] - a[0]*a[9]*a[6]*a[15] - a[4]*a[1]*a[10]*a[15] + a[0]*a[5]*a[10]*a[15];
	if ( fabs(det) < eps ) return false;
	gReal idet = 1/det;
	re[0] = (-(a[13]*a[10]*a[7]) + a[9]*a[14]*a[7] + a[13]*a[6]*a[11] - a[5]*a[14]*a[11] - a[9]*a[6]*a[15] + a[5]*a[10]*a[15])*idet;
	re[4] = (a[12]*a[10]*a[7] - a[8]*a[14]*a[7] - a[12]*a[6]*a[11] + a[4]*a[14]*a[11] + a[8]*a[6]*a[15] - a[4]*a[10]*a[15])*idet;
	re[8] = (-(a[12]*a[9]*a[7]) + a[8]*a[13]*a[7] + a[12]*a[5]*a[11] - a[4]*a[13]*a[11] - a[8]*a[5]*a[15] + a[4]*a[9]*a[15])*idet;
	re[12] = (a[12]*a[9]*a[6] - a[8]*a[13]*a[6] - a[12]*a[5]*a[10] + a[4]*a[13]*a[10] + a[8]*a[5]*a[14] - a[4]*a[9]*a[14])*idet;
	re[1] = (a[13]*a[10]*a[3] - a[9]*a[14]*a[3] - a[13]*a[2]*a[11] + a[1]*a[14]*a[11] + a[9]*a[2]*a[15] - a[1]*a[10]*a[15])*idet;
	re[5] = (-(a[12]*a[10]*a[3]) + a[8]*a[14]*a[3] + a[12]*a[2]*a[11] - a[0]*a[14]*a[11] - a[8]*a[2]*a[15] + a[0]*a[10]*a[15])*idet;
	re[9] = (a[12]*a[9]*a[3] - a[8]*a[13]*a[3] - a[12]*a[1]*a[11] + a[0]*a[13]*a[11] + a[8]*a[1]*a[15] - a[0]*a[9]*a[15])*idet;
	re[13] = (-(a[12]*a[9]*a[2]) + a[8]*a[13]*a[2] + a[12]*a[1]*a[10] - a[0]*a[13]*a[10] - a[8]*a[1]*a[14] + a[0]*a[9]*a[14])*idet;
	re[2] = (-(a[13]*a[6]*a[3]) + a[5]*a[14]*a[3] + a[13]*a[2]*a[7] - a[1]*a[14]*a[7] - a[5]*a[2]*a[15] + a[1]*a[6]*a[15])*idet;
	re[6] = (a[12]*a[6]*a[3] - a[4]*a[14]*a[3] - a[12]*a[2]*a[7] + a[0]*a[14]*a[7] + a[4]*a[2]*a[15] - a[0]*a[6]*a[15])*idet;
	re[10] = (-(a[12]*a[5]*a[3]) + a[4]*a[13]*a[3] + a[12]*a[1]*a[7] - a[0]*a[13]*a[7] - a[4]*a[1]*a[15] + a[0]*a[5]*a[15])*idet;
	re[14] = (a[12]*a[5]*a[2] - a[4]*a[13]*a[2] - a[12]*a[1]*a[6] + a[0]*a[13]*a[6] + a[4]*a[1]*a[14] - a[0]*a[5]*a[14])*idet;
	re[3] = (a[9]*a[6]*a[3] - a[5]*a[10]*a[3] - a[9]*a[2]*a[7] + a[1]*a[10]*a[7] + a[5]*a[2]*a[11] - a[1]*a[6]*a[11])*idet;
	re[7] = (-(a[8]*a[6]*a[3]) + a[4]*a[10]*a[3] + a[8]*a[2]*a[7] - a[0]*a[10]*a[7] - a[4]*a[2]*a[11] + a[0]*a[6]*a[11])*idet;
	re[11] = (a[8]*a[5]*a[3] - a[4]*a[9]*a[3] - a[8]*a[1]*a[7] + a[0]*a[9]*a[7] + a[4]*a[1]*a[11] - a[0]*a[5]*a[11])*idet;
	re[15] = (-(a[8]*a[5]*a[2]) + a[4]*a[9]*a[2] + a[8]*a[1]*a[6] - a[0]*a[9]*a[6] - a[4]*a[1]*a[10] + a[0]*a[5]*a[10])*idet;
	return true;
}

inline bool matSet_inv44sym(gReal *re, const gReal *a, gReal eps)
{
	gReal det = a[12]*a[9]*a[6]*a[3] - a[8]*a[13]*a[6]*a[3] - a[12]*a[5]*a[10]*a[3] + a[4]*a[13]*a[10]*a[3] + a[8]*a[5]*a[14]*a[3] - a[4]*a[9]*a[14]*a[3] 
	   - a[12]*a[9]*a[2]*a[7] + a[8]*a[13]*a[2]*a[7] + a[12]*a[1]*a[10]*a[7] - a[0]*a[13]*a[10]*a[7] - a[8]*a[1]*a[14]*a[7] + a[0]*a[9]*a[14]*a[7] 
	   + a[12]*a[5]*a[2]*a[11] - a[4]*a[13]*a[2]*a[11] - a[12]*a[1]*a[6]*a[11] + a[0]*a[13]*a[6]*a[11] + a[4]*a[1]*a[14]*a[11] - a[0]*a[5]*a[14]*a[11] 
	   - a[8]*a[5]*a[2]*a[15] + a[4]*a[9]*a[2]*a[15] + a[8]*a[1]*a[6]*a[15] - a[0]*a[9]*a[6]*a[15] - a[4]*a[1]*a[10]*a[15] + a[0]*a[5]*a[10]*a[15];
	if ( fabs(det) < eps ) return false;
	gReal idet = 1/det;
	re[0] = (-(a[13]*a[10]*a[7]) + a[9]*a[14]*a[7] + a[13]*a[6]*a[11] - a[5]*a[14]*a[11] - a[9]*a[6]*a[15] + a[5]*a[10]*a[15])*idet;
	re[4] = (a[12]*a[10]*a[7] - a[8]*a[14]*a[7] - a[12]*a[6]*a[11] + a[4]*a[14]*a[11] + a[8]*a[6]*a[15] - a[4]*a[10]*a[15])*idet;
	re[8] = (-(a[12]*a[9]*a[7]) + a[8]*a[13]*a[7] + a[12]*a[5]*a[11] - a[4]*a[13]*a[11] - a[8]*a[5]*a[15] + a[4]*a[9]*a[15])*idet;
	re[12] = (a[12]*a[9]*a[6] - a[8]*a[13]*a[6] - a[12]*a[5]*a[10] + a[4]*a[13]*a[10] + a[8]*a[5]*a[14] - a[4]*a[9]*a[14])*idet;
	re[1] = re[4];
	re[5] = (-(a[12]*a[10]*a[3]) + a[8]*a[14]*a[3] + a[12]*a[2]*a[11] - a[0]*a[14]*a[11] - a[8]*a[2]*a[15] + a[0]*a[10]*a[15])*idet;
	re[9] = (a[12]*a[9]*a[3] - a[8]*a[13]*a[3] - a[12]*a[1]*a[11] + a[0]*a[13]*a[11] + a[8]*a[1]*a[15] - a[0]*a[9]*a[15])*idet;
	re[13] = (-(a[12]*a[9]*a[2]) + a[8]*a[13]*a[2] + a[12]*a[1]*a[10] - a[0]*a[13]*a[10] - a[8]*a[1]*a[14] + a[0]*a[9]*a[14])*idet;
	re[2] = re[8];
	re[6] = re[9];
	re[10] = (-(a[12]*a[5]*a[3]) + a[4]*a[13]*a[3] + a[12]*a[1]*a[7] - a[0]*a[13]*a[7] - a[4]*a[1]*a[15] + a[0]*a[5]*a[15])*idet;
	re[14] = (a[12]*a[5]*a[2] - a[4]*a[13]*a[2] - a[12]*a[1]*a[6] + a[0]*a[13]*a[6] + a[4]*a[1]*a[14] - a[0]*a[5]*a[14])*idet;
	re[3] = re[12];
	re[7] = re[13];
	re[11] = re[14];
	re[15] = (-(a[8]*a[5]*a[2]) + a[4]*a[9]*a[2] + a[8]*a[1]*a[6] - a[0]*a[9]*a[6] - a[4]*a[1]*a[10] + a[0]*a[5]*a[10])*idet;
	return true;
}

inline bool matSet_invNN(gReal *re, const gReal *a, int n)
{
	int i, info, nn = n*n;
	int *ipvt = new int[n];
	gReal *a2 = new gReal[nn];
	matSet(a2, a, nn); // a2=a
	for (i=0; i<nn; i++) { re[i] = 0.0; }
	__dgefa(a2, n, n, ipvt, info);
	if ( info != 0 ) { delete ipvt; delete a2; return false; }
	for (i=0; i<n; i++ ) {
		re[i+n*i] = 1.0;
		__dgesl(a2, n, n, ipvt, &re[n*i], 0);
	}
	delete ipvt; delete a2;
	return true;
}

inline bool matSet_invNNfast(gReal *re, gReal *a, int n)
{
	if ( n > _MAX_SIZE_INVNNFAST ) return false;
	int i, info, nn = n*n;
	int ipvt[_MAX_SIZE_INVNNFAST];
	for (i=0; i<nn; i++) { re[i] = 0.0; }
	__dgefa(a, n, n, ipvt, info);
	if ( info != 0 ) return false;
	for (i=0; i<n; i++ ) {
		re[i+n*i] = 1.0;
		__dgesl(a, n, n, ipvt, &re[n*i], 0);
	}
	return true;
}

inline void matSet_ad(gReal *re, const gReal *s1, const gReal *s2)
{
	if ( s1 == NULL || s2 == NULL ) { matSet_zero(re, 6); return; }

	re[0] = s1[1] * s2[2] - s1[2] * s2[1];
	re[1] = s1[2] * s2[0] - s1[0] * s2[2];
	re[2] = s1[0] * s2[1] - s1[1] * s2[0];
	re[3] = s1[1] * s2[5] - s1[2] * s2[4] - s2[1] * s1[5] + s2[2] * s1[4];
	re[4] = s1[2] * s2[3] - s1[0] * s2[5] - s2[2] * s1[3] + s2[0] * s1[5];
	re[5] = s1[0] * s2[4] - s1[1] * s2[3] - s2[0] * s1[4] + s2[1] * s1[3];
}

inline void matAdd_ad(gReal *re, const gReal *s1, const gReal *s2)
{
	if ( s1 == NULL || s2 == NULL ) { return; }

	re[0] += s1[1] * s2[2] - s1[2] * s2[1];
	re[1] += s1[2] * s2[0] - s1[0] * s2[2];
	re[2] += s1[0] * s2[1] - s1[1] * s2[0];
	re[3] += s1[1] * s2[5] - s1[2] * s2[4] - s2[1] * s1[5] + s2[2] * s1[4];
	re[4] += s1[2] * s2[3] - s1[0] * s2[5] - s2[2] * s1[3] + s2[0] * s1[5];
	re[5] += s1[0] * s2[4] - s1[1] * s2[3] - s2[0] * s1[4] + s2[1] * s1[3];
}

