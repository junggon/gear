//====================================================================================================
//
//      Title :   rmatrix3j.cpp
//
//      Note  :   This is a modified version of the original work by Jinwook Kim (see below).
//                Modified and maintained by Junggon Kim (junggon@gmail.com)
//
//====================================================================================================

//////////////////////////////////////////////////////////////////////////////////
//
//		title		:	rmatrix3.cpp
//						
//		version		:	v2.89
//		author		:	Jinwook Kim (zinook@plaza.snu.ac.kr)
//		last update	:	2001.7.23
//
//		Note		:
//
//////////////////////////////////////////////////////////////////////////////////



#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "rmatrix3j.h"

#ifdef _MSC_VER
#include <windows.h>
static LARGE_INTEGER _tstart, _tend, _tfreq;
#else
#include <time.h>
static clock_t _start;
#endif

void tic() 
{ 
#ifdef _MSC_VER
	::QueryPerformanceCounter(&_tstart);
#else
	_start = clock(); 
#endif
}

gReal toc() 
{ 
#ifdef _MSC_VER
	::QueryPerformanceCounter(&_tend);
	::QueryPerformanceFrequency(&_tfreq);
	return gReal(_tend.QuadPart-_tstart.QuadPart)/gReal(_tfreq.QuadPart);
#else
	return (gReal)( clock() - _start ) / (gReal) CLOCKS_PER_SEC;
#endif
}


RMatrix Zeros(int r, int c)
{
	RMatrix re(r, c);
	int n = r * c;
	gReal *_r = re.element;
	while ( n-- ) *(_r++) = 0.0;
	return re;
}

RMatrix Ones(int r, int c)
{
	RMatrix re(r, c);
	int n = r * c;
	gReal *_r = re.element;
	while ( n-- ) *(_r++) = 1.0;
	return re;
}

RMatrix Rand(int r, int c)
{
	//srand( (unsigned)time( NULL ) );
	RMatrix re(r, c);
	int n = r * c;
	gReal *_r = re.element;
	while ( n-- ) *(_r++) = drand(1.0);
	return re;		
}

RMatrix Eye(int r, int c )
{
	RMatrix re(r,c);
	int n = r * c;
	gReal *_r = re.element;
	while ( n-- ) *(_r++) = 0.0;
	
	c = _MIN(r++,c);
	_r = re.element;
	while ( c-- )
	{
		*_r = 1.0;
		_r += r;
	}		
	return re;
}

RMatrix Eye(int r)
{
	return Eye(r,r);
}

RMatrix Inv(const RMatrix &m)
{
	return m % Eye(m.row, m.row);
}

int idamax(int n, gReal *dx)
{
	gReal dmax;
	int i, idamax = 0;
      
	if ( n < 1 ) return 0;
	if ( n == 1 ) return 1;
	
	dmax = fabs(dx[0]);
	for ( i = 1; i < n; i++ )
	{
		if ( fabs(dx[i]) > dmax )
		{
			idamax = i;
			dmax = fabs(dx[i]);
		}
	}
	return idamax;
}

void _dgefa(gReal *x, int lda, int n, int *jpvt, int &info)
{
	gReal t, *xk = x, *xj;
	int i, j, k, l;
	// gaussian elimination with partial pivoting
	info = 0;
	if ( n > 1 )
	{
		for ( k = 0; k < n - 1; k++, xk += lda )
		{
			// find l = pivot index
			l = idamax(n-k, xk+k) + k;
			jpvt[k] = l;
			// zero pivot implies this column already triangularized
			if ( xk[l] == 0.0 ) info = k;
			else
			{
				// interchange if necessary
				if ( l != k )
				{
					t = xk[l];
					xk[l] = xk[k];
					xk[k] = t;
				}
				// compute multipliers
				t = (gReal)-1.0 / xk[k];
				//for ( j = 1; j < n-k; j++ ) xk[k+j] *= t;
				for ( j = 1+k; j < n; j++ ) xk[j] *= t;
				// row elimination with column indexing
				for ( j = k+1, xj = xk+lda; j < n; j++, xj += lda )
				{
					t = xj[l];
					if ( l != k )
					{
						xj[l] = xj[k];
						xj[k] = t;
					}
					for ( i = 1+k; i < n; i++ ) xj[i] += t * xk[i];
				}
	        }
		}
	} else k = 0;

	jpvt[k] = k;
	if ( xk[k] == 0.0 ) info = k;
	return;
}

void _dgesl(gReal *x, int lda, int n, int *jpvt, gReal *b, int job)
{
	gReal t, *xk = x;
	int k, l;

	if ( job == 0 ) 
	{
		// job = 0 , solve  a * x = b
		// first solve  l*y = b
		if ( n >= 2 )
		{
			for( k = 0; k < n-1; k++ )
			{
				l = jpvt[k];
				t = b[l];
				if ( l != k )
				{
					b[l] = b[k];
					b[k] = t;
				}
				for ( l = k+1; l < n; l++ ) b[l] += t * xk[l];
				xk += lda;				
			}
		}
		// now solve  u*x = y
		for ( k = n-1; k >= 0; k-- )
		{
			b[k] /= xk[k];
			t = -b[k];
			for ( l = 0; l < k; l++ ) b[l] += t * xk[l];
			xk -= lda;			
		}
		return;
	}

	// job = nonzero, solve  trans(a) * x = b
	// first solve  trans(u)*y = b
	for ( k = 0; k < n; k++ )
	{
		t = 0.0;
		for ( l = 0; l < k; l++ ) t += xk[l] * b[l];
		b[k] = (b[k] - t) / xk[k];
		xk += lda;
	}
	// now solve trans(l)*x = y
  	if ( n >= 2 )
	{
		xk--;
		for ( k = n-1; k >= 0; k-- )
		{
			t = 0.0;			
			for ( l = 1; l < n-k; l++ ) t += xk[l] * b[k+l];
			b[k] += t;

			l = jpvt[k];
			if ( l != k )
			{
				t = b[l];
				b[l] = b[k];
				b[k] = t;
			}
			xk -= lda + 1;
		}
	}
	return;
}

// pivot x and b
// r, c : size of x
// k : starting pivot
// ipvt, jpvt : pivot index
void _fullpivoting(gReal *x, int r, int c, int k, int *ipvt, int *jpvt)
{
	int i, j, imax, jmax, itmp;
	gReal _max, dtmp, *xk = x + k + k * r, *xj;
	
	_max = fabs(*xk);
	imax = jmax = k;

	for ( j = k; j < c; j++ )
	{
		for ( i = k; i < r; i++ )
		{
			dtmp = fabs(*(xk++));
			if ( dtmp > _max ) { _max = dtmp; imax = i; jmax = j; }
		}
		xk += k;
	}
		
	// row swapping
	if ( imax != k )
	{
		itmp = ipvt[imax];	ipvt[imax] = ipvt[k]; 	ipvt[k] = itmp;
		for ( i = k, xk = x + k * r; i < c; i++, xk += r ) 
		{ 
			dtmp = xk[imax];	
			xk[imax] = xk[k];
			xk[k] = dtmp;			
		}
	}
	
	// column swapping
	if ( jmax != k ) 
	{
		itmp = jpvt[jmax];	jpvt[jmax] = jpvt[k];	jpvt[k] = itmp;
		for ( i = 0, xk = x + k * r, xj = x + jmax * r; i < r; i++, xk++, xj++ ) 
		{	
			dtmp = *xj;
			*xj = *xk;
			*xk = dtmp;			
		}
	}
}

// x a = b
// elimination process on x and b
// full pivoting on x and row pivoting on b
// if u want to ignore b, set bc = 0
// return value : rank of x
//
// x [ r X c ]
// b [ r X bc ]
// ipvt [ r ]
// jpvt [ c ]
// zero_tolerance : criterion for determining zero, 1e-8 will be good for usual case

int _gauss_elimination(gReal *x, int r, int c, int *ipvt, int *jpvt, gReal zero_tolerance)
{
	int i, j, k;
	gReal t, *xi = x, *xk;

	for ( i = 0; i < r; i++ ) ipvt[i] = i;
	for ( j = 0; j < c; j++ ) jpvt[j] = j;

	for ( i = 0; i < _MIN(r,c); i++, xi += r )
	{
		_fullpivoting(x, r, c, i, ipvt, jpvt);
		if ( fabs(xi[i]) < zero_tolerance ) return i;
		else
		{
			for ( j = i+1; j < r; j++ )
			{
				t =  - xi[j] / xi[i];
				xi[j] = 0.0;
				for ( k = i+1, xk = xi + r; k < c; k++, xk += r ) xk[j] += t * xk[i];
			}
		}
	}
	return _MIN(r,c);
}

int _dpofa(gReal *a, int n)
{
	gReal s, t;
	int i, j, k, info;
	// begin block with ...exits to 40
	for ( j = 0; j < n; j++ )
	{
		info = j;
		s = 0.0;
		for ( k = 0; k < j; k++ )
		{
			for ( t = 0.0, i = 0; i < k; i++ ) t += a[i+k*n] * a[i+j*n];
			t = a[k+j*n] - t;
			t /= a[k+k*n];
			a[k+j*n] = t;
			s += t * t;
		}
		s = a[j+j*n] - s;
		// exit
		if ( s <= 0.0 ) return info;
		a[j+j*n] = sqrt(s);
	}	
	return 0;
}

void _dposl(gReal *a, int n, gReal *b)
{
	gReal t;
	int i, k;

	// solve trans(r)*y = b
	for ( k = 0; k < n; k++ )
	{
		for ( t = 0.0, i = 0; i < k; i++ ) t += a[i+k*n] * b[i];
		b[k] = (b[k] - t) / a[k+k*n];
	}

	// solve r*x = y
	for ( k = n-1; k >= 0; k-- )
	{
		b[k] /= a[k+k*n];
		t = -b[k];
		for ( i = 0; i < k; i++ ) b[i] += t * a[i+k*n];		
	}
	return;
}

void _conv(gReal re[], const gReal u[], int m, const gReal v[], int n)
{
	int i, j;
	gReal sum;

	for ( i = 0; i < m+n-1; i++ )
	{
		sum = 0.0;
		for ( j = _MAX(0,1+i-n); j <= _MIN(i,m-1); j++ ) sum += u[j] * v[i-j];
		re[i] = sum;
	}
}

RMatrix Conv(const RMatrix &u, const RMatrix &v)
{
	int m = u.row*u.col, n = v.row*v.col;
	RMatrix re(m+n-1,1);
	_conv(re.element, u.element, m, v.element, n);
	return re;
}

bool SolveAxEqualB(const RMatrix &A, RMatrix &x, const RMatrix &B)
{
	if ( A.row * A.col == 1 ) 
	{
		if ( A[0] == 0.0 ) return false;
		x = B / A[0];
		return true;
	}
	if ( A.row != B.row ) return false;
	if ( A.row != A.col )
	{
		if ( A.row > A.col ) return SolveAxEqualB(A ^ A, x, A ^ B);
		bool flag = SolveAxEqualB(A | A, x, B);
		x = A ^ x;
		return flag;
		//return QRSolveAxEqualB(A, x, B);
	}
	int info, i;
	static IMatrix _ipvt_SlvAxB;
	static RMatrix _A_SlvAxB;
	
	_A_SlvAxB = A;
	x = B;

	if ( _ipvt_SlvAxB.row < _A_SlvAxB.row ) _ipvt_SlvAxB.ReNew(_A_SlvAxB.row, 1);
	
	_dgefa(_A_SlvAxB.element, _A_SlvAxB.row, _A_SlvAxB.col, _ipvt_SlvAxB.element, info);
	if ( info != 0 ) return false;
	for ( i = 0; i < x.col; i++ ) _dgesl(_A_SlvAxB.element, _A_SlvAxB.row, _A_SlvAxB.col, _ipvt_SlvAxB.element, x.element+x.row*i, 0);
	
	return true;
}

//_rmatrix operator & (const _rmatrix &m) const
bool SolveAtxEqualB(const RMatrix &A, RMatrix &x, const RMatrix &B)
{
	if ( A.row * A.col == 1 )
	{
		if ( A[0] == 0.0 ) return false;
		x = B / A[0];
		return true;
	}
	if ( A.col != B.row ) return false;
	if ( A.row != A.col )
		return SolveAxEqualB(~A, x, B);
		// return QRSolveAtxEqualB(A, x, B);
	
	int info, i;
	static IMatrix _ipvt_SlvAxB;
	static RMatrix _A_SlvAxB;
	
	_A_SlvAxB = A;
	x = B;

	if ( _ipvt_SlvAxB.row < _A_SlvAxB.row ) _ipvt_SlvAxB.ReNew(_A_SlvAxB.row, 1);
	
	_dgefa(_A_SlvAxB.element, _A_SlvAxB.row, _A_SlvAxB.col, _ipvt_SlvAxB.element, info);
	if ( info != 0 ) return false;
	for ( i = 0; i < x.col; i++ ) _dgesl(_A_SlvAxB.element, _A_SlvAxB.row, _A_SlvAxB.col, _ipvt_SlvAxB.element, x.element+x.row*i, 1);
	
	return true;
}

gReal Det(RMatrix A)
{
	int info, i;
	static IMatrix _ipvt_Det;
	
	if ( _ipvt_Det.row < A.row ) _ipvt_Det.ReNew(A.row, 1);
	gReal re = 1.0;

	_dgefa(A.element, A.row, A.col, _ipvt_Det.element, info);
	
	for ( i = 0; i < A.row; i++ )
	{
		if ( i != _ipvt_Det.element[i] ) re = -re;
		re *= A.element[i+i*A.row];
	}
	return re;
}

bool SolvePosDefAxEqualB(RMatrix &A, RMatrix &x, const RMatrix &B)
{
	if ( A.row != A.col ) return false;
	x = B;
	if ( _dpofa(A.element, A.row) != 0 ) return false;
	for ( int i = 0; i < B.col; i++ )
		_dposl(A.element, A.row, x.element+i*x.row);
	return true;
}

RMatrix Companion(const RMatrix &m)
{
	if ( m.row > 1 && m.col > 1 ) std::cerr << "Companion : not vector";
	int n = _MAX(m.row, m.col) - 1;
	RMatrix re(n, n);
	gReal div = (gReal)1.0 / m.element[0];
	for ( int i = 0 ; i < n; i++ )
	for ( int j = 0 ; j < n; j++ )
	{
		if ( i == 0 ) re.element[j*n] = - (gReal)m.element[j+1] * div;
		else if ( i == j+1 ) re.element[i+j*n] = 1.0;
		else re.element[i+j*n] = 0.0;
	}
	return re;
}

RMatrix Roots(const RMatrix& cof)
{
	int r = cof.row, c = cof.col, n = _MAX(r,c);
	if ( r > 1 && c > 1 ) std::cerr << "Roots : not vector" << std::endl;
	int i = 0;
	while ( 1 )
		if ( cof.element[i] != 0.0 || ++i >= n ) break;
	
	if ( i == 0 ) return Eig(Companion(cof));
	if ( i == n ) return RMatrix();
	RMatrix cof2(n-i,1);
	for ( int j = 0; j < n-i; j++ ) cof2.element[j] = cof.element[j+i];
	return Eig(Companion(cof2));
}

void _change(gReal *v, const int i, const int j, int *idx = NULL)
{
	gReal tmp;
	int itmp;
	tmp = v[i];
	v[i] = v[j];
	v[j] = tmp;
	if ( idx != NULL ) 
	{
		itmp = idx[i];
		idx[i] = idx[j];
		idx[j] = itmp;
	}
}

void _qsort(gReal *v, const int n, const int left, const int right, int *idx = NULL)
{
	int i, last; 
	if ( left >= right ) return;
	_change(v, left, (left+right)/2, idx );
	last = left;
	for ( i = left + 1; i <= right; i++ )
	{
		if ( v[i] < v[left] ) _change(v, ++last, i, idx);
	}
	_change(v, left, last, idx);
	_qsort(v, n, left, last-1, idx);
	_qsort(v, n, last+1, right, idx);
}

// need srand((unsigned)time(NULL));
gReal drand(gReal range)
{
	return (gReal)2.0 * range * ( (gReal)rand() / (gReal)RAND_MAX - (gReal)0.5 );
}

gReal prand(gReal range)
{
	return range * ( (gReal)rand() / (gReal)RAND_MAX );
}

int GaussElimination(RMatrix &A, IMatrix &row_pivot, IMatrix &column_pivot)
{
	return GaussElimination(A, row_pivot, column_pivot, (gReal)1e-6);
}

int GaussElimination(RMatrix &A, IMatrix &row_pivot, IMatrix &column_pivot, gReal eps)
{
	row_pivot.ReNew(A.row); column_pivot.ReNew(A.col);
	return _gauss_elimination(A.element, A.row, A.col, row_pivot.element, column_pivot.element, eps);
}

//----------------------------------------------------------------------
//
//		title		:	EISPACK C Version
//						
//		version		:	v1.5
//		author		:	Jinwook Kim (zinook@plaza.snu.ac.kr)
//		last update	:	2001.7.23
//
//----------------------------------------------------------------------
//     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
//     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
//     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
//     OF A gReal GENERAL MATRIX.
//
//     ON INPUT
//
//        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
//        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
//        DIMENSION STATEMENT.
//
//        N  IS THE ORDER OF THE MATRIX  A.
//
//        A  CONTAINS THE gReal GENERAL MATRIX.
//
//        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
//        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
//        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
//
//     ON OUTPUT
//
//        WR  AND  WI  CONTAIN THE gReal AND IMAGINARY PARTS,
//        RESPECTIVELY, OF THE EIGENVALUES.  COMPLEX CONJUGATE
//        PAIRS OF EIGENVALUES APPEAR CONSECUTIVELY WITH THE
//        EIGENVALUE HAVING THE POSITIVE IMAGINARY PART FIRST.
//
//        Z  CONTAINS THE gReal AND IMAGINARY PARTS OF THE EIGENVECTORS
//        IF MATZ IS NOT ZERO.  IF THE J-TH EIGENVALUE IS gReal, THE
//        J-TH COLUMN OF  Z  CONTAINS ITS EIGENVECTOR.  IF THE J-TH
//        EIGENVALUE IS COMPLEX WITH POSITIVE IMAGINARY PART, THE
//        J-TH AND (J+1)-TH COLUMNS OF  Z  CONTAIN THE gReal AND
//        IMAGINARY PARTS OF ITS EIGENVECTOR.  THE CONJUGATE OF THIS
//        VECTOR IS THE EIGENVECTOR FOR THE CONJUGATE EIGENVALUE.
//
//        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
//        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
//        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO.
//
//        IV1  AND  FV1  ARE TEMPORARY STORAGE ARRAYS.
//----------------------------------------------------------------------

#define dsign(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void _balanc(int nm, int n, gReal *a, int& low, int& igh, gReal* scale)
{
	int i, j, k, l, m, jj, iexc;
	gReal c, f, g, r, s, b2, radix;
	int noconv;
	radix = 16.0;
	b2 = radix * radix;
	k = 0;
	l = n-1;
	goto L100;
L20: ;
	scale[m] = (gReal)(j+1);
	if ( j == m ) goto L50;
	for ( i = 0; i <= l; i++ )
	{
		f = a[i+j*nm];
		a[i+j*nm] = a[i+m*nm];
		a[i+m*nm] = f;
	}
	for ( i = k; i <= n-1; i++ )
	{
		f = a[j+i*nm];
		a[j+i*nm] = a[m+i*nm];
		a[m+i*nm] = f;
	}
L50: ;
	if ( iexc == 1 ) goto L80;
	else goto L130;
L80: ;
	if ( l == 0 ) goto L280;
	l--;
L100: ;
	for ( jj = 1; jj <= l+1; jj++ )
	{
		j = l+1 - jj;
		for ( i = 0; i <= l; i++ )
		{
			if ( i == j ) goto L110;
			if ( a[j+i*nm] != 0.0 ) goto L120;
L110: ;
		}
		m = l;
		iexc = 1;
		goto L20;
L120: ;
	}
    goto L140;
L130: ;
	k++;
L140: ;
	for ( j = k; j <= l; j++)
	{
		for ( i = k; i <= l; i++ )
		{
			if ( i == j ) goto L150;
			if ( a[i+j*nm] != 0.0 ) goto L170;
L150: ;
		}
		m = k;
		iexc = 2;
		goto L20;
L170: ;
	  }
	for( i = k; i <= l; i++ ) scale[i] = 1.0;
	
L190: ;
	noconv = 0;
	for ( i = k; i <= l; i++ )
	{
		c = 0.0;
		r = 0.0;
		for ( j = k; j <= l; j++ )
		{
			if ( j == i ) goto L200;
			c += fabs(a[j+i*nm]);
			r += fabs(a[i+j*nm]);
L200: ;
		}
		if ( c == 0.0 || r == 0.0 ) goto L270;
		g = r / radix;
		f = 1.0;
		s = c + r;
L210: ;
		if ( c >= g ) goto L220;
		f *= radix;
		c *= b2;
		goto L210;
L220: ;
		g = r * radix;
L230: ;
		if ( c < g ) goto L240;
		f /= radix;
		c /= b2;
		goto L230;
L240: ;
		if ( (c + r) / f >= 0.95 * s ) goto L270;
		g = (gReal)1.0 / f;
		scale[i] *= f;
		noconv = 1;
		for ( j = k; j <= n-1; j++ ) a[i+j*nm] *= g;
		for ( j = 0; j <= l; j++ ) a[j+i*nm] *= f;
L270: ;
	}
	if ( noconv ) goto L190;
L280: ;
	low = k;
	igh = l;
	return;
}

void _balbak(int nm, int n, int low, int igh, gReal* scale, int m, gReal* z)
{
	int i, j, k, ii;
	gReal s;
	if ( m == 0 ) return;
	if ( igh != low ) 
	{
		for ( i = low; i <= igh; i++ )
		{
			s = scale[i];
			for ( j = 0; j <= m-1; j++) z[i+j*nm] *= s;
		}
	}
	for ( ii = 1; ii <= n; ii++)
	{
		i = ii-1;
		if ( i >= low && i <= igh ) goto L140;
		if ( i < low ) i = low - ii;
		k = (int)scale[i] - 1;
		if ( k == i ) goto L140;
		for ( j = 0; j <= m-1; j++ )
		{
			s = z[i+j*nm];
			z[i+j*nm] = z[k+j*nm];
			z[k+j*nm] = s;
		}
L140: ;
	}
	return;
}

void _elmhes(int nm, int n, int low, int igh, gReal* a, int* inf)
{
	int im1, jm1, m, mm1;
	gReal x, y;

	if ( igh < low + 2 ) goto L200;
	for ( m = low+1; m < igh; m++ )
	{
		mm1 = m - 1;
		x = 0.0;
		im1 = m;
		for ( jm1 = m; jm1 <= igh; jm1++)
		{
			if ( fabs(a[jm1+mm1*nm]) <= fabs(x) ) goto L100;
			x = a[jm1+mm1*nm];
            im1 = jm1;
L100: ;
		}
		inf[m] = im1+1;
		if ( im1 == m ) goto L130;
		for ( jm1 = mm1; jm1 < n; jm1++ )
		{
			y = a[im1+jm1*nm];
			a[im1+jm1*nm] = a[m+jm1*nm];
			a[m+jm1*nm] = y;
		}
		for ( jm1 = 0; jm1 <= igh; jm1++ )
		{
			y = a[jm1+im1*nm];
			a[jm1+im1*nm] = a[jm1+m*nm];
			a[jm1+m*nm] = y;
		}
L130: ;
		if ( x == 0.0 ) goto L180;
		
		for ( im1 = m+1; im1 <= igh; im1++ )
		{
			y = a[im1+mm1*nm];
			if (y == 0.0) goto L160;
			y /= x;
			a[im1+mm1*nm] = y;
			for ( jm1 = m; jm1 < n; jm1++ ) a[im1+jm1*nm] -= y * a[m+jm1*nm];
			for ( jm1 = 0; jm1 <= igh; jm1++ ) a[jm1+m*nm] += y * a[jm1+im1*nm];
L160: ;
		}
L180: ;
	}
L200: ;
	return;
}

void _eltran(int nm, int n, int low, int igh, gReal* a, int* inf, gReal* z)
{
	int i, j, kl, mm, m;

	for ( j = 0; j < n; j++ )
	{
		for ( i = 0; i < n; i++ ) z[i+j*nm] = 0.0;
		z[j+j*nm] = 1.0;
	}
	kl = igh - low - 1;
	if ( kl < 1 ) return;
	for ( mm = 1; mm <= kl; mm++ )
	{
		m = igh - mm;
		for ( i = m+1; i <= igh; i++ ) z[i+m*nm] = a[i+(m-1)*nm];
		i = inf[m]-1; 
		if ( i != m )
		{
			for ( j = m; j <= igh; j++ )
			{
				z[m+j*nm] = z[i+j*nm];
				z[i+j*nm] = 0.0;
			}
			z[i+m*nm] = 1.0;
		}
	}
	return;
}

void _hqr(int nm, int n, int low, int igh, gReal* h, gReal* wr, gReal* wi, int& ierr)
{
	int i, j, k, l, ll, m, mm, en, na, itn, its, mp2, enm2;
	gReal p, q, r, s, t, w, x, y, zz, norm, tst1, tst2;
	int notlas;
	ierr = 0;
	norm = 0.0;
	k = 0;
	for ( i = 0; i < n; i++ )
	{
		for ( j = k; j < n; j++ ) norm += fabs(h[i+j*nm]);
		k = i;
		if ( i >= low && i <= igh ) goto L50;
		wr[i] = h[i+i*nm];
		wi[i] = 0.0;
L50: ;
	}	
	en = igh;
	t = 0.0;
	itn = 30*n;
L60: ;
	if ( en < low ) goto L1001;
	its = 0;
	na = en-1;
	enm2 = en-2;	
L70: ;
	for ( ll = low+1; ll <= en+1; ll++ )
	{
		l = en + low+1 - ll;
		if ( l == low ) goto L100;
		s = fabs(h[l-1+(l-1)*nm]) + fabs(h[l+l*nm]);
		if ( s == 0.0 ) s = norm;
		tst1 = s;
		tst2 = tst1 + fabs(h[l+(l-1)*nm]);
		if ( tst2 == tst1 ) goto L100;
	}
L100: ;
	x = h[en+en*nm];
	if ( l == en ) goto L270;
	y = h[na+na*nm];
	w = h[en+na*nm] * h[na+en*nm];
	if ( l == na ) goto L280;
	if ( itn == 0 ) goto L1000;
	if ( its != 10 && its != 20 ) goto L130;
	t += x;
	for ( i = low; i <= en; i++ ) h[i+i*nm] -= x;
	
	s = fabs(h[en+na*nm]) + fabs(h[na+(enm2)*nm]);
	x = (gReal)0.75 * s;
	y = x;
	w = (gReal)-0.4375 * s * s;
L130: ;
	its++;
	itn--;
	
	for ( mm = l; mm <= enm2; mm++ )
	{
		m = enm2 + l - mm;

		zz = h[m+m*nm];
		r = x - zz;
		s = y - zz;
		p = (r * s - w) / h[m+1+m*nm] + h[m+(m+1)*nm];
		q = h[m+1+(m+1)*nm] - zz - r - s;
		r = h[m+2+(m+1)*nm];
		s = fabs(p) + fabs(q) + fabs(r);
		p /= s;
		q /= s;
		r /= s;
		if ( m == l ) goto L150;
		tst1 = fabs(p)*(fabs(h[m-1+(m-1)*nm]) + fabs(zz) + fabs(h[m+1+(m+1)*nm]));
		tst2 = tst1 + fabs(h[m+(m-1)*nm])*(fabs(q) + fabs(r));
		if ( tst2 == tst1 ) goto L150;
	}
L150: ;
	mp2 = m+2;
	for ( i = mp2; i <= en; i++ )
	{
		h[i+(i-2)*nm] = 0.0;
		if ( i == mp2 ) goto L160;
		h[i+(i-3)*nm] = 0.0;
L160: ;
	}
	for ( k = m; k <= na; k++)
	{
		notlas = (k != na);
		if ( k == m ) goto L170;
		p = h[k+(k-1)*nm];
		q = h[k+1+(k-1)*nm];
		r = 0.0;
		if ( notlas ) r = h[k+2+(k-1)*nm];
		x = fabs(p) + fabs(q) + fabs(r);
		if ( x == 0.0 ) goto L260;
		p /= x;
		q /= x;
		r /= x;
L170: ;
		s = dsign(sqrt(p*p+q*q+r*r), p);
		if ( k == m ) goto L180;
		h[k+(k-1)*nm] = -s * x;
		goto L190;
L180: ;
		if ( l != m ) h[k+(k-1)*nm] = -h[k+(k-1)*nm];
L190: ;
		p += s;
		x = p / s;
		y = q / s;
		zz = r / s;
		q /= p;
		r /= p;
		if ( notlas ) goto L225;
		for ( j = k; j <= en; j++ )
		{
			p = h[k+j*nm] + q * h[k+1+j*nm];
			h[k+j*nm] -= p * x;
			h[k+1+j*nm] -= p * y;
		}
		j = _MIN(en, k+3);
		for ( i = l; i <= j; i++ )
		{
            p = x * h[i+k*nm] + y * h[i+(k+1)*nm];
            h[i+k*nm] -= p;
            h[i+(k+1)*nm] -= p * q;
		}
		goto L255;
L225: ;
	 	for ( j = k; j <= en; j++ )
		{
            p = h[k+j*nm] + q * h[k+1+j*nm] + r * h[k+2+j*nm];
            h[k+j*nm] -= p * x;
            h[k+1+j*nm] -= p * y;
            h[k+2+j*nm] -= p * zz;
		}
		j = _MIN(en,k+3);
		for ( i = l; i <= j; i++ )
		{
            p = x * h[i+k*nm] + y * h[i+(k+1)*nm] + zz * h[i+(k+2)*nm];
            h[i+k*nm] -= p;
            h[i+(k+1)*nm] -= p * q;
            h[i+(k+2)*nm] -= p * r;
		}
L255: ;
L260: ;
	}
	goto L70;
L270: ;
	wr[en] = x + t;
	wi[en] = 0.0;
	en = na;
	goto L60;
L280: ;
	p = (y - x) / (gReal)2.0;
	q = p * p + w;
	zz = sqrt(fabs(q));
	x += t;
	if ( q < 0.0 ) goto L320;
	zz = p + dsign(zz,p);
	wr[na] = x + zz;
	wr[en] = wr[na];
	if ( zz != 0.0 ) wr[en] = x - w / zz;
	wi[na] = 0.0;
	wi[en] = 0.0;
	goto L330;
L320: ;
	wr[na] = x + p;
	wr[en] = x + p;
	wi[na] = zz;
	wi[en] = -zz;
L330: ;
	en = enm2;
	goto L60;
L1000: ;
	ierr = en+1;
L1001: ;
	return;
}

void _cdiv(gReal ar, gReal ai, gReal br, gReal bi, gReal& cr, gReal& ci)
{
	gReal is, ars, ais, brs, bis;
	is = (gReal)1.0 / ( fabs(br) + fabs(bi) );
	ars = ar * is;
	ais = ai * is;
	brs = br * is;
	bis = bi * is;
	is = (gReal)1.0 / ( brs * brs + bis * bis );
	cr = (ars * brs + ais * bis) * is;
	ci = (ais * brs - ars * bis) * is;
	return;
}

void _hqr2(int nm, int n, int low, int igh, gReal* h, gReal* wr, gReal* wi, gReal* z, int& ierr)
{
	int i, j, k, l, m, en, na, ii, jj, ll, mm, nn, itn, its, mp2, enm2;
	gReal p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, norm, tst1, tst2;
	int notlas;
	ierr = 0;
	norm = 0.0;
	k = 0;
	for ( i = 0; i < n; i++ )
	{
		for ( j = k; j <= n-1; j++ ) norm += fabs(h[i+j*nm]);
		
		k = i;
		if ( i >= low && i <= igh ) goto L50;
		wr[i] = h[i+i*nm];
		wi[i] = 0.0;
L50: ;
	}
	en = igh;
	t = 0.0;
	itn = 30*n;
L60: ;
	if ( en < low ) goto L340;
	its = 0;
	na = en - 1;
	enm2 = en - 1;
L70: ;
	for ( ll = low+1; ll <= en+1; ll++ )
	{
		l = en + low+1 - ll;
		if ( l == low ) goto L100;
		s = fabs(h[l-1+(l-1)*nm]) + fabs(h[l+l*nm]);
		if ( s == 0.0 ) s = norm;
		tst1 = s;
		tst2 = tst1 + fabs(h[l+(l-1)*nm]);
		if ( tst2 == tst1 ) goto L100;
	}
L100: ;
	x = h[en+en*nm];
	if ( l == en ) goto L270;
	y = h[na+na*nm];
	w = h[en+na*nm] * h[na+en*nm];
	if ( l == na ) goto L280;
	if ( itn == 0 ) goto L1000;
	if ( its != 10 && its != 20 ) goto L130;
	t += x;
	for (  i = low; i <= en; i++ ) h[i+i*nm] -= x;
	
	s = fabs(h[en+na*nm]) + fabs(h[na+(enm2-1)*nm]);
	x = (gReal)0.75 * s;
	y = x;
	w = (gReal)-0.4375 * s * s;
L130: ;
	its++;
	itn--;
	for ( mm = l+1; mm <= enm2; mm++ )
	{
		m = enm2 + l - mm;
		zz = h[m+m*nm];
		r = x - zz;
		s = y - zz;
		p = (r * s - w) / h[m+1+m*nm] + h[m+(m+1)*nm];
		q = h[m+1+(m+1)*nm] - zz - r - s;
		r = h[m+2+(m+1)*nm];
		s = fabs(p) + fabs(q) + fabs(r);
		p /= s;
		q /= s;
		r /= s;
		if ( m == l ) goto L150;
		tst1 = fabs(p)*(fabs(h[m-1+(m-1)*nm]) + fabs(zz) + fabs(h[m+1+(m+1)*nm]));
		tst2 = tst1 + fabs(h[m+(m-1)*nm])*(fabs(q) + fabs(r));
		if ( tst2 == tst1 ) goto L150;
	}
L150: ;
	mp2 = m+2;
	for ( i = mp2; i <= en; i++ )
	{
		h[i+(i-2)*nm] = 0.0;
		if ( i == mp2 ) goto L160;
		h[i+(i-3)*nm] = 0.0;
L160: ;
	}
	for ( k = m; k <= na; k++ )
	{
		notlas = (k != na);
		if ( k == m ) goto L170;
		p = h[k+(k-1)*nm];
		q = h[k+1+(k-1)*nm];
		r = 0.0;
		if ( notlas ) r = h[k+2+(k-1)*nm];
		x = fabs(p) + fabs(q) + fabs(r);
		if ( x == 0.0 ) goto L260;
		p /= x;
		q /= x;
		r /= x;
L170: ;
		s = dsign(sqrt(p*p+q*q+r*r),p);
		if ( k == m ) goto L180;
		h[k+(k-1)*nm] = -s * x;
		goto L190;
L180: ;
		if ( l != m ) h[k+(k-1)*nm] = -h[k+(k-1)*nm];
L190: ;
		p += s;
		x = p / s;
		y = q / s;
		zz = r / s;
		q /= p;
		r /= p;
		if ( notlas ) goto L225;
		for ( j = k; j <= n-1; j++ )
		{
			p = h[k+j*nm] + q * h[k+1+j*nm];
			h[k+j*nm] -= p * x;
			h[k+1+j*nm] -= p * y;
		}
		j = _MIN(en,k+3);
		for ( i = 0; i <= j; i++ )
		{
			p = x * h[i+k*nm] + y * h[i+(k+1)*nm];
			h[i+k*nm] -= p;
			h[i+(k+1)*nm] -= p * q;
		}
		for ( i = low; i <= igh; i++ )
		{
			p = x * z[i+k*nm] + y * z[i+(k+1)*nm];
			z[i+k*nm] -= p;
			z[i+(k+1)*nm] -= p * q;
		}
		goto L255;
L225: ;
		for ( j = k; j <= n-1; j++ )
		{
			p = h[k+j*nm] + q * h[k+1+j*nm] + r * h[k+2+j*nm];
			h[k+j*nm] -= p * x;
			h[k+1+j*nm] -= p * y;
			h[k+2+j*nm] -= p * zz;
		}
		j = _MIN(en,k+3);
		for ( i = 0; i <= j; i++ )
		{
			p = x * h[i+k*nm] + y * h[i+(k+1)*nm] + zz * h[i+(k+2)*nm];
			h[i+k*nm] -= p;
			h[i+(k+1)*nm] -= p * q;
			h[i+(k+2)*nm] -= p * r;
		}
		for ( i = low; i <= igh; i++ )
		{
			p = x * z[i+k*nm] + y * z[i+(k+1)*nm] + zz * z[i+(k+2)*nm];
			z[i+k*nm] -= p;
			z[i+(k+1)*nm] -= p * q;
			z[i+(k+2)*nm] -= p * r;
		}
L255: ;
L260: ;
	}
	goto L70;
L270: ;
	h[en+en*nm] = x + t;
	wr[en] = h[en+en*nm];
	wi[en] = 0.0;
	en = na;
	goto L60;
L280: ;
	p = (y - x) / (gReal)2.0;
	q = p * p + w;
	zz = sqrt(fabs(q));
	h[en+en*nm] = x + t;
	x = h[en+en*nm];
	h[na+na*nm] = y + t;
	if ( q < 0.0 ) goto L320;
	zz = p + dsign(zz,p);
	wr[na] = x + zz;
	wr[en] = wr[na];
	if ( zz != 0.0 ) wr[en] = x - w / zz;
	wi[na] = 0.0;
	wi[en] = 0.0;
	x = h[en+na*nm];
	s = fabs(x) + fabs(zz);
	p = x / s;
	q = zz / s;
	r = sqrt(p*p+q*q);
	p /= r;
	q /= r;
	for ( j = na; j < n; j++ )
	{
		zz = h[na+j*nm];
		h[na+j*nm] = q * zz + p * h[en+j*nm];
		h[en+j*nm] = q * h[en+j*nm] - p * zz;
	}
	for ( i = 0; i <= en; i++ )
	{
		zz = h[i+na*nm];
		h[i+na*nm] = q * zz + p * h[i+en*nm];
		h[i+en*nm] = q * h[i+en*nm] - p * zz;
	}
	for ( i = low; i <= igh; i++ )
	{
		zz = z[i+na*nm];
		z[i+na*nm] = q * zz + p * z[i+en*nm];
		z[i+en*nm] = q * z[i+en*nm] - p * zz;
	}
	goto L330;
L320: ;
	wr[na] = x + p;
	wr[en] = x + p;
	wi[na] = zz;
	wi[en] = -zz;
L330: ;
	en = enm2-1;
	goto L60;
L340: ;
	if ( norm == 0.0 ) goto L1001;
	for ( nn = 1; nn <= n; nn++ )
	{
		en = n - nn;
		p = wr[en];
		q = wi[en];
		na = en-1;
		if ( q < 0 ) goto L710;
		else if ( q == 0 ) goto L600;
		else goto L800;
L600: ;
		m = en;
		h[en+en*nm] = 1.0;
		if ( na+1 == 0 ) goto L800;
		for ( ii = 1; ii <= na+1; ii++ )
		{
			i = en - ii ;
			w = h[i+i*nm] - p;
			r = 0.0;
            for ( j = m; j <= en; j++ ) r += h[i+j*nm] * h[j+en*nm];
			
			if ( wi[i] >= 0.0 ) goto L630;
            zz = w;
            s = r;
            goto L700;
L630: ;   
			m = i;
            if ( wi[i] != 0.0 ) goto L640;
            t = w;
            if ( t != 0.0 ) goto L635;
			tst1 = norm;
			t = tst1;
L632: ;
			t = (gReal)0.01 * t;
			tst2 = norm + t;
			if ( tst2 > tst1 ) goto L632;
L635: ;
			h[i+en*nm] = -r / t;
            goto L680;
L640: ;
			x = h[i+(i+1)*nm];
            y = h[i+1+i*nm];
            q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
            t = (x * s - zz * r) / q;
            h[i+en*nm] = t;
            if ( fabs(x) <= fabs(zz) ) goto L650;
            h[i+1+en*nm] = (-r - w * t) / x;
            goto L680;
L650: ;
			h[i+1+en*nm] = (-s - y * t) / zz;
L680: ;
			t = fabs(h[i+en*nm]);
            if ( t == 0.0 ) goto L700;
            tst1 = t;
            tst2 = tst1 + (gReal)1.0/tst1;
            if ( tst2 > tst1 ) goto L700;
            for ( j = i; j <= en; j++ ) h[j+en*nm] /= t;			
L700: ;
		}
		goto L800;
L710: ;
		m = na;
		if ( fabs(h[en+na*nm]) <= fabs(h[na+en*nm]) ) goto L720;
		h[na+na*nm] = q / h[en+na*nm];
		h[na+en*nm] = -(h[en+en*nm] - p) / h[en+na*nm];
		goto L730;
L720: ;
		_cdiv(0.0, -h[na+en*nm], h[na+na*nm]-p, q, h[na+na*nm], h[na+en*nm] );
L730: ;
		h[en+na*nm] = 0.0;
		h[en+en*nm] = 1.0;
		enm2 = na;
		if ( enm2 == 0 ) goto L800;
		for ( ii = 1; ii <= enm2; ii++ )
		{
			i = na - ii;
			w = h[i+i*nm] - p;
			ra = 0.0;
			sa = 0.0;
			for ( j = m; j <= en; j++ )
			{
				ra += h[i+j*nm] * h[j+na*nm];
				sa += h[i+j*nm] * h[j+en*nm];
			}
            if ( wi[i] >= 0.0 ) goto L770;
            zz = w;
            r = ra;
            s = sa;
            goto L795;
L770: ;
			m = i;
            if ( wi[i] != 0.0 ) goto L780;
			_cdiv( -ra, -sa, w, q, h[i+na*nm], h[i+en*nm]);
            goto L790;
L780: ;
			x = h[i+(i+1)*nm];
            y = h[i+1+i*nm];
            vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
            vi = (wr[i] - p) * (gReal)2.0 * q;
            if ( vr != 0.0 || vi != 0.0 ) goto L784;
			tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
			vr = tst1;
L783: ;
			vr = (gReal)0.01 * vr;
			tst2 = tst1 + vr;
			if ( tst2 > tst1 ) goto L783;
L784: ;
			_cdiv(x*r-zz*ra+q*sa, x*s-zz*sa-q*ra, vr, vi, h[i+na*nm], h[i+en*nm]);
			if ( fabs(x) <= fabs(zz) + fabs(q) ) goto L785;
			h[i+1+na*nm] = (-ra - w * h[i+na*nm] + q * h[i+en*nm]) / x;
			h[i+1+en*nm] = (-sa - w * h[i+en*nm] - q * h[i+na*nm]) / x;
			goto L790;
L785: ;
			_cdiv(-r-y*h[i+na*nm], -s-y*h[i+en*nm], zz, q, h[i+1+na*nm], h[i+1+en*nm]);
L790: ;
			t = _MAX(fabs(h[i+na*nm]), fabs(h[i+en*nm]));
			if ( t == 0.0 ) goto L795;
			tst1 = t;
			tst2 = tst1 + (gReal)1.0/tst1;
			if ( tst2 > tst1 ) goto L795;
			for ( j = i; j <= en; j++ )
			{
				h[j+na*nm] /= t;
				h[j+en*nm] /= t;
			}
L795: ;
		}
L800: ;
	}
	for ( i = 0; i < n; i++ )
	{
		if ( i >= low && i <= igh ) goto L840;
		for ( j = i; j < n; j++ ) z[i+j*nm] = h[i+j*nm];
L840: ;
	}
	for ( jj = low+1; jj <= n; jj++ )
	{
		j = n + low - jj;
		m = _MIN(j,igh);
		for ( i = low; i <= igh; i++ )
		{
			zz = 0.0;
			for ( k = low; k <= m; k++ ) zz += z[i+k*nm] * h[k+j*nm];
			
			z[i+j*nm] = zz;
		}
	}
	goto L1001;
L1000: ;
	ierr = en+1;
L1001: ;
	return;
}

void _rg(int nm, int n, gReal* a, gReal* wr, gReal* wi, int matz, gReal* z, int* iv1, gReal* fv1, int& ierr)
{
	int is1, is2;
	if ( n > nm )
	{
		ierr = 10 * n;
		return;
	}
	_balanc(nm,n,a,is1,is2,fv1);
	_elmhes(nm,n,is1,is2,a,iv1);
	if ( matz == 0 )
	{
		_hqr(nm,n,is1,is2,a,wr,wi,ierr);
		return;
	}
	_eltran(nm, n, is1, is2, a, iv1, z);
	_hqr2(nm, n, is1, is2, a, wr, wi, z, ierr);
	if ( ierr != 0 ) return;
	_balbak(nm,n,is1,is2,fv1,n,z);
	return;
}

//////////////////////////////////////////////////////////////////////////////
//
// End of EISPACK
//
//////////////////////////////////////////////////////////////////////////////

RMatrix Eig(RMatrix m)
{
	if ( m.row != m.col ) std::cerr << "RMatrix Eig : not sqare" << std::endl;
	
	int i, ierr, matz = 0, n = m.row;
	static IMatrix _iv1_Eig;
	static RMatrix _fv1_Eig;

	if ( _iv1_Eig.row < n ) _iv1_Eig.ReNew(n,1);
	if ( _fv1_Eig.row < n ) _fv1_Eig.ReNew(n,1);
	RMatrix re(n,2);
	
	_rg(n, n, m.element, re.element, (re.element+n), matz, 0, _iv1_Eig.element, _fv1_Eig.element, ierr);
	
	gReal sum = 0.0;
	for ( i = 0; i < n; i++ ) sum += re.element[i+n] * re.element[i+n];
	if ( sum < _EPS ) re.col = 1;
	
	return re;
}

void Eig(RMatrix &re, RMatrix &m)
{
	if ( m.row != m.col ) std::cerr << "RMatrix Eig : not sqare" << std::endl;
	
	int ierr, matz = 0, n = m.row;
	static IMatrix _iv1_Eig;
	static RMatrix _fv1_Eig;

	if ( _iv1_Eig.row < n ) _iv1_Eig.ReNew(n,1);
	if ( _fv1_Eig.row < n ) _fv1_Eig.ReNew(n,1);
	re.ReNew(n,2);
	
	_rg(n, n, m.element, re.element, (re.element+n), matz, 0, _iv1_Eig.element, _fv1_Eig.element, ierr);
}

void Eig(RMatrix m, RMatrix &v, RMatrix &d)
{
//	int i;
	int n = m.row;
	static IMatrix _iv1_Eig;
	static RMatrix _fv1_Eig;

	if ( n != m.col ) std::cerr << "RMatrix Eig : not sqare";
	
	if ( _iv1_Eig.row < n ) _iv1_Eig.ReNew(n,1);
	if ( _fv1_Eig.row < n ) _fv1_Eig.ReNew(n,1);
	
	int ierr, matz = 1;
	
	if ( v.row*v.col != n*n )
	{
		delete [] v.element;
		v.element = new gReal [n*n];
	}
	v.row = v.col = n;

	if ( d.row*d.col != 2*n )
	{
		delete [] d.element;
		d.element = new gReal [2*n];
	}
	d.row = n;
	d.col = 2;

	_rg(n, n, m.element, d.element, d.element+n, matz, v.element, _iv1_Eig.element, _fv1_Eig.element, ierr);
	
//	gReal sum = 0.0;
//	for ( i = 0; i < n; i++ ) sum += d.element[i+n] * d.element[i+n];
//	if ( sum < _EPS ) d.col = 1;
}

gReal _pythag(gReal a, gReal b)
{
	gReal absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if ( absa > absb ) return absa * sqrt((gReal)1.0 + (absb / absa) * (absb / absa));
	else return ( absb == 0 ? 0 : absb * sqrt((gReal)1.0 + (absa / absb) * (absa / absb)) );
}

// singular value decomposition in 'numerical recipe'
// modified with EISPACK 'dsvd.f'
void _svdcmp(int m, int n, gReal *a, gReal *w, gReal *v, int matu, int matv, gReal *rv1)
{
	int flag, i, its, j, jj, k, l, nm;
	gReal anorm, c, f, g, h, s, scale, x, y, z;
	//gReal *rv1 = new gReal [n];
	g = scale = anorm = 0.0; // Householder reduction to bidiagonal form
	for ( i = 0; i < n; i++ ) 
	{
		l = i+1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if ( i < m ) 
		{
			for ( k = i; k < m; k++ ) scale += fabs(a[k+m*i]);
			if ( scale ) 
			{
				for ( k = i; k < m; k++ ) 
				{
					a[k+m*i] /= scale;
					s += a[k+m*i] * a[k+m*i];
				}
				f = a[i+m*i];
				g = -dsign(sqrt(s), f);
				h = f * g - s;
				a[i+m*i] = f - g;
				for ( j = l; j < n; j++ ) 
				{
					for ( s = 0.0, k = i; k < m; k++ ) s += a[k+m*i] * a[k+m*j];
					f = s / h;
					for ( k = i; k < m; k++ ) a[k+m*j] += f * a[k+m*i];
				}
				for ( k = i; k < m; k++ ) a[k+m*i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if ( i < m && i != n-1 )
		{
			for ( k = l; k < n; k++ ) scale += fabs(a[i+m*k]);
			if ( scale ) 
			{
				for ( k = l; k < n; k++ ) 
				{
					a[i+m*k] /= scale;
					s += a[i+m*k] * a[i+m*k];
				}
				f = a[i+m*l];
				g = -dsign(sqrt(s), f);
				h = f * g - s;
				a[i+m*l] = f - g; 
				for ( k = l; k < n; k++ ) rv1[k] = a[i+m*k] / h;
				for ( j = l; j < m; j++ ) 
				{
					for ( s = 0.0, k = l; k < n; k++ ) s += a[j+m*k] * a[i+m*k];
					for ( k = l ; k < n; k++ ) a[j+m*k] += s * rv1[k];
				}
				for ( k = l ; k < n; k++ ) a[i+m*k] *= scale;
			}
		}
		anorm = _MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	if ( matv ) 
	{
		for ( i = n-1 ; i >= 0; i-- ) // Accumulation of right-system transformations
		{
			if ( i < n ) 
			{
				if ( g )
				{
					for ( j = l; j < n; j++ ) // gReal devision to avoid possible underflow
						v[j+n*i] = (a[i+m*j] / a[i+m*l]) / g;
					for ( j = l; j < n; j++ ) 
					{
						for ( s = 0.0, k = l; k < n; k++ ) s += a[i+m*k] * v[k+n*j];
						for ( k = l; k < n; k++ ) v[k+n*j] += s * v[k+n*i];
					}
				}
				for ( j = l; j < n; j++ ) v[i+n*j] = v[j+n*i] = 0.0;
			}
			v[i+n*i] = 1.0;
			g = rv1[i];
			l = i;
		}
	}
	if ( matu )
	{
		for ( i = _MIN(m-1, n-1); i >= 0; i-- ) // Accumulation of left-system transformations
		{
			l = i + 1;
			g = w[i];
			for ( j = l; j < n; j++ ) a[i+m*j] = 0.0;
			if ( g )
			{
				g = (gReal)1.0 / g;
				for ( j = l; j < n; j++ ) 
				{
					for ( s = 0.0, k = l; k < m; k++ ) s += a[k+m*i] * a[k+m*j];
					f = (s / a[i+m*i]) * g;
					for ( k = i; k < m; k++ ) a[k+m*j] += f * a[k+m*i];
				}
				for ( j = i; j < m; j++ ) a[j+m*i] *= g;
			} else for ( j = i; j < m; j++ ) a[j+m*i] = 0.0;
			++a[i+m*i];
		}
	}
	for ( k = n-1; k >= 0; k-- ) // Diagonalization of the bidiagonal form: Loop over
	{
		for ( its = 1; its <= 30; its++ ) // singular values, and over allowed iterations
		{
			flag = 1;
			for ( l = k; l >= 0; l-- ) // Test for splitting
			{
				nm = l-1; // Note that rv1(1) is always zero
				if ( (fabs(rv1[l]) + anorm) == anorm )
				{
					flag = 0;
					break;
				}
				if ( (fabs(w[nm]) + anorm) == anorm ) break;
			}
			if ( flag )
			{
				c = 0.0;
				s = 1.0;
				for ( i = l; i <= k; i++) 
				{
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ( (fabs(f) + anorm) == anorm ) break;
					g = w[i];
					h = _pythag(f, g);
					w[i] = h;
					h = (gReal)1.0 / h;
					c = g * h;
					s = -f * h;
					if ( matu ) 
					{
						for(j = 0; j < m; j++) 
						{
							y = a[j+m*nm];
							z = a[j+m*i];
							a[j+m*nm] = y * c + z * s;
							a[j+m*i] = z * c - y * s;
						}
					}
				}
			}
			z = w[k];
			if ( l == k ) // Convergence
			{
				if ( z < 0.0 ) // Singular value is made nonnegative
				{
					w[k] = -z;
					if ( matv ) for(j = 0; j < n; j++) v[j+n*k] = -v[j+n*k];
				}
				break;
			}
			if ( its == 30 ) std::cerr << "svdcmp : no convergence in 30 iterations" << std::endl;
			x = w[l]; // Shift from bottom 2-by-2 minor
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z) + (g-h)*(g+h)) / ((gReal)2.0 * h * y);
			g = _pythag(f, 1.0);
			f = ((x-z)*(x+z) + h * ((y / (f + dsign(g, f))) - h)) / x;
			c = s = 1.0; // Next QR transformation
			for ( j = l; j <= nm; j++ )
			{
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = _pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				if ( matv ) 
				{
					for ( jj = 0; jj < n; jj++ )
					{
						x = v[jj+n*j];
						z = v[jj+n*i];
						v[jj+n*j] = x * c + z * s;
						v[jj+n*i] = z * c - x * s;
					}
				}
				z = _pythag(f, h);
				w[j] = z; // Rotation can be arbitary if z = 0
				if(z) 
				{
					z = (gReal)1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				if ( matu ) 
				{
					for(jj = 0; jj < m; jj++) 
					{
						y = a[jj+m*j];
						z = a[jj+m*i];
						a[jj+m*j] = y * c + z * s;
						a[jj+m*i] = z * c - y * s;
					}
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	//delete [] rv1;
}

RMatrix SVD(const RMatrix &m)
{
	static RMatrix _rv_SVD, _U_SVD;
	
	// debugged by junggon, 8/28/2007
	//if ( m.row < m.col )
	//{
	//	_U_SVD.ReNew(m.col, m.row);
	//	int i = 0, r = m.row, c;
	//	gReal *a = _U_SVD.element, *at;
	//	while ( r-- )
	//	{
	//		at = m.element + (i++);
	//		c = m.col;
	//		while ( c-- )
	//		{
	//			*(a++) = *at;
	//			at += m.row;
	//		}
	//	}
	//} else _U_SVD = m;
	
	_U_SVD = m;

	RMatrix S(_U_SVD.col, 1);
	
	if ( _rv_SVD.row < _U_SVD.col ) _rv_SVD.ReNew(_U_SVD.col, 1);
	_svdcmp(_U_SVD.row, _U_SVD.col, _U_SVD.element, S.element, NULL, 0, 0, _rv_SVD.element);

	return S;
}

void SVD(const RMatrix &M, RMatrix &U, RMatrix &S, RMatrix &V)
{
	static RMatrix _rv_SVD, _U_SVD;

	// debugged by junggon, 8/28/2007
	//if ( M.row > M.col ) U = ~M;		
	//else U = M;

	U = M;

	if ( M.col > _rv_SVD.row ) _rv_SVD.ReNew(M.col,1);
	S.ReNew(M.col,1);
	V.ReNew(M.col, M.col); 
	_svdcmp(M.row, M.col, U.element, S.element, V.element, 1, 1, _rv_SVD.element);
}

RMatrix SVD33(const RMatrix &M)
{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( M.row != 3 || M.col != 3 ) std::cerr << "size is not compatible" << std::endl;
#endif	
	gReal S[3];
	SVD33(M.element, S);
	return RMatrix(3,1,S);
}

void SVD33(gReal *M, gReal *S)
{
	gReal m, detC, q, p, phi, co, si, sqrt_p, sqrt_3, tmp, eps = (gReal)_EPS;
	gReal A[9], C[9], E[3];

	// if M=0, return [0;0;0]
	if ( fabs(M[0])<eps && fabs(M[1])<eps && fabs(M[2])<eps && fabs(M[3])<eps && fabs(M[4])<eps && fabs(M[5])<eps && fabs(M[6])<eps && fabs(M[7])<eps && fabs(M[8])<eps ) {
		S[0] = 0.; S[1] = 0.; S[2] = 0.; 
	}

	// A = ~M*M
	A[0] = M[0]*M[0] + M[1]*M[1] + M[2]*M[2];
	A[1] = M[0]*M[3] + M[1]*M[4] + M[2]*M[5];
	A[2] = M[0]*M[6] + M[1]*M[7] + M[2]*M[8];
	A[4] = M[3]*M[3] + M[4]*M[4] + M[5]*M[5];
	A[5] = M[3]*M[6] + M[4]*M[7] + M[5]*M[8];
	A[8] = M[6]*M[6] + M[7]*M[7] + M[8]*M[8];
	A[3] = A[1];
	A[6] = A[2];
	A[7] = A[5];

	// E(3,1) = the eigenvalues of the symmetric A (closed-form solution!)
	m = (A[0]+A[4]+A[8])/(gReal)3.;
	C[0] = A[0]-m; C[1] = A[1]; C[2] = A[2]; C[3] = A[3]; C[4] = A[4]-m; C[5] = A[5]; C[6] = A[6]; C[7] = A[7]; C[8] = A[8]-m; // C = A-m*eye(3,3)
	detC = -(C[2]*C[4]*C[6]) + C[1]*C[5]*C[6] + C[2]*C[3]*C[7] - C[0]*C[5]*C[7] - C[1]*C[3]*C[8] + C[0]*C[4]*C[8];
	q = (gReal)0.5*detC;
	p = (C[0]*C[0] + C[1]*C[1] + C[2]*C[2] + C[3]*C[3] + C[4]*C[4] + C[5]*C[5] + C[6]*C[6] + C[7]*C[7] + C[8]*C[8])/(gReal)6.;
	if ( p == 0. && q == 0. ) { //if ( fabs(p) < eps && fabs(q) < eps ) {
		E[0] = m; E[1] = m; E[2] = m;
	} else {
		if ( q == 0. ) {
			phi = (gReal)3.141592653589793/(gReal)6.;
		} else {
			phi = atan(sqrt(fabs(p*p*p-q*q))/q) / (gReal)3.;	// in theroy p^3-q^2>=0, but make sure it by using fabs()
		}
		co = cos(phi);
		si = sin(phi);
		sqrt_p = sqrt(p);
		sqrt_3 = sqrt((gReal)3.);
		if ( q >= 0 ) {
			E[0] = m + (gReal)2. * sqrt_p * co;
			E[1] = m - sqrt_p * ( co + sqrt_3 * si );
			E[2] = m - sqrt_p * ( co - sqrt_3 * si );
		} else {
			E[0] = m - (gReal)2. * sqrt_p * co;
			E[1] = m + sqrt_p * ( co - sqrt_3 * si );
			E[2] = m + sqrt_p * ( co + sqrt_3 * si );
		}
	}

	// S(3,1) = the singular values of M (closed-form solution!)
	S[0] = sqrt(fabs(E[0]));	// though E[i]>=0 in theory, make sure it even for singular M by using fabs().
	S[1] = sqrt(fabs(E[1]));
	S[2] = sqrt(fabs(E[2]));

	// sort S in descending order
	if ( S[0] < S[1] ) { 
		tmp = S[0]; S[0] = S[1]; S[1] = tmp; 
	}
	if ( S[0] < S[2] ) { 
		tmp = S[0]; S[0] = S[2]; S[2] = tmp; 
	}
	if ( S[1] < S[2] ) { 
		tmp = S[1]; S[1] = S[2]; S[2] = tmp; 
	}
}

void SVD33(const RMatrix &M, RMatrix &U, RMatrix &S, RMatrix &V, gReal tol)
{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( M.row != 3 || M.col != 3 ) std::cerr << "size is not compatible" << std::endl;
#endif	
	S.ReNew(3,1); U.ReNew(3,3); V.ReNew(3,3);
	SVD33(M.element, U.element, S.element, V.element, tol);
}

void SVD33(gReal *M, gReal *U, gReal *S, gReal *V, gReal tol)
{
	int idx, rnk;
	gReal m, detC, q, p, phi, co, si, sqrt_p, sqrt_3, tmp, detM, tmpVec[3], eps = (gReal)_EPS;
	gReal A[9], C[9], E[3], pinvS[3];

	if ( tol <= eps ) { tol = eps; }

	// if M=0, S = [0;0;0], U = V = Eye(3,3)
	if ( fabs(M[0])<eps && fabs(M[1])<eps && fabs(M[2])<eps && fabs(M[3])<eps && fabs(M[4])<eps && fabs(M[5])<eps && fabs(M[6])<eps && fabs(M[7])<eps && fabs(M[8])<eps ) {
		S[0] = 0.; S[1] = 0.; S[2] = 0.;
		U[0] = 1.; U[1] = 0.; U[2] = 0.; U[3] = 0.; U[4] = 1.; U[5] = 0.; U[6] = 0.; U[7] = 0.; U[8] = 1.;
		V[0] = 1.; V[1] = 0.; V[2] = 0.; V[3] = 0.; V[4] = 1.; V[5] = 0.; V[6] = 0.; V[7] = 0.; V[8] = 1.;
		return;
	}

	// detM = det(M)
	detM = -(M[2]*M[4]*M[6]) + M[1]*M[5]*M[6] + M[2]*M[3]*M[7] - M[0]*M[5]*M[7] - M[1]*M[3]*M[8] + M[0]*M[4]*M[8];

	// A = ~M*M
	A[0] = M[0]*M[0] + M[1]*M[1] + M[2]*M[2];
	A[1] = M[0]*M[3] + M[1]*M[4] + M[2]*M[5];
	A[2] = M[0]*M[6] + M[1]*M[7] + M[2]*M[8];
	A[4] = M[3]*M[3] + M[4]*M[4] + M[5]*M[5];
	A[5] = M[3]*M[6] + M[4]*M[7] + M[5]*M[8];
	A[8] = M[6]*M[6] + M[7]*M[7] + M[8]*M[8];
	A[3] = A[1];
	A[6] = A[2];
	A[7] = A[5];

	// E(3,1) = the eigenvalues of the symmetric A (closed-form solution!)
	m = (A[0]+A[4]+A[8])/(gReal)3.;
	C[0] = A[0]-m; C[1] = A[1]; C[2] = A[2]; C[3] = A[3]; C[4] = A[4]-m; C[5] = A[5]; C[6] = A[6]; C[7] = A[7]; C[8] = A[8]-m; // C = A-m*eye(3,3)
	detC = -(C[2]*C[4]*C[6]) + C[1]*C[5]*C[6] + C[2]*C[3]*C[7] - C[0]*C[5]*C[7] - C[1]*C[3]*C[8] + C[0]*C[4]*C[8];
	q = (gReal)0.5*detC;
	p = (C[0]*C[0] + C[1]*C[1] + C[2]*C[2] + C[3]*C[3] + C[4]*C[4] + C[5]*C[5] + C[6]*C[6] + C[7]*C[7] + C[8]*C[8])/(gReal)6.;
	if ( p == 0. && q == 0. ) { //if ( fabs(p) < eps && fabs(q) < eps ) {
		E[0] = m; E[1] = m; E[2] = m;
	} else {
		if ( q == 0. ) {
			phi = (gReal)3.141592653589793/(gReal)6.;
		} else {
			phi = atan(sqrt(fabs(p*p*p-q*q))/q) / (gReal)3.;	// in theroy p^3-q^2>=0, but make sure it by using fabs()
		}
		co = cos(phi);
		si = sin(phi);
		sqrt_p = sqrt(p);
		sqrt_3 = sqrt((gReal)3.);
		if ( q >= 0 ) {
			E[0] = m + (gReal)2. * sqrt_p * co;
			E[1] = m - sqrt_p * ( co + sqrt_3 * si );
			E[2] = m - sqrt_p * ( co - sqrt_3 * si );
		} else {
			E[0] = m - (gReal)2. * sqrt_p * co;
			E[1] = m + sqrt_p * ( co - sqrt_3 * si );
			E[2] = m + sqrt_p * ( co + sqrt_3 * si );
		}
	}

	// make sure E[i]>=0 by applying fabs()
	E[0] = fabs(E[0]);
	E[1] = fabs(E[1]);
	E[2] = fabs(E[2]);

	// sort E in descending order
	if ( E[0] < E[1] ) { 
		tmp = E[0]; E[0] = E[1]; E[1] = tmp; 
	}
	if ( E[0] < E[2] ) { 
		tmp = E[0]; E[0] = E[2]; E[2] = tmp; 
	}
	if ( E[1] < E[2] ) { 
		tmp = E[1]; E[1] = E[2]; E[2] = tmp; 
	}

	// S(3,1) = the singular values of M
	S[0] = sqrt(E[0]);
	S[1] = sqrt(E[1]);
	S[2] = sqrt(E[2]);

	// rank
	rnk = 3;
	if ( S[0] < tol ) { rnk--; }
	if ( S[1] < tol ) { rnk--; }
	if ( S[2] < tol ) { rnk--; }

	// V = eigenvectors of the symmetric A
	switch (rnk) {
		case 0:
			V[0] = 1.; V[1] = 0.; V[2] = 0.; V[3] = 0.; V[4] = 1.; V[5] = 0.; V[6] = 0.; V[7] = 0.; V[8] = 1.;
			break;
		case 1:
			// the 1st eigenvector
			V[0] = A[3]*A[7]-A[6]*(A[4]-E[0]);
			V[1] = A[3]*A[6]-A[7]*(A[0]-E[0]);
			V[2] = (A[0]-E[0])*(A[4]-E[0])-A[3]*A[3];
			// check and normalize it
			tmp = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
			if ( tmp == 0. ) {
				V[0] = 1.; V[1] = 0.; V[2] = 0.;
			} else {
				V[0] /= tmp; V[1] /= tmp; V[2] /= tmp; 
			}
			// tmpVec = a unit axis vector
			idx = ( fabs(V[0])>fabs(V[1]) ? ( fabs(V[1])>fabs(V[2]) ? 2 : 1 ) : ( fabs(V[0])>fabs(V[2]) ? 2 : 0 ) );
			tmpVec[0] = 0.; tmpVec[1] = 0.; tmpVec[2] = 0.; tmpVec[idx] = 1.;
			// the direction of the 2nd eigenvector = Cross(1st eigenvector, tmpVec)
			V[3] = tmpVec[2]*V[1] - tmpVec[1]*V[2];
			V[4] = -(tmpVec[2]*V[0]) + tmpVec[0]*V[2];
			V[5] = tmpVec[1]*V[0] - tmpVec[0]*V[1];
			// normalize it
			tmp = sqrt(V[3]*V[3] + V[4]*V[4] + V[5]*V[5]);
			V[3] /= tmp; V[4] /= tmp; V[5] /= tmp; 
			// the 3rd eigenvector = cross product of the 1st and 2nd
			V[6] = V[1]*V[5] - V[2]*V[4];
			V[7] = V[2]*V[3] - V[0]*V[5];
			V[8] = V[0]*V[4] - V[1]*V[3];
			break;
		case 2:
			// the 1st eigenvector
			V[0] = A[3]*A[7]-A[6]*(A[4]-E[0]);
			V[1] = A[3]*A[6]-A[7]*(A[0]-E[0]);
			V[2] = (A[0]-E[0])*(A[4]-E[0])-A[3]*A[3];
			// check and normalize it
			tmp = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
			if ( tmp == 0. ) {
				V[0] = 1.; V[1] = 0.; V[2] = 0.;
			} else {
				V[0] /= tmp; V[1] /= tmp; V[2] /= tmp; 
			}
			// the 2nd eigenvector
			V[3] = A[3]*A[7]-A[6]*(A[4]-E[1]);
			V[4] = A[3]*A[6]-A[7]*(A[0]-E[1]);
			V[5] = (A[0]-E[1])*(A[4]-E[1])-A[3]*A[3];
			// check and normalize it
			tmp = sqrt(V[3]*V[3] + V[4]*V[4] + V[5]*V[5]);
			if ( tmp == 0. ) {
				// tmpVec = a unit axis vector
				idx = ( fabs(V[0])>fabs(V[1]) ? ( fabs(V[1])>fabs(V[2]) ? 2 : 1 ) : ( fabs(V[0])>fabs(V[2]) ? 2 : 0 ) );
				tmpVec[0] = 0.; tmpVec[1] = 0.; tmpVec[2] = 0.; tmpVec[idx] = 1.;
				// the direction of the 2nd eigenvector = Cross(1st eigenvector, tmpVec)
				V[3] = tmpVec[2]*V[1] - tmpVec[1]*V[2];
				V[4] = -(tmpVec[2]*V[0]) + tmpVec[0]*V[2];
				V[5] = tmpVec[1]*V[0] - tmpVec[0]*V[1];
				// normalize it
				tmp = sqrt(V[3]*V[3] + V[4]*V[4] + V[5]*V[5]);
				V[3] /= tmp; V[4] /= tmp; V[5] /= tmp; 
			} else {
				V[3] /= tmp; V[4] /= tmp; V[5] /= tmp; 
			}
			// the 3rd eigenvector = cross product of the 1st and 2nd
			V[6] = V[1]*V[5] - V[2]*V[4];
			V[7] = V[2]*V[3] - V[0]*V[5];
			V[8] = V[0]*V[4] - V[1]*V[3];
			break;
		case 3: // same as case 2:
			// the 1st eigenvector
			V[0] = A[3]*A[7]-A[6]*(A[4]-E[0]);
			V[1] = A[3]*A[6]-A[7]*(A[0]-E[0]);
			V[2] = (A[0]-E[0])*(A[4]-E[0])-A[3]*A[3];
			// check and normalize it
			tmp = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
			if ( tmp == 0. ) {
				V[0] = 1.; V[1] = 0.; V[2] = 0.;
			} else {
				V[0] /= tmp; V[1] /= tmp; V[2] /= tmp; 
			}
			// the 2nd eigenvector
			V[3] = A[3]*A[7]-A[6]*(A[4]-E[1]);
			V[4] = A[3]*A[6]-A[7]*(A[0]-E[1]);
			V[5] = (A[0]-E[1])*(A[4]-E[1])-A[3]*A[3];
			// check and normalize it
			tmp = sqrt(V[3]*V[3] + V[4]*V[4] + V[5]*V[5]);
			if ( tmp == 0. ) {
				// tmpVec = a unit axis vector
				idx = ( fabs(V[0])>fabs(V[1]) ? ( fabs(V[1])>fabs(V[2]) ? 2 : 1 ) : ( fabs(V[0])>fabs(V[2]) ? 2 : 0 ) );
				tmpVec[0] = 0.; tmpVec[1] = 0.; tmpVec[2] = 0.; tmpVec[idx] = 1.;
				// the direction of the 2nd eigenvector = Cross(1st eigenvector, tmpVec)
				V[3] = tmpVec[2]*V[1] - tmpVec[1]*V[2];
				V[4] = -(tmpVec[2]*V[0]) + tmpVec[0]*V[2];
				V[5] = tmpVec[1]*V[0] - tmpVec[0]*V[1];
				// normalize it
				tmp = sqrt(V[3]*V[3] + V[4]*V[4] + V[5]*V[5]);
				V[3] /= tmp; V[4] /= tmp; V[5] /= tmp; 
			} else {
				V[3] /= tmp; V[4] /= tmp; V[5] /= tmp; 
			}
			// the 3rd eigenvector = cross product of the 1st and 2nd
			V[6] = V[1]*V[5] - V[2]*V[4];
			V[7] = V[2]*V[3] - V[0]*V[5];
			V[8] = V[0]*V[4] - V[1]*V[3];
			break;
	}

	// U = M*V*pInv(Diag(S))
	if ( fabs(S[0]) < eps ) { pinvS[0] = 0.; } else { pinvS[0] = (gReal)1./S[0]; }
	if ( fabs(S[1]) < eps ) { pinvS[1] = 0.; } else { pinvS[1] = (gReal)1./S[1]; }
	if ( fabs(S[2]) < eps ) { pinvS[2] = 0.; } else { pinvS[2] = (gReal)1./S[2]; }
	switch (rnk) {
		case 0:
			U[0] = 1.; U[1] = 0.; U[2] = 0.; U[3] = 0.; U[4] = 1.; U[5] = 0.; U[6] = 0.; U[7] = 0.; U[8] = 1.;
			break;
		case 1:
			// 1st column of U = M*(1st column of V)*pinvS[0]
			U[0] = (M[0]*V[0] + M[3]*V[1] + M[6]*V[2]) * pinvS[0];
			U[1] = (M[1]*V[0] + M[4]*V[1] + M[7]*V[2]) * pinvS[0];
			U[2] = (M[2]*V[0] + M[5]*V[1] + M[8]*V[2]) * pinvS[0];
			// normalize it
			tmp = sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
			U[0] /= tmp; U[1] /= tmp; U[2] /= tmp; 
			// 2nd column of U = a unit vector which is perpendicular to the 1st column vector
			idx = ( fabs(U[0])>fabs(U[1]) ? ( fabs(U[1])>fabs(U[2]) ? 2 : 1 ) : ( fabs(U[0])>fabs(U[2]) ? 2 : 0 ) );
			tmpVec[0] = 0.; tmpVec[1] = 0.; tmpVec[2] = 0.; tmpVec[idx] = 1.;
			U[3] = tmpVec[2]*U[1] - tmpVec[1]*U[2];
			U[4] = -(tmpVec[2]*U[0]) + tmpVec[0]*U[2];
			U[5] = tmpVec[1]*U[0] - tmpVec[0]*U[1];
			// normalize it
			tmp = sqrt(U[3]*U[3] + U[4]*U[4] + U[5]*U[5]);
			U[3] /= tmp; U[4] /= tmp; U[5] /= tmp; 
			// 3rd column of U = cross product of the 1st and 2nd columns
			U[6] = U[1]*U[5] - U[2]*U[4];
			U[7] = U[2]*U[3] - U[0]*U[5];
			U[8] = U[0]*U[4] - U[1]*U[3];
			if ( detM < 0 ) { U[6] *= -1.; U[7] *= -1.; U[8] *= -1.; }
			break;
		case 2:
			// 1st column of U = M*(1st column of V)*pinvS[0]
			U[0] = (M[0]*V[0] + M[3]*V[1] + M[6]*V[2]) * pinvS[0];
			U[1] = (M[1]*V[0] + M[4]*V[1] + M[7]*V[2]) * pinvS[0];
			U[2] = (M[2]*V[0] + M[5]*V[1] + M[8]*V[2]) * pinvS[0];
			// normalize it
			tmp = sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
			U[0] /= tmp; U[1] /= tmp; U[2] /= tmp; 
			// 2nd column of U = M*(2nd column of V)*pinvS[1]
			U[3] = (M[0]*V[3] + M[3]*V[4] + M[6]*V[5]) * pinvS[1];
			U[4] = (M[1]*V[3] + M[4]*V[4] + M[7]*V[5]) * pinvS[1];
			U[5] = (M[2]*V[3] + M[5]*V[4] + M[8]*V[5]) * pinvS[1];
			// normalize it
			tmp = sqrt(U[3]*U[3] + U[4]*U[4] + U[5]*U[5]);
			U[3] /= tmp; U[4] /= tmp; U[5] /= tmp; 
			// 3rd column of U = cross product of the 1st and 2nd columns
			U[6] = U[1]*U[5] - U[2]*U[4];
			U[7] = U[2]*U[3] - U[0]*U[5];
			U[8] = U[0]*U[4] - U[1]*U[3];
			if ( detM < 0 ) { U[6] *= -1.; U[7] *= -1.; U[8] *= -1.; }
			break;
		case 3: // same as the case 2:
			// 1st column of U = M*(1st column of V)*pinvS[0]
			U[0] = (M[0]*V[0] + M[3]*V[1] + M[6]*V[2]) * pinvS[0];
			U[1] = (M[1]*V[0] + M[4]*V[1] + M[7]*V[2]) * pinvS[0];
			U[2] = (M[2]*V[0] + M[5]*V[1] + M[8]*V[2]) * pinvS[0];
			// normalize it
			tmp = sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
			U[0] /= tmp; U[1] /= tmp; U[2] /= tmp; 
			// 2nd column of U = M*(2nd column of V)*pinvS[1]
			U[3] = (M[0]*V[3] + M[3]*V[4] + M[6]*V[5]) * pinvS[1];
			U[4] = (M[1]*V[3] + M[4]*V[4] + M[7]*V[5]) * pinvS[1];
			U[5] = (M[2]*V[3] + M[5]*V[4] + M[8]*V[5]) * pinvS[1];
			// normalize it
			tmp = sqrt(U[3]*U[3] + U[4]*U[4] + U[5]*U[5]);
			U[3] /= tmp; U[4] /= tmp; U[5] /= tmp; 
			// 3rd column of U = cross product of the 1st and 2nd columns
			U[6] = U[1]*U[5] - U[2]*U[4];
			U[7] = U[2]*U[3] - U[0]*U[5];
			U[8] = U[0]*U[4] - U[1]*U[3];
			if ( detM < 0 ) { U[6] *= -1.; U[7] *= -1.; U[8] *= -1.; }
			break;
	}
}

bool MultAbCt(RMatrix &M, const RMatrix &A, const RMatrix &b, const RMatrix &C)
{
	int i, j, k;
	gReal tmp;

	if ( b.col != 1 ) return false;

	if ( A.col != b.row || C.col != b.row || b.col != 1 ) return false;

	M.ReNew(A.row, C.row);

	for (i=0; i<A.row; i++) {
		for (j=0; j<C.row; j++) {
			tmp = 0;
			for (k=0; k<b.row; k++) {
				tmp += A.element[i+A.row*k] * C.element[j+C.row*k] * b.element[k];
			}
			M.element[i+M.row*j] = tmp;
		}
	}

	return true;
}

bool MultAbCt(RMatrix &M, const RMatrix &A, const gReal *b, const RMatrix &C)
{
	int i, j, k;
	gReal tmp;

	if ( A.col != C.col ) return false;

	M.ReNew(A.row, C.row);

	for (i=0; i<A.row; i++) {
		for (j=0; j<C.row; j++) {
			tmp = 0;
			for (k=0; k<A.col; k++) {
				tmp += A.element[i+A.row*k] * C.element[j+C.row*k] * b[k];
			}
			M.element[i+M.row*j] = tmp;
		}
	}

	return true;
}

bool MultAb(RMatrix &M, const RMatrix &A, const RMatrix &b)
{
	int i, j;

	if ( A.col != b.row || b.col != 1) return false;

	M.ReNew(A.row, A.col);

	for (i=0; i<A.row; i++) {
		for (j=0; j<A.col; j++) {
			M.element[i+A.row*j] = A.element[i+A.row*j] * b.element[j];
		}
	}

	return true;
}

bool MultaB(RMatrix &M, const RMatrix &a, const RMatrix &B)
{
	int i, j;

	if ( a.row != B.row || a.col != 1) return false;

	M.ReNew(B.row, B.col);

	for (i=0; i<B.row; i++) {
		for (j=0; j<B.col; j++) {
			M.element[i+B.row*j] = B.element[i+B.row*j] * a.element[i];
		}
	}

	return true;
}

void _qsort(gReal *v, const int n, const int left, const int right, int *idx);
gReal Cond(const RMatrix &m)
{
	RMatrix s = SVD(m);
	_qsort(s.element, s.row, 0, s.row-1, NULL);
	return s.element[s.row-1] / s.element[0];
}

int Rank(const RMatrix &m, gReal eps)
{
	int i, rank = 0;
	RMatrix s = SVD(m);
	gReal imx = (gReal)1.0 / MaxVec(s);
	for ( i = 0; i < s.row; i++ ) if ( s.element[i] * imx > eps ) rank++;
	return rank;
}

int Rank(const RMatrix &m) { return Rank(m, (gReal)1E-9); }


//----------------------------------------------------------------------
//
//		title		:	LINPACK C Version
//						
//		version		:	v1.5
//		author		:	Jinwook Kim (zinook@plaza.snu.ac.kr)
//		last update	:	2001.7.23
//
//----------------------------------------------------------------------

void _dqrdc(gReal *x, int ldx, int n, int p, gReal *qraux, int *jpvt, gReal *work, int job)
/*	
	dqrdc uses householder transformations to compute the qr factorization of an n by p matrix x.  column pivoting
	based on the 2-norms of the reduced columns may be performed at the users option.

	on entry

	x		gReal precision(ldx,p), where ldx .ge. n. x contains the matrix whose decomposition is to be computed.
	ldx		int. ldx is the leading dimension of the array x.
	n		int. n is the number of rows of the matrix x.
	p		int. p is the number of columns of the matrix x.
	jpvt	int(p). jpvt contains ints that control the selection of the pivot columns.  the k-th column x(k) of x
			is placed in one of three classes according to the value of jpvt(k).
			if jpvt(k) .gt. 0, then x(k) is an initial column.
			if jpvt(k) .eq. 0, then x(k) is a free column.
			if jpvt(k) .lt. 0, then x(k) is a final column.
			before the decomposition is computed, initial columns are moved to the beginning of the array x and final
			columns to the end.  both initial and final columns are frozen in place during the computation and only
			free columns are moved.  at the k-th stage of the reduction, if x(k) is occupied by a free column
			it is interchanged with the free column of largest reduced norm.  jpvt is not referenced if
			job .eq. 0.
	work	gReal precision(p). work is a work array.  work is not referenced if job .eq. 0.
	job		int. job is an int that initiates column pivoting.
			if job .eq. 0, no pivoting is done.
			if job .ne. 0, pivoting is done.

	on return

	x		x contains in its upper triangle the upper triangular matrix r of the qr factorization.
			below its diagonal x contains information from which the orthogonal part of the decomposition
			can be recovered.  note that if pivoting has been requested, the decomposition is not that
			of the original matrix x but that of x with its columns permuted as described by jpvt.
	qraux	gReal precision(p). qraux contains further information required to recover the orthogonal part of the decomposition.
	jpvt	jpvt(k) contains the index of the column of the original matrix that has been interchanged into the k-th column, 
			if pivoting was requested.
*/
{
	int j, l, pl, pu, maxj, jp, lup, k;
	gReal maxnrm, tt, nrmxl, t, tmp;
	pl = 0;
	pu = -1;
	if ( job != 0 )
	{
		// pivoting has been requested.  rearrange the columns according to jpvt.
		for ( j = 0; j < p; j++ )
		{
			jpvt[j] = j+1;
			if ( jpvt[j] < 0 ) jpvt[j] = -j+1;
			if ( jpvt[j] > 0 )
			{
				if ( j != pl ) for ( k = 0; k < n; k++ ) { tmp = x[k+pl*ldx]; x[k+pl*ldx] = x[k+j*ldx]; x[k+j*ldx] = tmp; }
				jpvt[j] = jpvt[pl];
				jpvt[pl++] = j+1;
			} 
		}
		pu = p-1;
		for ( j = p-1; j >= 0; j-- )
		{
			if ( jpvt[j] < 0 )
			{
				jpvt[j] = -jpvt[j];
				if ( j != pu )
				{
					for ( k = 0; k < n; k++ ) { tmp = x[k+pu*ldx]; x[pu+pl*ldx] = x[k+j*ldx]; x[k+j*ldx] = tmp; }
					jp = jpvt[pu];
					jpvt[pu] = jpvt[j];
					jpvt[j] = jp;
				}
				pu--;
			}
		}
	}
	// compute the norms of the free columns.
	if ( pu >= pl )
	{
		for ( j = pl; j <= pu; j++ )
		{
			tmp = 0.0; 
			for ( k = 0; k < n; k++ ) tmp += x[k+j*ldx]*x[k+j*ldx];
			qraux[j] = sqrt(tmp);
			work[j] = qraux[j];
		}
	}
	// perform the householder reduction of x.
	lup = ( n > p ? p : n );
	for ( l = 0; l < lup; l++ )
	{
		if ( l >= pl && l < pu )
		{
			// locate the column of largest norm and bring it into the pivot position.
			maxnrm = 0.0;
			maxj = l;
            for ( j = l; j <= pu; j++ )
			{
				if ( qraux[j] > maxnrm )
				{
					maxnrm = qraux[j];
					maxj = j;
				}
			}
            if ( maxj != l )
			{
				for ( k = 0; k < n; k++ ) { tmp = x[k+l*ldx]; x[k+l*ldx] = x[k+maxj*ldx]; x[k+maxj*ldx] = tmp; }
				qraux[maxj] = qraux[l];
				work[maxj] = work[l];
				jp = jpvt[maxj];
				jpvt[maxj] = jpvt[l];
				jpvt[l] = jp;
			}
		}
		qraux[l] = 0.0;
		if ( l != n-1 )
		{
			// compute the householder transformation for column (l+1).
			tmp = 0.0; 
			for ( k = 0; k < n-l; k++ ) tmp += x[l+k+l*ldx]*x[l+k+l*ldx];
			nrmxl = sqrt(tmp);
            if ( nrmxl != 0.0 )
			{
				if ( x[l+l*ldx] != 0.0 ) nrmxl = dsign(nrmxl,x[l+l*ldx]);
				for ( k = 0; k < n-l; k++ ) x[l+k+l*ldx] *= (gReal)1.0/nrmxl;
				x[l+l*ldx] += 1.0;
				// apply the transformation to the remaining columns,
				// updating the norms.
				if ( p >= l+2 )
				{
					for ( j = l+1; j < p; j++ )
					{
						tmp = 0.0; 
						for ( k = 0; k < n-l; k++ ) tmp += x[l+k+l*ldx]*x[l+k+j*ldx];
						t = -tmp/x[l+l*ldx];
						for ( k = 0; k < n-l; k++ ) x[l+k+j*ldx] += t*x[l+k+l*ldx];
						if ( j >= pl && j <= pu && qraux[j] != 0.0 )
						{
							tt = (gReal)1.0 - pow(fabs(x[l+j*ldx])/qraux[j], 2);
							tt = ( tt > 0.0 ? tt : (gReal)0.0 );
							t = tt;
							tt = (gReal)1.0 + (gReal)0.05*tt*pow(qraux[j]/work[j],2);
							if ( tt == 1.0 )
							{
								tmp = (gReal)0.0;
								for ( k = 0; k < n-l-1; k++ ) tmp += x[l+1+k+j*ldx]*x[l+1+k+j*ldx];
								qraux[j] = sqrt(tmp);
								work[j] = qraux[j];
							} else qraux[j] *= sqrt(t);
						}
					}
				}
				// save the transformation.
				qraux[l] = x[l+l*ldx];
				x[l+l*ldx] = -nrmxl;
			}
		}
	}
	return;
}

void _dqrsl(gReal *x, int ldx, int n, int k, gReal *qraux, const gReal *y, gReal *qy, gReal *qty, gReal *b, gReal *rsd, gReal *xb, int job, int &info)
/*
	dqrsl applies the output of dqrdc to compute coordinate transformations, projections, and least squares solutions.
	for k .le. _MIN(n,p), let xk be the matrix 

		xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))

	formed from columnns jpvt(1), ... ,jpvt(k) of the original n x p matrix x that was input to dqrdc (if no pivoting was
	done, xk consists of the first k columns of x in their original order).  dqrdc produces a factored orthogonal matrix q
	and an upper triangular matrix r such that

		xk = q * (r)
				 (0)

	this information is contained in coded form in the arrays x and qraux.

	on entry

	x		gReal precision(ldx,p). x contains the output of dqrdc.
	ldx		integer. ldx is the leading dimension of the array x.
	n		integer. n is the number of rows of the matrix xk.  it must have the same value as n in dqrdc.
	k		integer. k is the number of columns of the matrix xk.  k must nnot be greater than _MIN(n,p), where p is the
			same as in the calling sequence to dqrdc.
	qraux	gReal precision(p). qraux contains the auxiliary output from dqrdc.
	y		gReal precision(n). y contains an n-vector that is to be manipulated by dqrsl.
	job		integer. job specifies what is to be computed.  job has the decimal expansion abcde, with the following meaning.
				if a.ne.0, compute qy.
				if b,c,d, or e .ne. 0, compute qty.
				if c.ne.0, compute b.
				if d.ne.0, compute rsd.
				if e.ne.0, compute xb.
			note that a request to compute b, rsd, or xb automatically triggers the computation of qty, for
			which an array must be provided in the calling sequence.

     on return

	qy		gReal precision(n). qy conntains q*y, if its computation has been requested.
	qty		gReal precision(n). qty contains trans(q)*y, if its computation has been requested.  here trans(q) is the
			transpose of the matrix q.
	b		gReal precision(k) b contains the solution of the least squares problem
			minimize norm2(y - xk*b),

	if its computation has been requested.  (note that if pivoting was requested in dqrdc, the j-th
	component of b will be associated with column jpvt(j) of the original matrix x that was input into dqrdc.)
	
	rsd    gReal precision(n). rsd contains the least squares residual y - xk*b, if its computation has been requested.  rsd is
	also the orthogonal projection of y onto the orthogonal complement of the column space of xk.
	
	xb     gReal precision(n). xb contains the least squares approximation xk*b, if its computation has been requested.  xb is also
	the orthogonal projection of y onto the column space of x.
	
	info   integer. info is zero unless the computation of b has been requested and r is exactly singular.  
	in this case, info is the index of the first zero diagonal element of r and b is left unaltered.

	the parameters qy, qty, b, rsd, and xb are not referenced if their computation is not requested and in this case
	can be replaced by dummy variables in the calling program. to save storage, the user may in some cases use the same
	array for different parameters in the calling sequence.  a frequently occuring example is when one wishes to compute
	any of b, rsd, or xb and does not need y or qty.  in this case one may identify y, qty, and one of b, rsd, or xb, while
	providing separate arrays for anything else that is to be computed.  thus the calling sequence

		call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)

	will result in the computation of b and rsd, with rsd overwriting y.  more generally, each item in the following
	list contains groups of permissible identifications for a single callinng sequence.
		1. (y,qty,b) (rsd) (xb) (qy)
		2. (y,qty,rsd) (b) (xb) (qy)
		3. (y,qty,xb) (b) (rsd) (qy)
		4. (y,qy) (qty,b) (rsd) (xb)
		5. (y,qy) (qty,rsd) (b) (xb)
		6. (y,qy) (qty,xb) (b) (rsd)
	in any group the value returned in the array allocated to the group corresponds to the last member of the group.
*/
{
	int i, j, ju, l;
	gReal t,temp;
	int cb, cqy, cqty, cr, cxb;
	info = 0;
	// determine what is to be computed.

	cqy = ( job / 10000 != 0 );
	cqty = ( ( job % 10000 ) != 0 );
	cb = ( ( job % 1000 ) / 100 != 0 );
	cr = ( ( job % 100 ) / 10 != 0 );
	cxb = ( ( job % 10 ) !=  0 );
	ju = ( k > n-1 ? n-1 : k );

	// special action when n=1.

	if ( ju == 0 )
	{
		if ( cqy ) qy[0] = y[0];
		if ( cqty ) qty[0] = y[0];
		if ( cxb ) xb[0] = y[0];
		if ( cb )
		{
			if ( x[0] != 0.0 ) 
			{
				b[0] = y[0] / x[0];
				for ( j = 1; j < k; j++) b[j] = 0.0;
			}
			else info = 1;
		}
		if ( cr ) rsd[0] = 0.0;
		return;
	}

	// set up to compute qy or qty.

	if ( cqy ) for ( l = 0; l < n; l++ ) qy[l] = y[l];
	if ( cqty ) for ( l = 0; l < n; l++ ) qty[l] = y[l];
	if ( cqy )
	{
		// compute qy.
		for ( j = ju-1; j >= 0; j-- )
		{
			if ( qraux[j] != 0.0 )
			{
				temp = x[j+j*ldx];
				x[j+j*ldx] = qraux[j];
				t = 0.0; 
				for ( l = 0; l < n-j; l++ ) t += x[j+l+j*ldx]*qy[j+l];
				t /= -x[j+j*ldx];
				for ( l = 0; l < n-j; l++ ) qy[j+l] += t*x[j+l+j*ldx];
				x[j+j*ldx] = temp;
			}
		}
	}
	if ( cqty )
	{
		// compute trans(q)*y
		for ( j = 0; j+1 <= ju; j++ )
		{
			if ( qraux[j] != 0.0 )
			{
				temp = x[j+j*ldx];
				x[j+j*ldx] = qraux[j];
				t = 0.0; 
				for ( l = 0; l < n-j; l++ ) t += x[j+l+j*ldx]*qty[j+l];
				t /= -x[j+j*ldx];
				for ( l = 0; l < n-j; l++ ) qty[j+l] += t*x[j+l+j*ldx];
				x[j+j*ldx] = temp;
			}
		}
	}
	// set up to compute b, rsd, or xb.
	if ( cb ) for ( l = 0; l < k; l++ ) b[l] = qty[l];
	if ( cxb ) for ( l = 0; l < k; l++ ) xb[l] = qty[l];
	if ( cr && k < n) for ( l = 0; l < n-k; l++ ) rsd[l] = qty[l];
	if ( cxb && k < n )
		for ( i = k; i < n; i++ ) xb[i] = 0.0;
	
	if ( cr )
		for ( i = 0; i < k; i++ ) rsd[i] = 0.0;

	if ( cb )
	{
		// compute b.
		for ( j = k-1; j >= 0; j-- )
		{
			if ( j < _MIN(k,n) )
			{
				if ( fabs(x[j+j*ldx]) < _EPS )
				{
					b[j] = 0.0;
					info = j+1;
					//break;
				} else b[j] /= x[j+j*ldx];
				
				if ( j != 0 )
				{
					t = -b[j];
					for ( l = 0; l < j; l++ ) b[l] += t*x[l+j*ldx];
				}
			} else b[j] = 0.0;
		}
	}
	if ( cr || cxb )
	{
		// compute rsd or xb as required.
		for ( j = ju-1; j >= 0; j-- )
		{
			if ( qraux[j] != 0.0 )
			{
				temp = x[j+j*ldx];
				x[j+j*ldx] = qraux[j];
				if ( cr )
				{
					t = 0.0; 
					for ( l = 0; l < n-j; l++ ) t += x[j+l+j*ldx]*rsd[j+l];
					t /= -x[j+j*ldx];
					for ( l = 0; l < n-j; l++ ) rsd[j+l] += t*x[j+l+j*ldx];
				}
				if ( cxb )
				{
					t = 0.0; 
					for ( l = 0; l < n-j; l++ ) t += x[j+l+j*ldx]*xb[j+l];
					t /= -x[j+j*ldx];
					for ( l = 0; l < n-j; l++ ) xb[j+l] += t*x[j+l+j*ldx];
				}
				x[j+j*ldx] = temp;
			}
		}
	}
	return;
}

// END OF LINPACK DQRDC

bool QRSolveAxEqualB(const RMatrix &A, RMatrix &x, const RMatrix &B)
{
	int i, j;
	static IMatrix _jpvt_QR;
	static RMatrix _qraux_QR, _qty_QR, _qy_QR, _work_QR, _A_QR;
	_A_QR.ReNew(A.row, A.col);
	for ( i = 0; i < A.row * A.col; i++ ) _A_QR.element[i] = A.element[i];

	int n = _A_QR.row, p = _A_QR.col, info;
	
	x.ReNew(p, B.col);

	if ( _qraux_QR.row < p ) _qraux_QR.ReNew(p);
	if ( _qty_QR.row < n ) _qty_QR.ReNew(n);
	if ( _work_QR.row < _MAX(p,n) ) _work_QR.ReNew(_MAX(p,n));
	if ( _jpvt_QR.row < p ) _jpvt_QR.ReNew(p);
	
	for ( i = 0; i < p; i++ ) _jpvt_QR.element[i] = 0;

	_dqrdc(_A_QR.element, n, n, p, _qraux_QR.element, _jpvt_QR.element, _work_QR.element, 1);
	
	for ( j = 0; j < B.col; j++ ) 
	{
		_dqrsl(_A_QR.element, n, n, p, _qraux_QR.element, B.element+j*B.row, NULL, _qty_QR.element, _work_QR.element, NULL, NULL, 100, info);
		for ( i = 0; i < p; i++ ) x.element[_jpvt_QR.element[i]-1+j*x.row] = _work_QR.element[i];
	}
	
	return true;
}

bool QRSolveAtxEqualB(const RMatrix &At, RMatrix &x, const RMatrix &B)
{
	int i, j;
	static IMatrix _jpvt_QR;
	static RMatrix _qraux_QR, _qty_QR, _qy_QR, _work_QR, _A_QR;

	_A_QR.ReNew(At.col, At.row);
	gReal *a = _A_QR.element, *at;
	for ( i = 0; i < _A_QR.col; i++ )
	{
		at = At.element + i;
		for ( j = 0; j < _A_QR.row; j++, at += At.row )	*(a++) = *at;
	}
	
	int n = _A_QR.row, p = _A_QR.col, info;
	x.ReNew(p, B.col);

	if ( _qraux_QR.row < p ) _qraux_QR.ReNew(p);
	if ( _qty_QR.row < n ) _qty_QR.ReNew(n);
	if ( _work_QR.row < _MAX(p,n) ) _work_QR.ReNew(_MAX(p,n));
	if ( _jpvt_QR.row < p ) _jpvt_QR.ReNew(p);
	
	for ( i = 0; i < p; i++ ) _jpvt_QR.element[i] = 0;

	_dqrdc(_A_QR.element, n, n, p, _qraux_QR.element, _jpvt_QR.element, _work_QR.element, 1);
	
	for ( j = 0; j < B.col; j++ ) 
	{
		_dqrsl(_A_QR.element, n, n, p, _qraux_QR.element, B.element+j*B.row, NULL, _qty_QR.element, _work_QR.element, NULL, NULL, 100, info);
		for ( i = 0; i < p; i++ ) x.element[_jpvt_QR.element[i]-1+j*x.row] = _work_QR.element[i];
	}
	
	return true;
}

RMatrix Nullspace(const RMatrix &A, gReal eps)
{
	int i, cnt;
	RMatrix N, U, s, V;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) return Eye(A.ColSize());
	gReal imx = (gReal)1.0 / mx;
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) cnt++; }
	N.SetZero(V.row, cnt);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) N.Push(0,cnt++,V.Sub(0,V.row-1,i,i)); }
	return N;	
}

RMatrix pInv(const RMatrix &A, gReal eps)
{
	RMatrix re, U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) return Zeros(A.ColSize(),A.RowSize());
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(A.col,1);
	for (int i=0; i<s.row; i++) { if ( s.element[i] * imx > eps ) is.element[i] = (gReal)1. / s.element[i]; }
	MultAbCt(re, V, is, U);	// re = V*Diag(is)*~U
	return re;	
}

RMatrix pInv(const RMatrix &A, RMatrix &N, gReal eps)
{
	int i, cnt;
	RMatrix re, U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { N = Eye(A.ColSize()); return Zeros(A.ColSize(),A.RowSize()); }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(A.col,1);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx > eps ) { is.element[i] = (gReal)1./s.element[i]; } else cnt++; }
	N.SetZero(V.row, cnt);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) N.Push(0,cnt++,V.Sub(0,V.row-1,i,i)); }
	MultAbCt(re, V, is, U);	// re = V*Diag(is)*~U
	return re;	
}

RMatrix srInv(const RMatrix &A, gReal alpha)
{
	if ( alpha <= 0.0 ) return pInv(A);

	RMatrix re, U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { return Zeros(A.ColSize(),A.RowSize()); }
	is.SetZero(A.col,1);
	for (int i=0; i<s.row; i++) { is.element[i] = s.element[i] / (s.element[i] * s.element[i] + alpha); }
	MultAbCt(re, V, is, U);	// re = V*Diag(is)*~U
	return re;	
}

RMatrix srInv(const RMatrix &A, RMatrix &N, gReal alpha, gReal eps)
{
	if ( alpha <= 0.0 ) return pInv(A, N, eps);

	int i, cnt;
	RMatrix re, U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { N = Eye(A.ColSize()); return Zeros(A.ColSize(),A.RowSize()); }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(A.col,1);
	for (i=0, cnt=0; i<s.row; i++) { is.element[i] = s.element[i]/(s.element[i]*s.element[i] + alpha); if ( s.element[i] * imx <= eps ) cnt++; }
	N.SetZero(V.row, cnt);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) N.Push(0,cnt++,V.Sub(0,V.row-1,i,i)); }
	MultAbCt(re, V, is, U);	// re = V*Diag(is)*~U
	return re;	
}

RMatrix srInv2(const RMatrix &A, gReal alpha)
{
	if ( A.row <= A.col ) {
		return ~A * Inv( A * ~A + alpha*Eye(A.row) );
	} else {
		return Inv( ~A * A + alpha*Eye(A.col) ) * ~A;
	} 
}

int solve_Ax_b_pInv(RMatrix &x, const RMatrix &A, const RMatrix &b, gReal eps)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return 0; }
	int rank=0;
	RMatrix U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { x = Zeros(A.ColSize(),1); return 0; }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(s.row,1);
	for (int i=0; i<s.row; i++) { if ( s.element[i] * imx > eps ) is.element[i] = 1/s.element[i]; rank++; }
	MultAbCt(x, V, is, ~b*U);	// x = V*Diag(is)*~U*b
	return rank;
}

int solve_Ax_b_pInv(RMatrix &x, const RMatrix &A, const RMatrix &b, RMatrix &N, gReal eps)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return 0; }
	int i, cnt, rank=0;
	RMatrix U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { N = Eye(A.ColSize()); x = Zeros(A.ColSize(),1); return 0; }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(s.row,1);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx > eps ) { is.element[i] = 1/s.element[i]; rank++; } else cnt++; }
	N.SetZero(V.row, cnt);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) N.Push(0,cnt++,V.Sub(0,V.row-1,i,i)); }
	MultAbCt(x, V, is, ~b*U);	// x = V*Diag(is)*~U*b
	return rank;
}

void solve_Ax_b_srInv(RMatrix &x, const RMatrix &A, const RMatrix &b, gReal alpha)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return; }
	if ( alpha <= 0.0 ) { solve_Ax_b_pInv(x, A, b); return; }
	RMatrix U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { x = Zeros(A.ColSize(),1); return; }
	is.SetZero(s.row,1);
	for (int i=0; i<s.row; i++) { is.element[i] = s.element[i]/(s.element[i]*s.element[i] + alpha); }
	MultAbCt(x, V, is, ~b*U);	// x = V*Diag(is)*~U*b
}

void solve_Ax_b_srInv(RMatrix &x, const RMatrix &A, const RMatrix &b, RMatrix &N, gReal alpha, gReal eps)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return; }
	if ( alpha <= 0.0 ) { solve_Ax_b_pInv(x, A, b, N, eps); return; }
	int i, cnt;
	RMatrix U, s, V, is;
	SVD(A, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { N = Eye(A.ColSize()); x = Zeros(A.ColSize(),1); return; }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(s.row,1);
	for (i=0, cnt=0; i<s.row; i++) { is.element[i] = s.element[i]/(s.element[i]*s.element[i] + alpha); if ( s.element[i] * imx <= eps ) cnt++; }
	N.SetZero(V.row, cnt);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) N.Push(0,cnt++,V.Sub(0,V.row-1,i,i)); }
	MultAbCt(x, V, is, ~b*U);	// x = V*Diag(is)*~U*b
}

int solve_Ax_b_pInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, gReal eps)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return 0; }
	int rank=0;
	RMatrix A_, U, s, V, is, x_;
	MultAb(A_, A, iw); // A_ = A * Diag(iw);
	SVD(A_, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { x = Zeros(A.ColSize(),1); return 0; }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(s.row,1);
	for (int i=0; i<s.row; i++) { if ( s.element[i] * imx > eps ) is.element[i] = 1/s.element[i]; rank++; }
	MultAbCt(x_, V, is, ~b*U);	// x_ = V*Diag(is)*~U*b
	MultaB(x, iw, x_); // x = Diag(iw) * x_;
	return rank;
}

int solve_Ax_b_pInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, RMatrix &Nw, gReal eps)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return 0; }
	int i, cnt, rank=0;
	RMatrix A_, U, s, V, is, N_, x_;
	MultAb(A_, A, iw); // A_ = A * Diag(iw);
	SVD(A_, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { Nw = Eye(A.ColSize()); x = Zeros(A.ColSize(),1); return 0; }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(s.row,1);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx > eps ) { is.element[i] = 1/s.element[i]; rank++; } else cnt++; }
	N_.SetZero(V.row, cnt);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) N_.Push(0,cnt++,V.Sub(0,V.row-1,i,i)); }
	MultAbCt(x_, V, is, ~b*U);	// x_ = V*Diag(is)*~U*b
	MultaB(x, iw, x_); // x = Diag(iw) * x_;
	MultaB(Nw, iw, N_); // Nw = Diag(iw) * N_;
	return rank;
}

void solve_Ax_b_srInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, gReal alpha)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return; }
	if ( alpha <= 0.0 ) { solve_Ax_b_pInvW(x, A, b, iw); return; }
	RMatrix A_, U, s, V, is, x_;
	MultAb(A_, A, iw); // A_ = A * Diag(iw);
	SVD(A_, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { x = Zeros(A.ColSize(),1); return; }
//	gReal imx = (gReal)1.0 / mx;
	is.SetZero(s.row,1);
	for (int i=0; i<s.row; i++) { is.element[i] = s.element[i]/(s.element[i]*s.element[i] + alpha); }
	MultAbCt(x_, V, is, ~b*U);	// x_ = V*Diag(is)*~U*b
	MultaB(x, iw, x_); // x = Diag(iw) * x_;
}

void solve_Ax_b_srInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, RMatrix &Nw, gReal alpha, gReal eps)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return; }
	if ( alpha <= 0.0 ) { solve_Ax_b_pInvW(x, A, b, iw, Nw, eps); return; }
	int i, cnt;
	RMatrix A_, U, s, V, is, N_, x_;
	MultAb(A_, A, iw); // A_ = A * Diag(iw);
	SVD(A_, U, s, V);
	gReal mx = MaxVec(s);
	if ( mx < 1E-20 ) { Nw = Eye(A.ColSize()); x = Zeros(A.ColSize(),1); return; }
	gReal imx = (gReal)1.0 / mx;
	is.SetZero(s.row,1);
	for (i=0, cnt=0; i<s.row; i++) { is.element[i] = s.element[i]/(s.element[i]*s.element[i] + alpha); if ( s.element[i] * imx <= eps ) cnt++; }
	N_.SetZero(V.row, cnt);
	for (i=0, cnt=0; i<s.row; i++) { if ( s.element[i] * imx <= eps ) N_.Push(0,cnt++,V.Sub(0,V.row-1,i,i)); }
	MultAbCt(x_, V, is, ~b*U);	// x_ = V*Diag(is)*~U*b
	MultaB(x, iw, x_); // x = Diag(iw) * x_;
	MultaB(Nw, iw, N_); // Nw = Diag(iw) * N_;
}

void solve_Ax_b_srInv2(RMatrix &x, const RMatrix &A, const RMatrix &b, gReal alpha, RMatrix &s)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return; }
	RMatrix U, V, is;
	SVD(A, U, s, V);
	is.SetZero(s.row,1);
	for (int i=0; i<s.row; i++) { 
		is.element[i] = (s.element[i]*s.element[i])/(s.element[i]*s.element[i]*s.element[i] + alpha); 
	}
	MultAbCt(x, V, is, ~b*U);	// x = V*Diag(is)*~U*b
}

void solve_Ax_b_srInvW2(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, gReal alpha, RMatrix &s)
{
	if ( A.row != b.row || A.col == 0 ) { x = Zeros(A.row,1); return; }
	RMatrix A_, U, V, is, x_;
	MultAb(A_, A, iw); // A_ = A * Diag(iw);
	SVD(A_, U, s, V);
	is.SetZero(s.row,1);
	for (int i=0; i<s.row; i++) { is.element[i] = (s.element[i]*s.element[i])/(s.element[i]*s.element[i]*s.element[i] + alpha); }
	MultAbCt(x_, V, is, ~b*U);	// x_ = V*Diag(is)*~U*b
	MultaB(x, iw, x_); // x = Diag(iw) * x_;
}

