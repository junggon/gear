//====================================================================================================
//
//      Title :   rmatrix3j.h
//
//      Note  :   This is a modified version of the original rmatrix3.h by Jinwook Kim (see below).
//                Modified and maintained by Junggon Kim (junggon@gmail.com)
//
//====================================================================================================

//////////////////////////////////////////////////////////////////////////////////
//
//		title		:	rmatrix3.h
//						
//		version		:	v2.89
//		author		:	Jinwook Kim (zinook@plaza.snu.ac.kr)
//		last update	:	2001.7.23
//
//////////////////////////////////////////////////////////////////////////////////

/*
	RMatrix3 is designed using template class.
	However it is highly recommended to use only type of 'double'.
	In case of using type of 'int' or 'float', there may be some internal conflictions with type of 'double'.
	Hence if you find any errors or incompatibility in using _rmatrix <int>, please do not inform me that.

	If you want not to check index for high performance, use directive '#define RMATRIX3_DO_NOT_CHECK_INDEX' before including 'rmatrix3.h'
	If you want not to check size for high performance, use directive '#define RMATRIX3_DO_NOT_CHECK_SIZE' before including 'rmatrix3.h'
	Usually, it is recommended to use the directives after complete debugging.

	template <class TYPE> class _rmatrix
	method : 
	// constructors
		_rmatrix()
		_rmatrix(int r, int c)
		_rmatrix(const _rmatrix &m)
		_rmatrix(int r, int c, const Type d[])
		_rmatrix(const char * fmt)
	// destructor
		~_rmatrix()
	
	// ith element in column order - zero based index
		TYPE &operator [] (int i)
		const TYPE &operator [] (int i) const
	// ith row, jth column element - index depends on idxbase
		TYPE &operator () (int i, int j)
		const TYPE &operator () (int i, int j) const
	// unary plus operator
		const _rmatrix &operator + (void) const
	// unary minus operator
		_rmatrix operator - (void) const
	// transpose operator
		_rmatrix operator ~ (void) const
	// substitute operator
		const _rmatrix &operator = (const _rmatrix &m)
	// += operator
		const _rmatrix &operator += (const _rmatrix &m)
	// -= operator
		const _rmatrix &operator -= (const _rmatrix &m)
	// *= operator
		const _rmatrix &operator *= (TYPE c)
	// /= operator
		const _rmatrix &operator /= (TYPE c)
		_rmatrix operator + (const _rmatrix &m) const
		_rmatrix operator - (const _rmatrix &m) const
		_rmatrix operator * (const _rmatrix &m) const	:  (*this) *  m
		_rmatrix operator ^ (const _rmatrix &m) const	: ~(*this) *  m
		_rmatrix operator | (const _rmatrix &m) const	:  (*this) * ~m
		_rmatrix operator * (TYPE c) const
		_rmatrix operator / (TYPE c) const
		_rmatrix operator % (const _rmatrix &m) const	: Inv(*this) * m
		_rmatrix operator & (const _rmatrix &m) const	: Inv(~*this) * m

		int RowSize(void) const
		int ColSize(void) const
		const _rmatrix &ReSize(int r, int c)
		_rmatrix Sub(int rs, int re, int cs, int ce) const
		const _rmatrix &Push(int i, int j, const _rmatrix &m);
		_rmatrix &SetZero(void)
		_rmatrix &SetZero(int r, c)
		
		TYPE* GetPtr()
		const TYPE* GetPtr()

		void Format(const char * fmt)

		bool IsZero();

		friend _rmatrix operator + (TYPE c, _rmatrix m) 
		friend _rmatrix operator - (TYPE c, _rmatrix m) 
		friend _rmatrix operator * (TYPE c, _rmatrix m)
		friend _rmatrix operator / (TYPE c, _rmatrix m) 
		
		friend std::ostream &operator << (std::ostream &os, const _rmatrix &m)
		friend RMatrix Zeros(int r, int c);
		friend RMatrix Ones(int r, int c);
		friend RMatrix Rand(int r, int c);
		friend RMatrix Eye(int r, int c );
		friend RMatrix Eye(int r);
		friend Type SquareSum(_rmatrix m)
		friend Type SquareSumW(_rmatrix m, _rmatrix w);
		friend Type FNorm(_rmatrix m)
		friend void AMultB(_rmatrix &re, const _rmatrix &a, const _rmatrix &b)
		friend int GaussElimination(Rmatrix &A, IMatrix ipvt, IMatrix jpvt)
		friend bool SolveAxEqualB(const RMatrix &A, RMatrix &x, const RMmatrix &B)
		friend bool SolveAtxEqualB(const RMatrix &A, RMatrix &x, const RMmatrix &B)
		friend bool SolvePosDefAxEqualB(RMatrix &A, RMatrix &x, const RMatrix &B)
		friend void ConvColumn(RMatrix &re, const RMatrix &u, int i, const RMatrix &v, int j);
		friend RMatrix Conv(const RMatrix &u, const RMatrix &v);
		friend void Conv(RMatrix &re, const RMatrix &u, const RMatrix &v);
		friend gReal Det(RMatrix A);
		friend RMatrix Eig(RMatrix m);
		friend void Eig(RMatrix &re, RMatrix &m);
		friend void Eig(RMatrix m, RMatrix &v, RMatrix &d);
		friend RMatrix Companion(const RMatrix &m);
		friend void Companion(RMatrix &re, const RMatrix &m);
		friend RMatrix Roots(const RMatrix& cof);
		friend void Roots(RMatrix &re, const RMatrix& cof);
*/

#ifndef _RMatrix3j_
#define _RMatrix3j_

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "greal.h"

#define RMATRIX3_DO_NOT_CHECK_SIZE
#define RMATRIX3_DO_NOT_CHECK_INDEX

#define _EPS 2.2204E-16
#define _TINY 1E-20
#define _INF 1E100

#define _MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define _MIN(a, b)  (((a) < (b)) ? (a) : (b))

#define _ISSAMESIZE(m) ( row == m.row && col == m.col )
#define _ISEXIST0(i,j) ( i >= 0 && i <= row-1 && j >= 0 && j <= col-1 )
#define _ISEMPTY(m) ( m.row == 0 || m.col == 0 )

enum INDEXBASE { ZERO = 0, ONE = 1 };

template <class TYPE> class _rmatrix
{
public:
	int row, col;
	TYPE *element;
public:
	static INDEXBASE idxbase;

public:
	// default constructor
	_rmatrix() : row(0), col(0), element(NULL) { }
	
	// constructor with dimension
	_rmatrix(int r, int c) : row(r), col(c), element(new TYPE[r*c]) { }
	
	// copy constructor
	_rmatrix(const _rmatrix &m) : row(m.row), col(m.col), element(new TYPE[row*col]) 
	{ 
		TYPE *_t = element, *_m = m.element; 
		int n = row * col;	
		while ( n-- ) 
			*(_t++) = *(_m++); 
	}
	
	// constructor from single pointer, element is arranged along column
	_rmatrix(int r, int c, const TYPE d[]) : row(r), col(c), element(new TYPE[r*c]) { TYPE *_t = element; int n = row * col;	while ( n-- ) *(_t++) = *(d++); }

	// constructor from const char*
	_rmatrix(const char * fmt)
	{
		int i, j, m, n, count;
		int ch, idx, ddx, eol;
		static char buf[80]; 

		m = n = 0;
		idx = ddx = eol = i = j = 0;
		while( ! eol ) {
			ch = fmt[idx++];
			if( ch == '\t' || ch == '\r' || ch == '\n' || ch == ' ' || ch == ',' || ch == ';' || ch == '\0' ) {
				if( ddx > 0 ) { j++; ddx = 0; }
				if( ch == ';' || ch == '\0' ) { i++; if(n < j) n = j; j = 0; }
				if( ch == '\0' ) { m = i; eol = 1; }
			} else ddx++;
		}

		count = m * n;
		if( count > 0 ) {
			element = new TYPE [ count ];
			for(i = 0; i < count; i++) element[i] = (TYPE)0.0;
			row = m; col = n;
		} else {
			element = NULL;
			row = col = 0;
		}

		idx = ddx = eol = i = j = 0;
		while( ! eol ) {
			ch = fmt[idx++];
			if( ch == '\t' || ch == '\r' || ch == '\n' || ch == ' ' || ch == ',' || ch == ';' || ch == '\0' ) {
				if( ddx > 0 ) { buf[ddx] = '\0'; element[i + j*m] = (TYPE)(atof(buf)); j++; ddx = 0; }
				if( ch == ';' || ch == '\0' ) { i++; j = 0; }
				if( ch == '\0' ) eol = 1;
			} else buf[ddx++] = ch;
		}
	}

	// destructor
	~_rmatrix() { 
		delete [] element; 
	}

	////////////////////////////////////////////////////////////////
	//
	// operators
	//
	////////////////////////////////////////////////////////////////

	// ith element in column order : zero-base, i = 0 : row*col-1
#ifdef RMATRIX3_DO_NOT_CHECK_INDEX
	TYPE &operator [] (int i) { return element[i]; }
	const TYPE &operator [] (int i) const { return element[i]; }
#else
	TYPE &operator [] (int i) 
	{ 
		if ( i < 0 || i >= row*col ) 
			std::cerr << "RMatrix [int] : index over range" << std::endl; 
		return element[i]; 
	}
	const TYPE &operator [] (int i) const { if ( i < 0 || i >= row*col ) std::cerr << "RMatrix [int] : index out of range" << std::endl; return element[i]; }
#endif

	// (i, j)th element : index depends on 'idxbase'
#ifdef RMATRIX3_DO_NOT_CHECK_INDEX
	TYPE &operator () (int i, int j ) 
	{ i -= idxbase; j -= idxbase; return element[i+j*row]; }
	const TYPE &operator () (int i, int j ) const
	{ i -= idxbase; j -= idxbase; return element[i+j*row]; }
#else
	TYPE &operator () (int i, int j ) 
	{ i -= idxbase; j -= idxbase; if ( i < 0 || i >= row || j < 0 || j >= col ) std::cerr << "RMatrix (int, int) : index out of range" << std::endl; return element[i+j*row]; }
	const TYPE &operator () (int i, int j ) const
	{ i -= idxbase; j -= idxbase; if ( i < 0 || i >= row || j < 0 || j >= col ) std::cerr << "RMatrix (int, int) : index out of range" << std::endl; return element[i+j*row]; }
#endif
	
	// unary plus operator
	const _rmatrix &operator + (void) const { return *this; }
	
	// unary minus operator
	_rmatrix operator - (void) const
	{
		_rmatrix re(row, col);
		TYPE *_m = re.element, *_t = element;
		int n = row * col;
		while ( n-- ) *(_m++) = -*(_t++);
		return re; 
	}
	
	// transpose operator
	_rmatrix operator ~ (void) const 
	{ 
		_rmatrix re(col, row);
		int i = 0, r = row, c;
		TYPE *_m = re.element, *_mt;
		while ( r-- )
		{
			_mt = element + (i++);
			c = col;
			while ( c-- )
			{
				*(_m++) = *_mt;
				_mt += row;
			}			
		}	
		return re;		
	}
	
	// substitute operator
	const _rmatrix &operator = (const _rmatrix &m)
	{
		int n = m.row * m.col;
		if ( row * col != n ) 
		{
			delete [] element;
			element = new TYPE [ n ];
		}
		row = m.row;	col = m.col;
		TYPE *_t = element, *_m = m.element;
		while ( n-- ) *(_t++) = *(_m++);		
		return *this;
	}
	
	// += operator
	const _rmatrix &operator += (const _rmatrix &m)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( !_ISSAMESIZE(m) ) std::cerr << "RMatrix += : size not compatible" << std::endl;
#endif	
		int n = row * col;
		TYPE *_t = element, *_m = m.element;
		while ( n-- ) *(_t++) += *(_m++);
		return *this;
	}
	
	// -= operator
	const _rmatrix &operator -= (const _rmatrix &m)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( !_ISSAMESIZE(m) ) std::cerr << "RMatrix -= : size not compatible" << std::endl;
#endif		
		int n = row * col;
		TYPE *_t = element, *_m = m.element;
		while ( n-- ) *(_t++) -= *(_m++);
		return *this;
	}
	
	// *= operator
	const _rmatrix &operator *= (TYPE c)
	{
		int n = row * col;
		TYPE *_t = element;
		while ( n-- ) *(_t++) *= c;		
		return *this;
	}

	// -= operator
	const _rmatrix &operator /= (TYPE c)
	{
		int n = row * col;
		TYPE *_t = element, ci = (TYPE)1.0 / c;
		while ( n-- ) *(_t++) *= ci;		
		return *this;
	}

	// + operator 
	_rmatrix operator + (const _rmatrix &m) const
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( !_ISSAMESIZE(m) ) std::cerr << "RMatrix + : size not compatible" << std::endl;
#endif
		_rmatrix re(row, col);
		int n = row * col;
		TYPE *_t = element, *_r = re.element, *_m = m.element;
		while ( n-- ) *(_r++) = *(_t++) + *(_m++);
		return re;
	}

	// - operator 
	_rmatrix operator - (const _rmatrix &m) const
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( !_ISSAMESIZE(m) ) std::cerr << "RMatrix - : size not compatible" << std::endl;
#endif
		_rmatrix re(row, col);
		int n = row * col;
		TYPE *_t = element, *_r = re.element, *_m = m.element;
		while ( n-- ) *(_r++) = *(_t++) - *(_m++);
		return re;
	}

	// multiplication operator 
	_rmatrix operator * (const _rmatrix &m) const
	{
		_rmatrix re;
		AMultB(re, *this, m);
		return re;
	}

	// multiplication operator A * ~B
	_rmatrix operator | (const _rmatrix &m) const
	{
		_rmatrix re;
		AMultBt(re, *this, m);
		return re;	
	}

	// multiplication operator ~A * B
	_rmatrix operator ^ (const _rmatrix &m) const
	{
		_rmatrix re;
		AtMultB(re, *this, m);
		return re;
	}

	_rmatrix operator * (TYPE c) const
	{
		_rmatrix re(row, col);
		int n = row * col;
		TYPE *_t = element, *_r = re.element;
		while ( n-- ) *(_r++) = c * *(_t++);
		return re;
	}

	_rmatrix operator / (TYPE c) const
	{
		_rmatrix re(row, col);
		int n = row * col;
		TYPE *_t = element, *_r = re.element, ci = (TYPE)1.0 / c;
		while ( n-- ) *(_r++) = ci * *(_t++);
		return re;
	}

	// return Inv(*this) * m
	_rmatrix operator % (const _rmatrix &m) const
	{
		_rmatrix x;
		SolveAxEqualB(*this, x, m);
		return x;
	}

	// return Inv(~(*this)) * m
	_rmatrix operator & (const _rmatrix &m) const
	{
		_rmatrix x;
		SolveAtxEqualB(*this, x, m);
		return x;
	}

	////////////////////////////////////////////////////////////////
	//
	// member functions
	//
	////////////////////////////////////////////////////////////////

	// return number of row
	int RowSize(void) const { return row; }
	
	// return number of column
	int ColSize(void) const { return col; }
	
	// return row * col
	int Size(void) const { return row*col; }
	
	// resize matrix
	// elements of folding location are conserved.
	// Newly made elements are set zero.
	// In case of only modifying the size and not concerning the resulting values,
	// 'ReNew' is sufficient and faster than 'ReSize'.
	const _rmatrix &ReSize(int r, int c = 1)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_INDEX
		if ( r < 0 || c < 0 ) std::cerr << "RMatrix ReSize : invalid input" << std::endl;
#endif		
		if ( r != row || c != col )
		{
			_rmatrix tmp(r, c);
			int i, j;
			for ( i = 0; i < r; i++ )
			for ( j = 0; j < c; j++ )
			{
				if ( i < row && j < col ) tmp.element[i+j*r] = element[i+j*row];
				else tmp.element[i+j*r] = (TYPE)0.0;
			}
			TYPE *dum = element;
			element = tmp.element;
			tmp.element = dum;
			row = r;	col = c;
		}
		return *this;
	}

	// renew matrix
	// newly made elements are not initialized.
	const _rmatrix &ReNew(int r, int c = 1)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_INDEX
		if ( r < 0 || c < 0 ) std::cerr << "RMatrix ReNew : invalid input" << std::endl;
#endif		
		if ( r * c != row * col )
		{
			delete [] element;
			element = new TYPE [r * c];
		}
		row = r;	col = c;
		return *this;
	}

	// make sub matrix. index depends on 'idxbase'
#ifdef RMATRIX3_DO_NOT_CHECK_INDEX
	_rmatrix Sub(int rs, int re, int cs, int ce) const
	{
		int r = re-rs+1, c = ce-cs+1, i, j;
		_rmatrix Re(r, c);
		for ( i = 0; i < r; i++) 
		for ( j = 0; j < c; j++) Re.element[i+j*r] = element[(i+rs-idxbase)+(j+cs-idxbase)*row];
		return Re;
	}
#else
	_rmatrix Sub(int rs, int re, int cs, int ce) const
	{
		if ( rs > re || cs > ce || !_ISEXIST0(rs-idxbase, cs-idxbase) || !_ISEXIST0(re-idxbase, ce-idxbase) ) std::cerr << "RMatrix Sub : invalid input" << std::endl;
		int r = re-rs+1, c = ce-cs+1, i, j;
		_rmatrix Re(r, c);
		for ( i = 0; i < r; i++) 
		for ( j = 0; j < c; j++) Re.element[i+j*r] = element[(i+rs-idxbase)+(j+cs-idxbase)*row];
		return Re;
	}
#endif

	// push 'm' into the matrix
#ifdef RMATRIX3_DO_NOT_CHECK_INDEX
	const _rmatrix &Push(int r, int c, const _rmatrix &m)
	{
		int nr = _MAX(r-idxbase+m.row,row), nc = _MAX(c-idxbase+m.col, col), i, j;
			
		if ( r-idxbase+m.row > row || c-idxbase+m.col > col ) this->ReSize(nr, nc);		
		
		for ( i = r-idxbase; i < r-idxbase+m.row; i++)
		for ( j = c-idxbase; j < c-idxbase+m.col; j++) element[i+j*nr] = m.element[i-r+idxbase+(j-c+idxbase)*m.row];
		return *this;
	}
#else
	// push 'm' into the matrix
	const _rmatrix &Push(int r, int c, const _rmatrix &m)
	{
		int nr = _MAX(r-idxbase+m.row,row), nc = _MAX(c-idxbase+m.col, col), i, j;
			
		if ( _ISEMPTY(m) ) std::cerr << "RMatrix Push : empty argument" << std::endl;
		else if ( r < idxbase || c < idxbase ) std::cerr << "RMatrix Push : invalid input" << std::endl;
		else if ( r-idxbase+m.row > row || c-idxbase+m.col > col ) this->ReSize(nr, nc);		
		
		for ( i = r-idxbase; i < r-idxbase+m.row; i++)
		for ( j = c-idxbase; j < c-idxbase+m.col; j++) element[i+j*nr] = m.element[i-r+idxbase+(j-c+idxbase)*m.row];
		return *this;
	}
#endif


	_rmatrix &SetZero(void)
	{
		int n = row * col;
		TYPE *_t = element;
		while ( n-- ) *(_t++) = (TYPE)0.0;
		return *this;
	}

	_rmatrix &SetZero(int r, int c)
	{
		int n = r * c;
		if ( n != row * col ) { delete [] element; element = new TYPE [n]; }
		row = r; col = c;
		TYPE *_t = element;
		while ( n-- ) *(_t++) = (TYPE)0.0;
		return *this;
	}

	TYPE Normalize(void)
	{
		TYPE norm = FNorm(*this), inorm = 1.0 / norm;
		*this *= inorm;
		return norm;
	}

	TYPE* GetPtr() { return element; }
	const TYPE* GetPtr() const { return element; }

	void Format(const char * fmt) 
	{
		delete [] element;

		int i, j, m, n, count;
		int ch, idx, ddx, eol;
		static char buf[80]; 

		m = n = 0;
		idx = ddx = eol = i = j = 0;
		while( ! eol ) {
			ch = fmt[idx++];
			if( ch == '\t' || ch == '\r' || ch == '\n' || ch == ' ' || ch == ',' || ch == ';' || ch == '\0' ) {
				if( ddx > 0 ) { j++; ddx = 0; }
				if( ch == ';' || ch == '\0' ) { i++; if(n < j) n = j; j = 0; }
				if( ch == '\0' ) { m = i; eol = 1; }
			} else ddx++;
		}

		count = m * n;
		if( count > 0 ) {
			element = new TYPE [ count ];
			for(i = 0; i < count; i++) element[i] = (TYPE)0.0;
			row = m; col = n;
		} else {
			element = NULL;
			row = col = 0;
		}

		idx = ddx = eol = i = j = 0;
		while( ! eol ) {
			ch = fmt[idx++];
			if( ch == '\t' || ch == '\r' || ch == '\n' || ch == ' ' || ch == ',' || ch == ';' || ch == '\0' ) {
				if( ddx > 0 ) { buf[ddx] = '\0'; element[i + j*m] = (TYPE)(atof(buf)); j++; ddx = 0; }
				if( ch == ';' || ch == '\0' ) { i++; j = 0; }
				if( ch == '\0' ) eol = 1;
			} else buf[ddx++] = ch;
		}
	}

	bool IsZero()
	{
		int n = row * col;
		TYPE *_t = element;
		while ( n-- ) { if ( *(_t++) != (TYPE)0.0 ) return false; }
		return true;
	}

	////////////////////////////////////////////////////////////////
	//
	// friend functions
	//
	////////////////////////////////////////////////////////////////

	// concatenation of RMatrix
//	friend _rmatrix Concat(int r, int c, ...);

	// * operator : * c compoenet wise
	friend _rmatrix operator * (TYPE c, _rmatrix m)
	{ 
		int n = m.row * m.col;
		TYPE *_m = m.element;
		while ( n-- ) *(_m++) *= c;
		return m;		
	}
	
	// std::ostream standard output
	friend std::ostream &operator << (std::ostream &os, const _rmatrix &m)
	{
		os.setf(std::ios::scientific);
		//os.precision(16);
		os.precision(8);
		//os << "[" << std::endl;
		os << "[ % size = (" << m.row << ", " << m.col << ")" << std::endl;
		for ( int i = 0; i < m.row; i++ )
		{
			os << "  ";
			for ( int j = 0; j < m.col; j++ )
			{
				if ( m.element[i+j*m.row] >= (TYPE)0.0 ) os << " ";
				os << m.element[i+j*m.row];
				if ( j == m.col-1 ) os << " ;" << std::endl;
				else os << "  ";
			}
		}
		os << "];" << std::endl;
		os.unsetf(std::ios::scientific);
		return os;
	}
	
	friend TYPE SquareSum(const _rmatrix &m)
	{
		TYPE sum = (TYPE)0.0, *_m = m.element;
		int n = m.row * m.col;
		while ( n-- ) { sum += *_m * *_m; ++_m; }
		return sum;
	}

	friend TYPE SquareSumW(const _rmatrix &m, const _rmatrix &w)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( m.row != w.row ) std::cerr << "RMatrix SquareSumW : size mismatch" << std::endl;
		if ( m.col != w.col ) std::cerr << "RMatrix SquareSumW : size mismatch" << std::endl;
#endif
		TYPE sum = (TYPE)0.0, *_m = m.element, *_w = w.element;
		int n = m.row * m.col;
		while ( n-- ) sum += *_m * *(_m++) * *(_w++);
		return sum;
	}

//	friend TYPE Inner(const _rmatrix &a, const _rmatrix &b)
//	{
//#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
//		if ( _length(a) != _length(b) ) std::cerr << "size is not compatible" << std::endl;
//#endif
//		TYPE sum = (TYPE)0.0, *_tmpa = a.element, *_tmpb = b.element;
//		int n = _length(a);
//		while ( n-- ) sum += *(_tmpa++) * *(_tmpb++);
//		return sum;
//	}

	friend TYPE FNorm(const _rmatrix &m)
	{
		return sqrt(SquareSum(m));
	}
	
	friend void AMultB(_rmatrix &re, const _rmatrix &a, const _rmatrix &b)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( a.col != b.row ) std::cerr << "RMatrix AMultB : size not compatible" << std::endl;
#endif
		re.ReNew(a.row, b.col);
		int i, bc = b.col, k, ar;
		TYPE sum, *tmpa, *tmpb = b.element, *rij = re.element;
		while ( bc-- )
		{
			ar = a.row;
			i = 0;
			while ( ar-- )
			{
				tmpa = a.element + (i++);
				sum = (TYPE)0.0;
				k = a.col;
				while ( k-- )
				{
					sum += *tmpa * *tmpb;
					tmpa += a.row;
					tmpb++;
				}
				tmpb -= b.row;
				*(rij++) = sum;				
			}
			tmpb += b.row;
		}
	}

	friend void AMultBt(_rmatrix &re, const _rmatrix &a, const _rmatrix &b)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( a.col != b.col ) std::cerr << "RMatrix AMultBt : size not compatible" << std::endl;
#endif
		int i, j = 0, br = b.row, ar, ac;
		re.ReNew(a.row, br);
		TYPE sum, *tmpa, *tmpb, *rij = re.element;
		
		while ( br-- )
		{
			ar = a.row;
			i = 0;
			while ( ar-- )
			{
				tmpa = a.element + (i++);
				tmpb = b.element + j;
				sum = (TYPE)0.0;
				ac = a.col;
				while ( ac-- )
				{
					sum += *tmpa * *tmpb;
					tmpa += a.row;
					tmpb += b.row;
				}
				*(rij++) = sum;
			}
			j++;
		}
	}

	friend void AtMultB(_rmatrix &re, const _rmatrix &a, const _rmatrix &b)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( a.row != b.row ) std::cerr << "RMatrix AtMultB : size not compatible" << std::endl;
#endif
		re.ReNew(a.col, b.col);
		int ac, bc = b.col, ar;
		TYPE sum, *tmpa, *tmpb = b.element, *rij = re.element;
		while ( bc-- )
		{
			tmpa = a.element;
			ac = a.col;
			while ( ac-- )
			{
				sum = (TYPE)0.0;
				ar = a.row;
				while ( ar-- ) sum += *(tmpa++) * *(tmpb++);
					
				*(rij++) = sum;
				tmpb -= b.row;
			}
			tmpb += b.row;
		}
	}

//	friend TYPE Quadratic(const _rmatrix &x, const _rmatrix &A, const _rmatrix &y)
//	{
//#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
//		if ( _length(x) != A.row && _length(y) != A.col ) std::cerr << "size is not compatible" << std::endl;
//#endif
//		int r, c = A.col;
//		TYPE sum = (TYPE)0.0, xa, *tmpa = A.element, *tmpx, *tmpy = y.element;
//		while ( c-- )
//		{
//			xa = 0.0;
//			tmpx = x.element;
//			r = A.row;
//			while ( r-- ) xa += *(tmpx++) * *(tmpa++);			
//			sum += xa * *(tmpy++);
//		}
//		return sum;		
//	}

	friend TYPE MaxVec(const _rmatrix &m, int *idx = NULL)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( m.row * m.col <= 0 ) std::cerr << "RMatrix MaxVec : zero size" << std::endl;
#endif
		int n = m.row * m.col;
		TYPE mx = m.element[0];
		if ( idx != NULL ) *idx = 0;
		for ( int i = 1; i < n; i++ )
		{
			if ( m.element[i] > mx )
			{
				mx = m.element[i];
				if ( idx != NULL ) *idx = i;
			}
		}
		if ( idx != NULL ) *idx += idxbase;
		return mx;
	}

	friend TYPE MinVec(const _rmatrix &m, int *idx = NULL)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( m.row * m.col <= 0 ) std::cerr << "RMatrix MinVec : zero size" << std::endl;
#endif
		int n = m.row * m.col;
		TYPE mn = m.element[0];
		if ( idx != NULL ) *idx = 0;
		for ( int i = 1; i < n; i++ )
		{
			if ( m.element[i] < mn )
			{
				mn = m.element[i];
				if ( idx != NULL ) *idx = i;
			}
		}
		if ( idx != NULL ) *idx += idxbase;
		return mn;
	}

	friend TYPE Trace(_rmatrix &m)
	{
#ifndef RMATRIX3_DO_NOT_CHECK_SIZE
		if ( m.row != m.col ) std::cerr << "RMatrix Trace : not square" << std::endl;
#endif
		TYPE tr = (TYPE)0.0, *tmp = m.element;
		int n = m.row;
		while ( n-- )
		{
			tr += *tmp;
			tmp += (m.row + 1);
		}
		return tr;
	}
	
	friend _rmatrix Diag(_rmatrix &m)
	{
		int n, i;
		_rmatrix re;
		if ( m.row == 1 || m.col == 1 )
		{
			n = _MAX(m.row, m.col);
			re.SetZero(n,n);
			for ( i = 0; i < n; i++ ) re.element[i*n+i] = m.element[i];
		} else 
		{
			n = _MIN(m.row, m.col);
			re.ReNew(n,1);
			for ( i = 0; i < n; i++ ) re.element[i] = m.element[i*m.row+i];
		}
		return re;
	}
	
	friend _rmatrix Diag(const _rmatrix &m)
	{
		int n, i;
		_rmatrix re;
		if ( m.row == 1 || m.col == 1 )
		{
			n = _MAX(m.row, m.col);
			re.SetZero(n,n);
			for ( i = 0; i < n; i++ ) re.element[i*n+i] = m.element[i];
		} else 
		{
			n = _MIN(m.row, m.col);
			re.ReNew(n,1);
			for ( i = 0; i < n; i++ ) re.element[i] = m.element[i*m.row+i];
		}
		return re;
	}

	//friend _rmatrix Zeros(int r, int c)
	//{
	//	_rmatrix re(r, c);
	//	int n = r * c;
	//	TYPE *_r = re.element;
	//	while ( n-- ) *(_r++) = (TYPE)0.0;
	//	return re;
	//}

	//friend _rmatrix Ones(int r, int c)
	//{
	//	_rmatrix re(r, c);
	//	int n = r * c;
	//	TYPE *_r = re.element;
	//	while ( n-- ) *(_r++) = (TYPE)1.0;
	//	return re;
	//}

	//friend _rmatrix Rand(int r, int c)
	//{
	//	//srand( (unsigned)time( NULL ) );
	//	_rmatrix re(r, c);
	//	int n = r * c;
	//	TYPE *_r = re.element;
	//	while ( n-- ) *(_r++) = (TYPE)drand(1.0);
	//	return re;		
	//}

	//friend _rmatrix Eye(int r, int c )
	//{
	//	_rmatrix re(r,c);
	//	int n = r * c;
	//	TYPE *_r = re.element;
	//	while ( n-- ) *(_r++) = (TYPE)0.0;
	//	
	//	c = _MIN(r++,c);
	//	_r = re.element;
	//	while ( c-- )
	//	{
	//		*_r = (TYPE)1.0;
	//		_r += r;
	//	}		
	//	return re;
	//}

	//friend _rmatrix Eye(int r)
	//{
	//	return Eye(r,r);
	//}

};

typedef _rmatrix <gReal> RMatrix;
typedef _rmatrix <int> IMatrix;

template <class TYPE>
INDEXBASE _rmatrix <TYPE> ::idxbase = ZERO;

void tic();
gReal toc();
gReal drand(gReal range = 1.0);	// generates random value in [-range, range]
gReal prand(gReal range = 1.0);	// generates positive random value in [0, range]

RMatrix Zeros(int r, int c);
RMatrix Ones(int r, int c);
RMatrix Rand(int r, int c);
RMatrix Eye(int r, int c );
RMatrix Eye(int r);
RMatrix Inv(const RMatrix &m);

RMatrix Conv(const RMatrix &u, const RMatrix &v);
gReal Det(RMatrix A);
bool SolveAxEqualB(const RMatrix &A, RMatrix &x, const RMatrix &B);
bool SolveAtxEqualB(const RMatrix &A, RMatrix &x, const RMatrix &B);
bool SolvePosDefAxEqualB(RMatrix &A, RMatrix &x, const RMatrix &B);
int GaussElimination(RMatrix &A, IMatrix &row_pivot, IMatrix &column_pivot);
int GaussElimination(RMatrix &A, IMatrix &row_pivot, IMatrix &column_pivot, gReal eps);
RMatrix Eig(RMatrix m);
void Eig(RMatrix &re, RMatrix &m);
void Eig(RMatrix m, RMatrix &v, RMatrix &d);
RMatrix Companion(const RMatrix &m);
RMatrix Roots(const RMatrix& cof);
bool QRSolveAxEqualB(const RMatrix &A, RMatrix &x, const RMatrix &B);
bool QRSolveAtxEqualB(const RMatrix &A, RMatrix &x, const RMatrix &B);

// singular value decomposition of a general m by n matrix (the singular values are NOT sorted!)
RMatrix SVD(const RMatrix &M);										// return the unsorted singular values of M
void SVD(const RMatrix &M, RMatrix &U, RMatrix &S, RMatrix &V);		// M = U*Diag(S)*~V, (m,n) = (m,n)*(n,n)*(n,n), size of S = (n,1)

// singular value decomposition of a 3 by 3 matrix based on a closed-form solution (the singular values are SORTed!)
// SVD33() is faster than SVD(), but recommended for applications requiring moderate accuracy only! 
// implemented by Junggon Kim
void SVD33(gReal *M, gReal *S);													// input = M[9], output = S[3] = the sorted singular values of M
void SVD33(gReal *M, gReal *U, gReal *S, gReal *V, gReal tol = 0);				// input = M[9], output = U[9], S[3], V[9]
RMatrix SVD33(const RMatrix &M);													// return the sorted singular values of M
void SVD33(const RMatrix &M, RMatrix &U, RMatrix &S, RMatrix &V, gReal tol = 0);	// M is a 3x3 matrix, M = U*Diag(S)*~V, size(S) = (3,1)
																					// fabs(det(U))==1 may not be guaranteed especially when more than two singular values are (nearly) zeros.
																					// ==> increase tol to make fabs(det(U))==1 forcefully.
																					//     (This could decrease the accuracy of V and U even more and lead M != U*Diag(S)*~V.
																					//      So, use it at your own risk.)
																					// tol: the singular values less than tol will be considered as zeros 
																					//      when calculating corresponding columns of V and U internally.
																					//      default = 0 (internally 0 will be replaced with _EPS.)
																					//      This will not affect the singular values S.


gReal Cond(const RMatrix &m);
int Rank(const RMatrix &m, gReal singular_criteria);
int Rank(const RMatrix &m);

bool MultAbCt(RMatrix &M, const RMatrix &A, const RMatrix &b, const RMatrix &C);	// M = A*Diag(b)*~C
bool MultAbCt(RMatrix &M, const RMatrix &A, const gReal *b, const RMatrix &C);		// M = A*Diag(b)*~C
bool MultAb(RMatrix &M, const RMatrix &A, const RMatrix &b);						// M = A*Diag(b)
bool MultaB(RMatrix &M, const RMatrix &a, const RMatrix &B);						// M = Diag(a)*B

//-------------------------------------------------------
// Additional Linear Solver (added by Junggon Kim)
//-------------------------------------------------------

// (**) eps = tolerance ratio of a singular value with respect to the maximum singular value
//      Any singular values less than eps*_MAX(SVD) are treated as zero in calculating pseudo-inverse and null space.

// nullspace of A
RMatrix Nullspace(const RMatrix &A, gReal eps = 1E-9);							// return null space of A

// pseudo-inverse of A
RMatrix pInv(const RMatrix &A, gReal eps = 1E-9);								
RMatrix pInv(const RMatrix &A, RMatrix &N, gReal eps = 1E-9);					// N = orthonormal basis of the nullspace of A

// sr-inverse (singularity-robust inverse) of A with alpha > 0(See Y. Nakamura et al.)
// If alpha <= 0, pInv(...) will be called internally.
RMatrix srInv(const RMatrix &A, gReal alpha);									
RMatrix srInv(const RMatrix &A, RMatrix &N, gReal alpha, gReal eps = 1E-9);	// N = orthonormal basis of the nullspace of A

// solve A*x=b using pseudo-inverse of A and return rank of A
int solve_Ax_b_pInv(RMatrix &x, const RMatrix &A, const RMatrix &b, gReal eps = 1E-9);
int solve_Ax_b_pInv(RMatrix &x, const RMatrix &A, const RMatrix &b, RMatrix &N, gReal eps = 1E-9);
int solve_Ax_b_pInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, gReal eps = 1E-9);
int solve_Ax_b_pInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, RMatrix &Nw, gReal eps = 1E-9);

// solve A*x=b using sr-inverse of A
// If alpha <= 0, solve_Ax_b_pInv(...) will be called internally.
void solve_Ax_b_srInv(RMatrix &x, const RMatrix &A, const RMatrix &b, gReal alpha);
void solve_Ax_b_srInv(RMatrix &x, const RMatrix &A, const RMatrix &b, RMatrix &N, gReal alpha, gReal eps = 1E-9);
void solve_Ax_b_srInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, gReal alpha);
void solve_Ax_b_srInvW(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, RMatrix &Nw, gReal alpha, gReal eps = 1E-9);

// test functions
RMatrix srInv2(const RMatrix &A, gReal alpha);									// srInv() = srInv2()
void solve_Ax_b_srInv2(RMatrix &x, const RMatrix &A, const RMatrix &b, gReal alpha, RMatrix &s);
void solve_Ax_b_srInvW2(RMatrix &x, const RMatrix &A, const RMatrix &b, const RMatrix &iw, gReal alpha, RMatrix &s);

#endif