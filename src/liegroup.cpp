//====================================================================================================
//
//      Title :   liegroup.cpp
//
//      Note  :   This is a modified version of the original work by Jinwook Kim (see liegroup.h).
//                Modified and maintained by Junggon Kim (junggon@gmail.com)
//
//====================================================================================================

#include <math.h>
#include "liegroup.h"

int __idamax(int n, gReal *dx)
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

void __dgefa(gReal *x, int lda, int n, int *jpvt, int &info)
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
			l = __idamax(n-k, xk+k) + k;
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

void __dgesl(gReal *x, int lda, int n, int *jpvt, gReal *b, int job)
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

