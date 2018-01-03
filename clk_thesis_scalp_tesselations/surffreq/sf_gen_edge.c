/*
 * sf_gen_edge	- generate edge array (end point indices) for given locations
 *		- C MEX-file for MATLAB
 *		- started 6/96 by Clif Kussmaul
 *		- revised 5/97 for MATLAB 5.0

 * Copyright (c) 1997 Clif Kussmaul. All rights reserved.

 * NOTES:
 * - MATLAB manual says to use "gateway" and "computational" routines
 * - this function takes an (N x 2) intrinsic point array
 *   and returns a (M x 2) edge array of end point indices
 * - the output list is ordered by length, each edge is ordered by index
 * - MATLAB uses 1-based indexing,  C uses 0-based,  so add 1 at very end
 
 * THINGS TO DO:
 * - add jitter to avoid problems
 * - realloc edges as needed to minimize memory use
 * - wrapping on either/both axes for closed topologies?
 * - flag to control output messages (cnrrently commented)
 *   - use cpu time instead of tic to avoid interfering with scripts
 */
 
/***********************************************************************
 * declarations start here
 */

/*** include files ***/

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "mex.h"

/*** defined constants ***/

/*** global variables (for qsort) ***/
double	*x, *y;

/***********************************************************************
 * miscellaneous routines
 */

int qsort_func(a, b)    long *a, *b;
{
float	la, lb;

la = (x[a[0]] - x[a[1]]) * (x[a[0]] - x[a[1]]) +
     (y[a[0]] - y[a[1]]) * (y[a[0]] - y[a[1]]);
lb = (x[b[0]] - x[b[1]]) * (x[b[0]] - x[b[1]]) +
     (y[b[0]] - y[b[1]]) * (y[b[0]] - y[b[1]]);
return (la<lb) ? -1 : (la>lb);
}

/***********************************************************************
 * gateway routine
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *edge_p, *intr_p;
    float   ax, ay, bx, by, cx, cy, dx, dy, p, q, 
        // changed this from HUGE to FLT_MAX to fix compiler error - CLK 2004-02-23
        maxedge = FLT_MAX;
    long    i, j, k,
	    *list_p, *lp, 
	    edge_n, edge_np, intr_n;
    
    /*** check arguments and get useful pointers ***/
    if (nrhs == 0) {
	mexEvalString("help sf_gen_edge"); return;
    }
    if (nrhs == 0 || nrhs > 2)
	mexErrMsgTxt("One or two input arguments required");
    if (nlhs != 1)
	mexErrMsgTxt("One output argument required");
    if (nrhs == 2) {
	maxedge = *mxGetPr(prhs[1]);
	maxedge *= maxedge;
	}
    if (mxGetN(prhs[0]) != 2)
	mexErrMsgTxt("Input array must have 2 columns");

/*
    printf("generating edges... \n");
*/
/*
    mexEvalString("tic");
*/

    intr_n = mxGetM(prhs[0]);
    intr_p = mxGetPr(prhs[0]);
    x = &intr_p[0];
    y = &intr_p[intr_n];
	
    /*** make list of all possible edges and sort by length ***/
    lp = list_p = mxCalloc(intr_n*(intr_n - 1), sizeof(long));
    for (i=0; i<intr_n; i++) 
	for (j=i+1; j<intr_n; j++) {
	    if (((x[i] - x[j]) * (x[i] - x[j]) +
                 (y[i] - y[j]) * (y[i] - y[j])) < maxedge) {
		*lp++ = i; *lp++ = j; 
	    }
	}
    edge_np = (lp - list_p ) / 2;
    qsort(list_p, edge_np, 2*sizeof(long), qsort_func);

    /*** for each edge, add to triangulation if it doesn't intersect ***/
    edge_n = 0;
    for (i=0; i<edge_np; i++) {
	for (j=0; j<edge_n; j++) {
	    /* edges sharing endpoints can't intersect */
	    if ((list_p[2*i  ] == list_p[2*j  ]) ||
		(list_p[2*i  ] == list_p[2*j+1]) ||
		(list_p[2*i+1] == list_p[2*j  ]) ||
		(list_p[2*i+1] == list_p[2*j+1]))
		continue;
	    /* compute point of intersection */
	    ax = x[list_p[2*i  ]]; ay = y[list_p[2*i  ]];
	    bx = x[list_p[2*i+1]]; by = y[list_p[2*i+1]];
	    cx = x[list_p[2*j  ]]; cy = y[list_p[2*j  ]];
	    dx = x[list_p[2*j+1]]; dy = y[list_p[2*j+1]];
	    p = ((dx-cx)*(ay-cy) - (dy-cy)*(ax-cx)) /
		((bx-ax)*(dy-cy) - (by-ay)*(dx-cx));
	    q = ((bx-ax)*(cy-ay) - (by-ay)*(cx-ax)) /
		((dx-cx)*(by-ay) - (dy-cy)*(bx-ax));
	    if (1.0>p && p>0.0  &&  1.0>q && q>0.0) break;
	    
	}
	if (j == edge_n) {
	    list_p[2*j  ] = list_p[2*i  ]; 
	    list_p[2*j+1] = list_p[2*i+1];
	    edge_n++;
	}
    }
    
    /*** create Matrix of proper size and copy edges into it ***/
    plhs[0] = mxCreateDoubleMatrix(edge_n, 2, 0);
    edge_p = mxGetPr(plhs[0]);
    
    for (i=0; i<edge_n; i++) {
	edge_p[i]        = 1+ list_p[2*i];
	edge_p[i+edge_n] = 1+ list_p[2*i+1];
    }
    
    mxFree(list_p);
/*
    mexEvalString("disp(['  edge computation time = ', num2str(toc), ' secs'])");
*/
 
} /* end of mexFunction */


/***********************************************************************
 * that's all folks
 */

