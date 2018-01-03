/*
 * sf_gen_face	- generate face array (point indices) for given edge array
 *		- C MEX-file for MATLAB
 *		- started 6/96 by Clif Kussmaul
 *		- revised 5/97 for MATLAB 5.0

 * Copyright (c) 1997 Clif Kussmaul. All rights reserved.

 * NOTES:
 * - this function takes an (N x 2) edge endpoint array
 *   and returns a (M x 3) face array of end point indices
 * - number of faces is related to number of edges and vertices
 *
 * THINGS TO DO:
 * - clockwise vertex ordering for computing normals
 * - flag to control output messages (cnrrently commented)
 *   - use cpu time instead of tic to avoid interfering with scripts
 */

/***********************************************************************
 * declarations start here
 */

/*** include files ***/

#include <stdlib.h>
#include <math.h>
#include "mex.h"

/*** defined functions ***/
#define hplane(p, p1, p2)   \
    ((((pnt_p[p -1] - pnt_p[p1-1]) * (pnt_p[p2+pnt_n-1] - pnt_p[p1+pnt_n-1])) >	\
      ((pnt_p[p2-1] - pnt_p[p1-1]) * (pnt_p[p +pnt_n-1] - pnt_p[p1+pnt_n-1]))) ? 1 : -1)


/***********************************************************************
 * gateway routine
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *edge_p, *face_p, *pnt_p;
    long    i, j, k, m,
	    v1, v2, v3, 
	    *count, *list_p, 
	    edge_n, face_n, pnt_n;
    
    /*** check arguments and get useful pointers ***/
    if (nrhs == 0) {
	mexEvalString("help sf_gen_face"); return;
    }
    if (nrhs != 2)
	mexErrMsgTxt("Two input arguments required");
    if (nlhs > 1)
	mexErrMsgTxt("One output argument required");
    if (mxGetN(prhs[0]) != 2)
	mexErrMsgTxt("Point array must have 2 columns");
    if (mxGetN(prhs[1]) != 2)
	mexErrMsgTxt("Edge array must have 2 columns");


/*
    printf("generating faces... \n");
*/
/*
    mexEvalString("tic");
*/
    pnt_n  = mxGetM( prhs[0]);
    pnt_p  = mxGetPr(prhs[0]);
    edge_n = mxGetM( prhs[1]);
    edge_p = mxGetPr(prhs[1]);
     count = mxCalloc(edge_n, sizeof(long));
    face_n = 0;
    list_p = mxCalloc(edge_n*6, sizeof(long));


    /*** check each edge for adjacent faces ***/
    for (i=0; i<edge_n; i++) {
	
	/*** see if edges i and j share a vertex ***/
	for (j=i+1; j<edge_n; j++) {
	    /* skip edges on two faces */
	    if ((count[i] >= 2) || (count[j] >= 2)) continue;
	    
	    if      (edge_p[i       ]    == edge_p[j       ]) {
		v1 = edge_p[i+edge_n]; v2 = edge_p[j+edge_n];
		v3 = edge_p[i];
	    }
	    else if (edge_p[i       ]    == edge_p[j+edge_n]) {
		v1 = edge_p[i+edge_n]; v2 = edge_p[j       ];
		v3 = edge_p[i];
	    }
	    else if (edge_p[i+edge_n]    == edge_p[j       ]) {
		v1 = edge_p[i       ]; v2 = edge_p[j+edge_n];
		v3 = edge_p[i+edge_n];
	    }
	    else if (edge_p[i+edge_n]    == edge_p[j+edge_n]) {
		v1 = edge_p[i       ]; v2 = edge_p[j       ];
		v3 = edge_p[i+edge_n];
	    }
	    else continue;
	    
	    /*** see if edge k shares vertices with edges i and j ***/
	    for (k=j+1; k<edge_n; k++) {
		/* skip edges on two faces */
		if (count[k] >= 2) continue;
		
		if (((edge_p[k] == v1) && (edge_p[k+edge_n] == v2)) ||
		    ((edge_p[k] == v2) && (edge_p[k+edge_n] == v1))) {
		    /*** points v1,v2,v3 form a triangle, so ***/
		    /*** check for interior vertices ***/
		    for (m=0; m<pnt_n; m++) {
			if ((m != v1) && (m != v2) && (m != v3) &&
			    (hplane(m, v2, v3) == hplane(v1, v2, v3)) &&
			    (hplane(m, v1, v3) == hplane(v2, v1, v3)) &&
			    (hplane(m, v1, v2) == hplane(v3, v1, v2))) {
			    break;
			}
		    }
		    if (m >= pnt_n) {
			/*** points v1,v2,v3 form a face ***/
			list_p[3*face_n  ] = v1;
			list_p[3*face_n+1] = v2;
			list_p[3*face_n+2] = v3;
			count[i]++; count[j]++; count[k]++;
			face_n++;
			continue;
			}    
		}
	    } /* end of for k */
	} /* end of for j */
    } /* end of for i */
    
    /*** create Matrix of proper size and copy faces into it ***/
    plhs[0] = mxCreateDoubleMatrix(face_n, 3, 0);
    face_p = mxGetPr(plhs[0]);
    
    for (i=0; i<face_n; i++) {
	face_p[i         ] = list_p[3*i];
	face_p[i+face_n  ] = list_p[3*i+1];
	face_p[i+face_n*2] = list_p[3*i+2];
    }    
    
    mxFree(list_p);
/*
    mexEvalString("disp(['  face computation time = ', num2str(toc), ' secs'])");
*/

} /* end of mexFunction */


/***********************************************************************
 * that's all folks
 */

