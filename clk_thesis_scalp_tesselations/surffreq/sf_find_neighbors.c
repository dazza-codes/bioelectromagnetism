/*
 * sf_find_neighbors 	- find neighbors of src points in target points
 *  	    	    	- C MEX-file for MATLAB
 *  	    	    	- started 7/97 by Clif Kussmaul
 
 * Copyright (c) 1997 Clif Kussmaul. All rights reserved.

 * NOTES
 *
 
 * THINGS TO DO
 * - what to do if there aren't CNTMAX neighbors within DSTMAX?
 *   - currently returns index=0, dst=DSTMAX
 */

/***********************************************************************
 * declarations start here
 */

/*** include files ***/

#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "mex.h"

#define DEFAULT_DSTMAX	(100000.0)

/***********************************************************************
 * gateway routine starts here
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  	dstmax, dstuse,
    	    	dq, dx, dy, dz,
    	    	*dst,
    	    	*ncnt, *ndst, *nind, *npos,
    	    	       mx,  my,  mz,
    	    	*src,  sx,  sy,  sz, 
    	    	*tgt, *tx, *ty, *tz;
    int     	cntmax,
    	    	*ind,
    	    	src_nval, 
    	    	tgt_nval,
    	    	i,j,k, kmax, mcnt;
    
    /*** check input arguments and get useful values ***/
    if ((nrhs < 2) || (nlhs < 1)) {
    	mexEvalString("help sf_find_neighbors"); return;
    	}    
    if (mxGetNumberOfDimensions(prhs[0]) != 2)
    	mexErrMsgTxt("SRC must have two dimensions.");
    if (mxGetNumberOfDimensions(prhs[1]) != 2)
    	mexErrMsgTxt("TGT must have two dimensions.");
    if ((mxGetN(prhs[0]) != mxGetN(prhs[0])) ||
    	((mxGetN(prhs[0]) != 2) && (mxGetN(prhs[0]) != 3)))
    	mexErrMsgTxt("SRC and TGT must both have 2 or 3 columns.");
    
    src = mxGetPr(prhs[0]); src_nval = mxGetM(prhs[0]);
    tgt = mxGetPr(prhs[1]); tgt_nval = mxGetM(prhs[1]);

    /* get optional arguments */
    if (nrhs > 2) cntmax = mxGetScalar(prhs[2]); else cntmax = 1;
    if (nrhs > 3) dstmax = mxGetScalar(prhs[3]); else dstmax = DEFAULT_DSTMAX;

    /*** create output arguments and temporary storage ***/
    if (nlhs > 0) {
    	plhs[0] = mxCreateDoubleMatrix(src_nval, cntmax, mxREAL);
    	nind = mxGetPr(plhs[0]);
    	}
    if (nlhs > 1) {
	plhs[1] = mxCreateDoubleMatrix(src_nval, cntmax, mxREAL);
	ndst = mxGetPr(plhs[1]);
	}
    if (nlhs > 2) {
	plhs[2] = mxCreateDoubleMatrix(src_nval,      1, mxREAL);
	// commented this to remove compiler error - CLK 2004-02-23
//	mxSetLogical(plhs[2]);
	ncnt = mxGetPr(plhs[2]);
	}
    if (nlhs > 3) {
	plhs[3] = mxCreateDoubleMatrix(src_nval,      3, mxREAL);
	npos = mxGetPr(plhs[3]);
	}
    dst = mxCalloc(cntmax,sizeof(double));
    ind = mxCalloc(cntmax,sizeof(int   ));


    if (mxGetN(prhs[0]) == 2) {
	/*** search for neighbors of each 2D point ***/
	for (i=0; i<src_nval; i++) {
    	    sx = src[i]; 
    	    sy = src[i+src_nval]; 
    	    for (k=0; k<cntmax; k++) { ind[k] = 0; dst[k] = dstmax; }
    	    dstuse = dstmax;

    	    /*** search through all target points ***/
    	    tx = tgt; 
    	    ty = tgt + tgt_nval; 
    	    for (j=0; j<tgt_nval; j++, tx++, ty++) {
    		if ((fabs(dx = sx - *tx) < dstuse) && 
    	    	    (fabs(dy = sy - *ty) < dstuse) &&
    	    	    ((dq = dx*dx + dy*dy) < dstuse*dstuse)) {
    	    		/*** found valid neighbor, so update ncnt, ind, dst ***/
    	    		if (nlhs > 2) ncnt[i]++;
    	    		/* find index of most distant neighbor */
    	    		for (kmax=0,k=1; k<cntmax; k++) {
    	    	    	    if (dst[k] > dst[kmax]) kmax = k;
    	    	    	    }
    	    		/* update index and distance to neighbor */
    	    		ind[kmax] = j+1; dst[kmax] = sqrt(dq);
    	    		/* find distance to most distant neighbor */
    	    		for (dstuse=dst[0],k=1; k<cntmax; k++) {
    	    	    	    if (dst[k] > dstuse) dstuse = dst[k];
    	    	    	    }
    	    	    } /* end if */
    		} /* end for j */

    	    /*** save appropriate output arguments ***/
    	    if (nlhs > 0) for (k=0; k<cntmax; k++) nind[i+src_nval*k] = ind[k];
    	    if (nlhs > 1) for (k=0; k<cntmax; k++) ndst[i+src_nval*k] = dst[k];
    	    if (nlhs > 3) {
    		mcnt = 0; mx = my = 0;
    		for (k=0; k<cntmax; k++) {
    		    if (dst[k]<dstmax) {
    	    		mcnt++;
    	    		/* subtract 1 because MATLAB uses 1-based indexing */
    	    		mx += tgt[ind[k]-1           ];
    	    		my += tgt[ind[k]-1+tgt_nval  ];
    	    		}
    	    	    }
    		if (mcnt > 0) {
    	    	    npos[i           ] = mx/mcnt;
    	    	    npos[i+src_nval  ] = my/mcnt;
    	    	    }
    		}
    	    } /* end for i */
    	} /* if 2D points */

    else {
	/*** search for neighbors of each 3D point ***/
	for (i=0; i<src_nval; i++) {
    	    sx = src[i]; 
    	    sy = src[i+src_nval]; 
    	    sz = src[i+src_nval*2];
    	    for (k=0; k<cntmax; k++) { ind[k] = 0; dst[k] = dstmax; }
    	    dstuse = dstmax;

    	    /*** search through all target points ***/
    	    tx = tgt; 
    	    ty = tgt + tgt_nval; 
    	    tz = tgt + tgt_nval*2;
    	    for (j=0; j<tgt_nval; j++, tx++, ty++, tz++) {
    		if ((fabs(dx = sx - *tx) < dstuse) && 
    	    	    (fabs(dy = sy - *ty) < dstuse) &&
    	    	    (fabs(dz = sz - *tz) < dstuse) &&
    	    	    ((dq = dx*dx + dy*dy + dz*dz) < dstuse*dstuse)) {
    	    		/*** found valid neighbor, so update ncnt, ind, dst ***/
    	    		if (nlhs > 2) ncnt[i]++;
    	    		/* find index of most distant neighbor */
    	    		for (kmax=0,k=1; k<cntmax; k++) {
    	    	    	    if (dst[k] > dst[kmax]) kmax = k;
    	    	    	    }
    	    		/* update index and distance to neighbor */
    	    		ind[kmax] = j+1; dst[kmax] = sqrt(dq);
    	    		/* find distance to most distant neighbor */
    	    		for (dstuse=dst[0],k=1; k<cntmax; k++) {
    	    	    	    if (dst[k] > dstuse) dstuse = dst[k];
    	    	    	    }
    	    	    } /* end if */
    		} /* end for j */

    	    /*** save appropriate output arguments ***/
    	    if (nlhs > 0) for (k=0; k<cntmax; k++) nind[i+src_nval*k] = ind[k];
    	    if (nlhs > 1) for (k=0; k<cntmax; k++) ndst[i+src_nval*k] = dst[k];
    	    if (nlhs > 3) {
    		mcnt = 0; mx = my = mz = 0;
    		for (k=0; k<cntmax; k++) {
    		    if (dst[k]<dstmax) {
    	    		mcnt++;
    	    		/* subtract 1 because MATLAB uses 1-based indexing */
    	    		mx += tgt[ind[k]-1           ];
    	    		my += tgt[ind[k]-1+tgt_nval  ];
    	    		mz += tgt[ind[k]-1+tgt_nval*2];
    	    		}
    	    	    }
    		if (mcnt > 0) {
    	    	    npos[i           ] = mx/mcnt;
    	    	    npos[i+src_nval  ] = my/mcnt;
    	    	    npos[i+src_nval*2] = mz/mcnt;
    	    	    }
    		}
    	    } /* end for i */
        } /* end if 3D points */

    mxFree(dst);
    mxFree(ind);
    
} /* end of mexFunction */


/***********************************************************************
 * that's all folks
 */

