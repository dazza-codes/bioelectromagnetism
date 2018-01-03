/*
 * sf_find_vals     - find vals in larger list of values
 *  	    	    - C MEX-file for MATLAB
 *  	    	    - started 7/97 by Clif Kussmaul
 
 * Copyright (c) 1997 Clif Kussmaul. All rights reserved.

 * NOTES
 *
 
 * THINGS TO DO
 *
 
 */

/***********************************************************************
 * declarations start here
 */

/*** include files ***/

#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "mex.h"

/***********************************************************************
 * gateway routine starts here
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *vals, *list,
    	    *mtch,
    	    *idx;
    int     listM, listN,
    	    valsM, valsN,
    	    i1,i2,j1,j2;
    	    
    /*** check input arguments and get useful values ***/
    if ((nlhs < 1) || (nrhs < 2)) {
    	mexEvalString("help sf_find_vals"); return;
    	}
    if (mxGetNumberOfDimensions(prhs[0]) != 2)
    	mexErrMsgTxt("VALS must have two dimensions.");
    if (mxGetNumberOfDimensions(prhs[1]) != 2)
    	mexErrMsgTxt("LIST must have two dimensions.");
    valsM = mxGetM(prhs[0]);
    valsN = mxGetN(prhs[0]);
    listM = mxGetM(prhs[1]);
    listN = mxGetN(prhs[1]);
    if (valsN > listM)
    	mexErrMsgTxt("LIST cannot have fewer columns that VALS.");
    vals = mxGetPr(prhs[0]);
    list = mxGetPr(prhs[1]);


    /*** create temporary storage ***/
    mtch = mxCalloc(listM,sizeof(double));

    for (i1=0; i1<listM; i1++) {
    	for (j1=0; j1<valsM; j1++) {
    	    	for (j2=0; j2<valsN; j2++) {
    	    	    for (i2=0; i2<listN; i2++) {
    	    	    	/* if VALS[j1,j2] == list[i1,i2], try VALS[j1,j2+1] */
    	    	    	if (vals[j1+j2*valsM] == list[i1+i2*listM]) break;
    	    	    }
    	    	    /* if VALS[j1,j2] not in LIST[i1,:], try VALS[j1+1,:] */
    	    	    if (i2==listN) break;
    	    	    
    	    	}
    	    	/* if VALS[j1,:] in LIST[i1,:], save index value */
    	    	if (j2==valsN) mtch[i1] = j1+1;
    	    }
    	}
    

    /*** create output arguments ***/
    if (nlhs > 0) {
    	/*** count matches, allocate output, copy values ***/
    	for (i1=i2=0; i1<listM; i1++) if (mtch[i1] > 0) i2++;
    	plhs[0] = mxCreateDoubleMatrix(i2,1,mxREAL);
    	idx = mxGetPr(plhs[0]);
    	for (i1=i2=0; i1<listM; i1++) if (mtch[i1] > 0) idx[i2++] = i1+1;
    	}

} /* end of mexFunction */


/***********************************************************************
 * that's all folks
 */

