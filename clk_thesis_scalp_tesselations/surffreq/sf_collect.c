/*
 * sf_collect	    - collect indexed terms into a matrix
 *  	    	    - C MEX-file for MATLAB
 *  	    	    - started 7/97 by Clif Kussmaul
 
 * Copyright (c) 1997 Clif Kussmaul. All rights reserved.

 *** NOTES ***
 *

 *** THINGS TO DO ***
 - use to apply full(sparse) trick to sf_def_surf

 - options for duplicate values: std (save list of values)
 - use doubly-linked list (or tree) for speed
 - creates multidimensional arrays
 - allow complex values
 - handle vectors of values (x,y,z data for each index)

/***********************************************************************
 * declarations start here
 */

/*** include files ***/

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "matrix.h"
#include "mex.h"

#define MAX_STRING  20

#define max(A,B)    ((A>B)?(A):(B))
#define min(A,B)    ((A<B)?(A):(B))
// changed this from strcasecmp to strcmp to prevent compiler error - CLK 2004-02-23
#define streq(A,B)  (!strcmp(A,B))

typedef enum {
    MODE_MAX, MODE_MIN, MODE_NUM, MODE_SUM,
    MODE_FIRST, MODE_LAST, MODE_MEAN, MODE_STD
    } mode_t;

typedef struct row_s *row_p;
typedef struct row_s {
    double  val;
    int     num, r;
    row_p   next;
    } row_t;

typedef struct col_s *col_p;
typedef struct col_s {
    int     c;
    col_p   next;
    row_p   row;
    } col_t;


/***********************************************************************
 * supporting functions start here
 */

void free_row(row_p r)
{ if (r != NULL) { free_row(r->next); mxFree(r); } }

void free_col(col_p c)
{ if (c != NULL) { free_col(c->next); free_row(c->row); mxFree(c); } }

row_p	new_row(int r, row_p next, double val)
{
    row_p   nr;
    
    nr       = mxCalloc(1,sizeof(col_t));
    nr->r    = r;
    nr->next = next;
    nr->val  = val;
    nr->num  = 1;
    return(nr);
}

col_p	new_col(int c, col_p next, int row, double val)
{
    col_p   nc;
    
    nc       = mxCalloc(1,sizeof(col_t));
    nc->c    = c;
    nc->next = next;
    nc->row  = new_row(row,NULL,val);
    return(nc);
}


/***********************************************************************
 * gateway routine starts here
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char    	str1[MAX_STRING];
    double  	*arg, *idx, *val, *rslt, *rslt1,
    	    	*irp, *icp, *valp,
    	    	*ones;
    int     	*ir, *jc, *rnum, *rnum1,
    	    	idxM, idxN,
    	    	valM,
    	    	rsltM = 1, rsltN = 1,
    	    	i,
    	    	node_num,
    	    	col_flag = 1,
    	    	sparse_flag = 0;

    mode_t  	mode = MODE_SUM;
    col_p   	data,
    	    	col1, col2;   	
    row_p   	row1, row2;

    /*** check input arguments and get useful values ***/
    if ((nlhs < 1) || (nrhs < 2)) {
    	mexEvalString("help sf_collect"); return;
    	}
    if (mxGetNumberOfDimensions(prhs[0]) != 2)
    	mexErrMsgTxt("IDX must have two dimensions.");
    if (mxGetNumberOfDimensions(prhs[1]) != 2)
    	mexErrMsgTxt("VAL must have two dimensions.");
    if (mxIsComplex(prhs[1]))
    	mexErrMsgTxt("Complex arguments not allowed.");
    idx = mxGetPr(prhs[0]); idxM = mxGetM(prhs[0]); idxN = mxGetN(prhs[0]);
    val = mxGetPr(prhs[1]); valM = mxGetM(prhs[1])*mxGetN(prhs[1]);
    /*** if IDX is 1-D, treat as row vector ***/
    if (idxM == 1) idxM = idxN; idxN = 1;
    if (idxM != valM)
    	mexErrMsgTxt("IDX and VAL must have same number of values.");
    if (idxN == 1) {
    	col_flag = 1;
    	ones = mxCalloc(idxM,sizeof(double));
    	for (i=0; i<idxM; i++) ones[i] = 1;
    	}
    else if (idxN != 2) mexErrMsgTxt("IDX must have two (or one) columns.");
    	

    /*** process optional arguments ***/
    for (i=2; i<nrhs; i++) {
    	switch (mxGetClassID(prhs[i])) {
    	    case mxCHAR_CLASS:
    	    	if ((mxGetM(prhs[i]) != 1) || (mxGetN(prhs[i]) > MAX_STRING))
    	    	    mexErrMsgTxt("Oversize string argument.");
    	    	mxGetString(prhs[i],str1,mxGetN(prhs[i])+1);
    	    	/* ways to deal with optional arguments */
    	    	if      streq(str1, "max" ) 	mode = MODE_MAX;
    	    	else if streq(str1, "min" ) 	mode = MODE_MIN;
    	    	else if streq(str1, "num" ) 	mode = MODE_NUM;
    	    	else if streq(str1, "sum" ) 	mode = MODE_SUM;
    	    	else if streq(str1, "first")	mode = MODE_FIRST;
    	    	else if streq(str1, "last") 	mode = MODE_LAST;
    	    	else if streq(str1, "mean") 	mode = MODE_MEAN;
    	    	/* full or sparse result matrix */
    	    	else if streq(str1, "full")  	sparse_flag = 0;
    	    	else if streq(str1, "sparse") 	sparse_flag = 1;
    	    	else mexErrMsgTxt("Unknown mode.");
    	    	break;  
    	    case mxDOUBLE_CLASS:
    	    	if ((mxGetM(prhs[i]) == 1) && (mxGetN(prhs[i]) == 2)) {
    	    	    arg = mxGetPr(prhs[i]);
    	    	    rsltM = arg[0];
    	    	    rsltN = arg[1];
    	    	    }
    	    	else mexErrMsgTxt("Invalid numeric argument.");
    	    	break;
    	    default:	
    	    	mexErrMsgTxt("Invalid argument.");
            } /* end switch mxGetClassID */
    	} /* end for optional args */
    

    /*** get size of RESULT matrix if not specified ***/
    if (rsltM==1 && rsltN==1) {
    	for (irp=idx,      i=0; i<idxM; i++, irp++) {
    	    if      (*irp > rsltM)  rsltM = *irp;
    	    else if (*irp < 1)      mexErrMsgTxt("Row value < 1 in IDX.");
    	    }
    	if (col_flag) rsltN = 1;
    	else for (icp=idx+idxM, i=0; i<idxM; i++, icp++) {
    	    if      (*icp > rsltN)  rsltN = *icp;
    	    else if (*icp < 1)      mexErrMsgTxt("Col value < 1 in IDX.");
    	    }
    	} /* end get size of RESULT matrix ***/


    if (sparse_flag) {
    	/*** insert values in dynamic data structure ***/
    	irp = idx; 
    	icp = (col_flag)?ones:idx+idxM;
    	valp = val;
    	data = new_col(*icp, NULL, *irp, *valp);
    	irp++; icp++; valp++; 
    	node_num = 1;
    	for (; valp<val+valM; irp++, icp++, valp++) {
    	    /* search for column */
    	    for (col1 = data; 
    	    	 col1 != NULL && col1->c < *icp; 
    	    	 col2 = col1, col1 = col1->next);

    	    /* if column doesn't exist, create it */
    	    if (col1 == NULL || col1->c != *icp) {
    	    	col2->next = new_col(*icp, col1, *irp, *valp);
    	    	node_num++;
    	    	}
    	    
    	    /* if column does exist, search for row */
    	    else {
    	    	for (row1 = col1->row;
    	    	     row1 != NULL && row1->r < *irp;
    	    	     row2 = row1, row1 = row1->next);
    	    	
    	    	/* if row doesn't exist, create it */
    	    	if (row1 == NULL || row1->r != *irp) {
    	    	    row2->next = new_row(*irp, row1, *valp);
    	    	    node_num++;
    	    	    } 

    	    	/* if row does exist, update value */
    	    	else {
    	    	    switch (mode) {
    	    	    	case MODE_MAX: 
    	    	    	    row1->val = max(*valp, row1->val);
    	    	    	    break;
    	    	    	case MODE_MIN:
    	    	    	    row1->val = min(*valp, row1->val);
    	    	    	    break;
    	    	    	case MODE_NUM:
    	    	    	    row1->num++;
    	    	    	    break;
    	    	    	case MODE_SUM:
    	    	    	    row1->val += *valp;
    	    	    	    break;
    	    	    	case MODE_FIRST:
    	    	    	    break;
    	    	    	case MODE_LAST:
    	    	    	    row1->val = *valp;
    	    	    	    break;
    	    	    	case MODE_MEAN:
    	    	    	    row1->num++;
    	    	    	    row1->val += *valp;
    	    	    	    break;
    	    	    	default:
    		    	    mexErrMsgTxt("Unimplemented mode for sparse output.");
    		    	    break;
    	    	    	} /* end switch */
    	    	    } /* end update value */
    	    	} /* end search for row */

    	    } /* end for */
    	    
    	/*** convert data structure into sparse matrix ***/
    	plhs[0] = mxCreateSparse(rsltM,rsltN, node_num, mxREAL);
    	rslt = mxGetPr(plhs[0]);
    	  ir = mxGetIr(plhs[0]);
    	  jc = mxGetJc(plhs[0]);
    	node_num = 0;
    	col1 = data;
    	for (i=0; i<rsltN; i++) {
    	    *jc++ = node_num;
    	    if ((col1 == NULL) || (i < col1->c - 1)) continue;
    	    for (row1 = col1->row; row1 != NULL; row1=row1->next) {
    	    	  *ir++ = row1->r - 1;
    	    	switch (mode) {
    	    	    case MODE_NUM:
    	    	    	*rslt++ = row1->num;
    	    	    	break;
    	    	    case MODE_MEAN:
    	    	    	*rslt++ = row1->val / row1->num;
    	    	    	break;
    	    	    default:
    	    	    	*rslt++ = row1->val;
    	    	    	break;
    	    	    } /* end switch */
    	    	node_num++;
    	    	}
    	    col1 = col1->next;
    	    }
    	free_col(data);
    	} /* end if sparse */
    	

    else {
    	/*** create, initialize, and fill RESULT matrix ***/
    	plhs[0] = mxCreateDoubleMatrix(rsltM, rsltN, mxREAL);
    	rslt = mxGetPr(plhs[0]);
    	rnum = mxCalloc(rsltM*rsltN, sizeof(int));
    	irp = idx; 
    	icp = (col_flag)?ones:idx+idxM;

    	for (valp = val; valp < val+valM; valp++, irp++, icp++) {
    	    i = (*irp-1) + (*icp-1)*rsltM;
    	    switch (mode) {
    		case MODE_MAX:
    		    if (rnum[i] != 0)	rslt[i] = max(*valp, rslt[i]);
    		    else    	    	rslt[i] =     *valp;
    		    break;
    		case MODE_MIN:
    		    if (rnum[i] != 0) 	rslt[i] = min(*valp, rslt[i]);
    		    else    	    	rslt[i] =     *valp;
    		    break;
    		case MODE_NUM:
    		    break;
    		case MODE_SUM:
    		    rslt[i] += *valp;
    		    break;
    		case MODE_FIRST:
    		    if (rnum[i] == 0) 	rslt[i] = *valp;
    		    break;
    		case MODE_LAST:
    		    rslt[i] = *valp;
    		    break;
    		case MODE_MEAN:
    		    rslt[i] += *valp;
    		    break;
    		default:
    		    mexErrMsgTxt("Unimplemented mode for full output.");
    		    break;
    		} /* end switch */
    	    rnum[i]++;
    	    } /* end for */
    	
    	/*** perform post-processing if necessary ***/
    	switch (mode) {
    	    case MODE_NUM: 
		for (rslt1=rslt, rnum1=rnum, i=0; i<rsltM*rsltN; i++) 
		    *rslt1++ = *rnum1++;
    	    	break;
    	    case MODE_MEAN: 
		for (rslt1=rslt, rnum1=rnum, i=0; i<rsltM*rsltN; i++, rslt1++, rnum1++)
    	    	    if (*rnum1 > 0) *rslt1 = *rslt1 / *rnum1;
    	    	break;
    	    } /* end switch mode */
    	mxFree(rnum);
    	} /* end if full */

    if (col_flag) mxFree(ones);
    	
} /* end of mexFunction */


/***********************************************************************
 * that's all folks
 */

