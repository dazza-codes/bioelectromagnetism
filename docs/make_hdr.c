/*
 * (c) Copyright, 1986-1994
 * Biomedical Imaging Resource
 * Mayo Foundation
 *
 * to compile:
 *
 *    cc -o make_hdr make_hdr.c
 *
 */
 
#include <stdio.h>
#include "dbh.h"

main(argc,argv) /* file x y z t datatype max min */
int argc;
char **argv;
{
    int i;
    struct dsr hdr;
    FILE *fp;
    static char DataTypes[9][12] = {"UNKNOWN", "BINARY", "CHAR", "SHORT", "INT",
				    "FLOAT", "COMPLEX", "DOUBLE", "RGB"};
    							    
    static int DataTypeSizes[9] = {0,1,8,16,32,32,64,64,24};
    
    if(argc != 9)
    {
    	usage();
        exit(0);
    }
    memset(&hdr,0, sizeof(struct dsr));
    for(i=0;i<8;i++)
    	hdr.dime.pixdim[i]=0.0;
   
    hdr.dime.vox_offset = 0.0;
    hdr.dime.roi_scale   = 1.0;
    hdr.dime.funused1    = 0.0;
    hdr.dime.funused2    = 0.0;
    hdr.dime.cal_max     = 0.0;
    hdr.dime.cal_min     = 0.0;
  
    
    hdr.dime.datatype = -1;

    for(i=1;i<=8;i++)
    	if(!strcmp(argv[6],DataTypes[i]))
    	{
    		hdr.dime.datatype = (1<<(i-1));
    		hdr.dime.bitpix = DataTypeSizes[i];
    		break;
    	}
    		
    if(hdr.dime.datatype <= 0)
    {
    	printf("<%s> is an unacceptable datatype \n\n", argv[6]);
    	usage();
        exit(0);
    }
 
    if((fp=fopen(argv[1],"w"))==0)
    {
        printf("unable to create: %s\n",argv[1]);
        exit(0);
    }

    hdr.dime.dim[0] = 4;  /* all Analyze images are taken as 4 dimensional */
    hdr.hk.regular = 'r';
    hdr.hk.sizeof_hdr = sizeof(struct dsr);

    hdr.dime.dim[1] = atoi(argv[2]);  /* slice width  in pixels */
    hdr.dime.dim[2] = atoi(argv[3]);  /* slice height in pixels */
    hdr.dime.dim[3] = atoi(argv[4]);  /* volume depth in slices */
    hdr.dime.dim[4] = atoi(argv[5]);  /* number of volumes per file */

    hdr.dime.glmax  = atoi(argv[7]);  /* maximum voxel value  */
    hdr.dime.glmin  = atoi(argv[8]);  /* minimum voxel value */
    
/*	Set the voxel dimension fields: 
       A value of 0.0 for these fields implies that the value is unknown.
         Change these values to what is appropriate for your data
         or pass additional command line arguments     */      
         
    hdr.dime.pixdim[1] = 0.0; /* voxel x dimension */
    hdr.dime.pixdim[2] = 0.0; /* voxel y dimension */
    hdr.dime.pixdim[3] = 0.0; /* pixel z dimension, slice thickness */
    
/*   Assume zero offset in .img file, byte at which pixel
       data starts in the image file */

    hdr.dime.vox_offset = 0.0; 
    
/*   Planar Orientation;    */
/*   Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
     Movie flag ON:  3 = transverse, 4 = coronal, 5 = sagittal  */  

    hdr.hist.orient     = 0;  
    
/*   up to 3 characters for the voxels units label; i.e. 
    	mm., um., cm.               */

    strcpy(hdr.dime.vox_units," ");
   
/*   up to 7 characters for the calibration units label; i.e. HU */

    strcpy(hdr.dime.cal_units," ");  
    
/*     Calibration maximum and minimum values;  
       values of 0.0 for both fields imply that no 
       calibration max and min values are used    */

    hdr.dime.cal_max = 0.0; 
    hdr.dime.cal_min = 0.0;

    fwrite(&hdr,sizeof(struct dsr),1,fp);
    fclose(fp);
}

usage()
{
   printf("usage:  make_hdr name.hdr x y z t datatype max min \n\n");
   printf("  name.hdr = the name of the header file\n");
   printf("  x = width, y = height,  z = depth,  t = number of volumes\n");
   printf("  acceptable datatype values are: BINARY, CHAR, SHORT,\n");
   printf("                 INT, FLOAT, COMPLEX, DOUBLE, and RGB\n");
   printf("  max = maximum voxel value,  min = minimum voxel value\n");
}
