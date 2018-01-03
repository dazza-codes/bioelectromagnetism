function spm_diffusion_convert(startima,writeall)
% spm_diffusion_convert: part of the SPM diffusion tensor toolbox
% Russ Poldrack, 8/15/00
%
% - This routine converts Siemens .ima files into Analyze-format
%   (.img) files for subsequent SPM analysis.  It assumes that
%   files are laid out as follows
%      acquisition (changes most quickly)
%      diffusion weighting (7 weightings)
%      slice  (changes most slowly)
%   Based upon the first image, the number of images in the series is
%   counted and the number of acquisitions are determined based upon
%   the number of files (with the number of slices determined based 
%   upon the Siemens header).  
%   NOTE: This particular ordering is specific to the particular
%   pulse sequence used at the MGH-NMR Center.
%
% Arguments:
% startima: the filename of the first .ima image in the series
% writeall: if set to 0, only averaged images (over acquisitions)
%   will be created (defaults to 1, creating images for each acquisition
%   in addition to averaged images)
%
% Effects:
% Created .img files for each acquisition X diffusion weighting, along 
% with an average across acquisitions
%
% CVS repository and development information can be found at:
% http://spm-toolbox.sourceforge.net 
%
% The SPM diffusion toolbox is distributed freely under the GNU Public License.
%


if isempty(writeall),writeall=1;end;

n_weightings=7;
tmp=sscanf(startima,'%d-DTI-%d.ima');
exam=tmp(1);
first_img=tmp(2);
stem='dwitest';

hdr=read_ima_img(startima);
nacq=hdr.ImaHdr.h_G18_Acq_NumberOfAverages;
% nslices=hdr.ImaHdr.h_G21_Rel2_Mr_NumberOfSlicesCur;
teststem=sprintf('%d-DTI*',exam);
ima_files=dir(teststem);
nfiles=size(ima_files,1);
nslices=nfiles/(nacq*n_weightings);
min_loc=-10000;
max_loc=10000;
avg=zeros(n_weightings,hdr.DIM(1),hdr.DIM(2),nslices);
prog_counter=0;
spm_progress_bar('Init',n_weightings*nacq);
spm('Pointer','Watch');

for dwt=0:n_weightings-1,
    for acquisition=0:nacq-1,
      tmpvol=zeros(hdr.DIM(1),hdr.DIM(2),nslices);
      tmp2=zeros(hdr.DIM(1),hdr.DIM(2),nslices);
        
      for slices=0:nslices-1

  %NOTE - currenlty the slice thickness is incorrect because
  % it does not include the gap

        index=slices*nacq*n_weightings + dwt*nacq + acquisition + first_img;
        imgfile=sprintf('%d-DTI-%d.ima',exam,index);
%fprintf('dwt=%d,acq=%d,slice=%d: %s\n',dwt,acquisition,slices,imgfile);
        tmpdata=read_ima_img(imgfile);
        tmpvol(:,:,slices+1)=tmpdata.V;
        if tmpdata.ImaHdr.h_G21_Rel1_CM_ImagePosition_Tra<max_loc , max_loc=tmpdata.ImaHdr.h_G21_Rel1_CM_ImagePosition_Tra;end;
        if tmpdata.ImaHdr.h_G21_Rel1_CM_ImagePosition_Tra>min_loc , min_loc=tmpdata.ImaHdr.h_G21_Rel1_CM_ImagePosition_Tra;end;


      end;
      VOX=tmpdata.VOX;
      VOX(3)=tmpdata.ImaHdr.h_G18_Acq_SliceThickness;
      imgname=sprintf('%s-%d-%d.img',stem,dwt+1,acquisition+1); 
      DIM=tmpdata.DIM;
      DIM(3)=nslices;
      ORIGIN=tmpdata.ORIGIN;
      ORIGIN(3)=nslices/2;  
      tmpvol=flipdim(tmpvol,3);
      tmpvol=flipdim(tmpvol,2);
      tmpvol=flipdim(tmpvol,1);
      tmp2(:,:,:)=avg(dwt+1,:,:,:);
      tmp2=tmpvol+tmp2;
      avg(dwt+1,:,:,:)=tmp2;
      if writeall,
	  spm_hwrite(imgname,DIM,VOX,1,4,0,ORIGIN,'DWI');
    	  fid=fopen(imgname,'w'); fwrite(fid,tmpvol,'int16'); fclose(fid);
    	  fprintf('Writing %s: %d %d %d %d %d %d %d %d%d\n',imgname,DIM(1),DIM(2),DIM(3),ORIGIN(1),ORIGIN(2),ORIGIN(3),VOX(1),VOX(2),VOX(3));
      else,
	fprintf('Creating weighting #%d/%d, acq. #%d/%d\n',dwt+1,n_weightings,acquisition+1,nacq);
      end; %writeall
    prog_counter=prog_counter+1;
    spm_progress_bar('Set',prog_counter);
    end; %nacq
end; %dwt


% create average images


avg=avg/nacq;
fprintf('Writing averaged images...\n');
for dwt=1:n_weightings,
      imgname=sprintf('%s-%d-avg.img',stem,dwt);      
      spm_hwrite(imgname,DIM,VOX,1,4,0,ORIGIN,'DWI_AVG');
      fid=fopen(imgname,'w');
      fwrite(fid,avg(dwt,:,:,:),'int16');
      fclose(fid);
end;


%create average DWI image
grand_mean=zeros(hdr.DIM(1),hdr.DIM(2),nslices);
for dwt=2:n_weightings,
 grand_mean=grand_mean+squeeze(avg(dwt,:,:,:)); 
end;
grand_mean=fix(grand_mean/(n_weightings-1));

imgname=sprintf('%s-avg-dwi.img',stem);      
spm_hwrite(imgname,DIM,VOX,1,4,0,ORIGIN,'DWI_GRAND_AVG');
fid=fopen(imgname,'w');
fwrite(fid,grand_mean,'int16');
fclose(fid);


spm_progress_bar('Clear');

spm('Pointer','Arrow');




