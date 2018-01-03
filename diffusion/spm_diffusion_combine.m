function spm_diffusion_combine(imgs1,imgs2)
% spm_diffusion_combine: part of the SPM diffusion tensor toolbox
% Russ Poldrack, 11/19/00
%
% - This routine combines data from two separate DTI sequences acquired
%	sequentially.  A "between-scan interleaving" procedure, in which
%	contiguous slices are acquired in separate scans, seems to increase
%	SNR compared to within-scan interleaving.  In the procedure used at
%	MGH, two scans are acquired with a distance factor (i.e., skip) 
%	equal to the slice thickness, and the z offset is incremented by
%	the slice thickness between scans.  
%   Arguments:
%   imgs1: names of .ima images for the first series (in the format returned
%	by spm_get)
%   imgs2: names of .ima images for the second series
%
%   Effects:
%   Creates a renumbered set of .ima images, which will be read correctly
%   by spm_diffusion_convert()
%   
% CVS repository and development information can be found at:
% http://spm-toolbox.sourceforge.net 
%
% The SPM diffusion toolbox is distributed freely under the GNU Public License.
%

[ima_dir fname ext]=fileparts(imgs1(1,:));
ser1=sprintf('%s%s',fname,ext);
nimgs=size(imgs1,1);

[ima_dir fname ext]=fileparts(imgs2(1,:));
ser2=sprintf('%s%s',fname,ext);
nimgs=nimgs+size(imgs2,1);

% images are laid out as follows:
% acquisition (changes most quickly)
% diffusion weighting 
% slice  (changes most slowly)
 
tmp=sscanf(ser1,'%d-%d-%d.ima');
exam(1)=tmp(1);
series(1)=tmp(2);
first_img(1)=tmp(3);

tmp=sscanf(ser2,'%d-%d-%d.ima');
exam(2)=tmp(1);
series(2)=tmp(2);
first_img(2)=tmp(3);
filenames=[ser1 ser2];

hdr1=read_ima_img(ser1);
nslices1=hdr1.ImaHdr.h_G21_Rel2_Mr_NumberOfSlicesCur;
hdr2=read_ima_img(ser2);
nslices2=hdr2.ImaHdr.h_G21_Rel2_Mr_NumberOfSlicesCur;
nslices=nslices1+nslices2;
navg=hdr1.ImaHdr.h_G18_Acq_NumberOfAverages;

% find which one is on top
% mroe negative = higher up

if (hdr2.ImaHdr.h_G21_Rel1_CM_ImagePosition_Tra < ...
    hdr1.ImaHdr.h_G21_Rel1_CM_ImagePosition_Tra)
   % series 2 is on top - reverse the indices
   filenames=fliplr(filenames);
   exam=fliplr(exam);
   series=fliplr(series);
   first_img=fliplr(first_img);
end;

slice=1;
set=1;
counter=[0 0 1];
spm('Pointer','Watch');
spm_progress_bar('Init',nimgs,'converting .ima files...');

while slice<=nslices,

 for dwt=1:7,
  for acq=1:navg,
   cmd=sprintf('cp -f %d-%d-%d.ima %d-DTI-%d.ima',...
       exam(set),series(set),first_img(set)+counter(set),exam(set),counter(3));
   unix(cmd);
   spm_progress_bar('Set',counter(3));
   counter(3)=counter(3)+1; 
   counter(set)=counter(set)+1; 
%   fprintf('%s\n',cmd);
  end;   
 end; 
 if (set==1),set=2;else,set=1;end; 
 slice=slice+1;
end;
  

spm('Pointer','Arrow');
spm_progress_bar('Clear');


