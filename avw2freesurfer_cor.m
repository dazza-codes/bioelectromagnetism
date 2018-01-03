function avw2cor(avw,CORpath)

% avw2cor - converts an avw struct to FreeSurfer COR-* files
% 
% avw2cor(avw,CORpath)
%
% avw - Analyze volume, see avw_read
%
% CORpath - full path to FreeSurfer directory
%
% This function will translate an 8-bit
% uchar AVW file to a series of MGH-style 
% COR files that also have an 8-bit uchar 
% datatype.
%
% Example: avw2cor(avw,'/data/subjects/bert/mri/orig')
% will use the workspace avw struct and output files
% named /data/subjects/bert/mri/orig/COR-[0-256]
%
% However, it is probably best to use the FreeSurfer
% command line tool, mri_convert, as it will handle
% various datatypes and it will reslice Analyze
% volumes that are not 256^3 mm FOV and 1^3 mm 
% voxels.
%

% An AFNI BRIK volume may be formed from the resulting 
% COR files by using the AFNI command:
% to3d -anat -prefix raw -session $SUBJECTS_DIR/$subjname/afni \
%      -view orig -datum byte -orient LSP \
%      $SUBJECTS_DIR/$subjname/mri/orig/COR-\[0-256]\*


% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  1/21/99: Timothy M. Ellmore, LBC/NIMH
%           Aug 2003, Darren.Weber_at_radiology.ucsf.edu
%                     Adapted to mri_toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '[$Revision: 1.1 $]';
fprintf('\nAVW2COR [v%s]\n',version(12:16));  tic;

if ~exist('avw','var'),
  doc avw2cor;
  msg = sprintf('...no input avw\n');
  error(msg);
elseif isempty(avw),
  doc avw2cor;
  msg = sprintf('...empty input avw\n');
  error(msg);
end

% MGH COR files have these dimensions
xdim = 256;
ydim = 256;
zdim = 256;

% datatype is always 8-bit uchar for MGH COR files
type = 'uchar';

% this function used to read in an avwfile, but
% it now takes a workspace avw struct, so the
% conversion process requires avw_read first
%avw = avw_read(AVWfile);


% check that avw.img is the right dimensions
dim = avw.hdr.dime.dim(2:4);
if (dim(1) == xdim) & (dim(2) == ydim) & (dim(3) == zdim),
    % OK we have the right dimensions
else
    %should try to reslice the volume, but not sure if
    %this works properly, as meshgrid is a strange beast.
    %[xi,yi,zi] = meshgrid(1:1:256,1:1:256,1:1:256);
    %reslice = interp3(avw.img,xi,yi,zi);
    error('...avw.img is not 256x256x256 voxels, reslice it.\n');
end


% extract sliceplanes and write COR files

for i = 1:ydim,
    
    % COR files are coronal slices, so we assume
    % the analyze file was read correctly by avw_read
    % and we use the ydim to loop over coronal slices
    
    % create zero sliceplane
    CORfile = zeros(xdim,zdim);
    
    % extract the current sliceplane
    CORfile(:,:) = avw.img(:,i,:);
    
    % create the CORfile name
    CORfname = [CORpath,filesep,sprintf('COR-%03d',i)];
    
    if i > 1,
        backspaces = repmat('\b',1,7);
    else
        backspaces = '';
    end
    fprintf([backspaces,'%s image.'],sprintf('COR-%03d',i));
    
    fid = fopen(CORfname, 'w');
    fwrite(fid, CORfile(:), type);
    fclose(fid);
    
end
