function status = GE_convertC(inDir,outStem,starttime,nvols)
%
% status = GE_convertC(inDir,outStem)
%
%GE_convertC
%
% Version 1.0
%
% Makes two spm format images from complex GE data
% ie two phases (real and imag, or abs and phase)
% interleaved.
%
% Use the FSL functions avwcomplex to put them together
% for use with prelude to unwrap the phase and then fugue
% to unwarp the epi data.
%
% Souheil J. Inati  
% Dartmouth College
% September 2001
% souheil.inati@dartmouth.edu
%

% Create the name of the first file in inDir
firstfile = fullfile(inDir,'I.001');

% Read the Header from the first file
[su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(firstfile);

% The image Dimensions
Nx =  im_hdr.imatrix_X; % X Voxels
Ny =  im_hdr.imatrix_Y; % Y Voxels
Nz =  im_hdr.slquant;   % Z Voxels
volSize = [Nx Ny Nz];

% Is the first image the first or the last?
if (se_hdr.end_loc - se_hdr.start_loc) > 0
  scandir = 1;
else
  scandir = -1;
end

% Compute the M matrix
M = GE_createSPMmat(im_hdr,scandir);

% Create the template SPM volume structure
V.fname = '';
V.dim = [volSize spm_type('int16')]; % int_16
V.mat = M;
V.pinfo = [1 0 0]';

% Compute the number of files
D = dir(inDir);
numfiles = length(D) - 3;
nimtot = numfiles/2/Nz;

% Loop over the volumes until there is no more looping to be done.
for nim = 1:nimtot  

  % Create two analyze file
  vol_str = sprintf('000%d',nim);
  vol_str = vol_str(length(vol_str)-3:length(vol_str));
  outName = [outStem sprintf('_r%s.img',vol_str)];
  V.fname = outName;
  R = spm_create_image(V);
  outName = [outStem sprintf('_i%s.img',vol_str)];
  V.fname = outName;
  I = spm_create_image(V);
    
  for j = 1:Nz
    
    % Compute the image no.
    realsliceno = (nim-1)*2*Nz+2*j-1;
    imagsliceno = realsliceno + 1;
    
    % Read in the Real Data
    fname = fullfile(inDir,sprintf('I.%03d',realsliceno));
    % Open the file
    [fid,message] = fopen(fname,'r','b');
    if (fid == -1)
      error(sprintf('Cannot Open %s.\n',fname));
    end
    % Skip to the data
    fseek(fid,im_offset,-1);
    % Read the slice
    realbuff = fread(fid,Nx*Ny,sprintf('int%d',pix_hdr.img_depth));
    realbuff = reshape(realbuff, [Nx Ny]);
    % Close the file
    status = fclose(fid);

    % Read in the Image Data
    fname = fullfile(inDir,sprintf('I.%03d',imagsliceno));
    % Open the file
    [fid,message] = fopen(fname,'r','b');
    if (fid == -1)
      error(sprintf('Cannot Open %s.\n',fname));
    end
    % Skip to the data
    fseek(fid,im_offset,-1);
    % Read the slice
    imagbuff = fread(fid,Nx*Ny,sprintf('int%d',pix_hdr.img_depth));
    imagbuff = reshape(imagbuff, [Nx Ny]);
    % Close the file
    status = fclose(fid);

    % Write out the data
    R = spm_write_plane(R,realbuff,j);    
    I = spm_write_plane(I,imagbuff,j);
    
  end
  
end

% Done
fprintf('\nConverted %d volumes. \n',nim);

return
