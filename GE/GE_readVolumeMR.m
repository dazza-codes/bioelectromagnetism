function [imageVol] = GE_readVolumeMR(firstfile, volSize, depth, im_offset,volnum)
%
%GE_readVolume
% 
% [imageVol] = GE_readVolumeMR(firstfile, [nX nY nZ], depth, im_offset,volnum)
%
% reads volume no. volnum of data with the EXXXXSYYYYI1.MR as a first file
% volnum = 1 is default
%
% Souheil J. Inati
% Dartmouth College
% January 2002
% souheil.inati@dartmouth.edu
%

if (nargin < 4)
        error('Not enough input arguments.')
        return
end

if nargin < 5
  volnum = 1;  % default is LX
end

% initialize some variables
nX = volSize(1);
nY = volSize(2);
nZ = volSize(3);
sliceSize = nX*nY;
imageVol = zeros(nX, nY, nZ);
[path_stem,filestem,ext] = fileparts(firstfile);
filestem = filestem(1:end-1); % chop of the trailing 1.

for i = 1:nZ
  filenum = (volnum-1)*nZ + i;
  % Make the filenames
  imageFile = fullfile(path_stem,[filestem sprintf('%d.MR',filenum)]);
  
  % Open the file
  [fid,message] = fopen(imageFile,'r','b');
  if (fid == -1)
    fprintf('Cannot Open %s.\n',imageFile);
    break
  end
  
  % Skip to the data
  fseek(fid,im_offset,-1);
  % Read the slice
  buffer = fread(fid,sliceSize,sprintf('int%d',depth));
  % append the slice to the imageSet
  imageVol(:,:,i) = reshape(buffer, nX, nY);
  
  % Close the file
  status = fclose(fid);

end

return
