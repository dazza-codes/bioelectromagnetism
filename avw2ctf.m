function mri = avw2ctf(avw)

% mri = avw2ctf(avw)
%
% The purpose of this function is to convert an Analyze volume into a CTF
% .mri volume.  It currently requires that the Analyze volume is
% 256x256x256, 1x1x1 mm voxels.  The avw_read function can
% handle various Analyze orientations, so this function assumes that the
% avw.img volume has axial unflipped orientation (it always is when
% returned by avw_read, regardless of the format in the .img file). The
% returned mri struct in the matlab workspace can be saved using
% ctf_write_mri (available in the ctf module of the cvs repository, see
% http://eeg.sf.net/.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    - adapted from an appendex to CTF document
%                    MRIConverter.pdf, which is copied at the end of this
%                    function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mri = ctf_make_mri;

mri.file = [avw.fileprefix,'.mri'];


ver = '[$Revision: 1.1 $]';
fprintf('\nCTF_AVW2MRI [v%s]\n',ver(12:16));  tic;


% these checks for the volume dims/pixel size could be replaced with an
% interpolation function that could take any Analyze volume and convert it
% to 256x256x256, 1x1x1 mm
if avw.hdr.dime.dim(2) ~= 256,
    error('avw.hdr.dime.dim(2) ~= 256');
end
if avw.hdr.dime.dim(3) ~= 256,
    error('avw.hdr.dime.dim(3) ~= 256');
end
if avw.hdr.dime.dim(4) ~= 256,
    error('avw.hdr.dime.dim(4) ~= 256');
end

if avw.hdr.dime.pixdim(2) ~= 1,
    error('avw.hdr.dime.dim(2) ~= 256');
end
if avw.hdr.dime.pixdim(3) ~= 1,
    error('avw.hdr.dime.dim(3) ~= 256');
end
if avw.hdr.dime.pixdim(4) ~= 1,
    error('avw.hdr.dime.dim(4) ~= 256');
end


% mri.hdr.dataSize = 1 or 2 (bytes), 8 or 16 bits
if avw.hdr.dime.bitpix == 8,
    mri.hdr.dataSize = 1;
    mri.hdr.clippingRange = 255;
else
    mri.hdr.dataSize = 2;
    mri.hdr.clippingRange = 65536;
end

% This next step should always work correctly, given that avw_read always
% converts any Analyze orientation into axial unflipped in the matlab
% workspace variable avw.img; the axial unflipped orientation is 
% (+X left,  +Y anterior,  +Z superior); the ctf orientation is
% (+X right, +Y posterior, +Z inferior), the complete opposite.
temp = flipdim(avw.img,1);
temp = flipdim(temp,2);
temp = flipdim(temp,3);
mri.img = temp;

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
