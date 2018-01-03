function avw2brainstorm(avw,Segment,Scalp,PCS,Comment)

% avw2brainstorm - Convert Analyze struct into BrainStorm file
% 
% avw2brainstorm(avw,[segment],[scalp],[PCS],[comment])
% 
% This function saves a brainstorm file with a name that combines
% avw.fileprefix with _subjectimage.mat - for example, if
% avw.fileprefix is set to 'subject01', this function will
% output a brainstorm MRI file, in the current working
% directory, called 'subject01_subjectimage.mat'
% 
% avw       - Analyze data loaded by avw_img_read.
% segment   - a sparse 3D matrix the same size as avw.img, with
%             zeros (not brain) and other integer flags to indicate
%             different tissue types (as yet unknown values!)
% scalp     - an Nx3 series of vertex points on the scalp, with
%             coordinates in the Patient Coordinate System (PCS).
% PCS       - a struct with the fields described below
% comment   - a char array describing the image
% 
% PCS contains the following fields:
% 
% PCS.R % [3x3 double] rotations from PCS to MRI
% PCS.t % [3x1 double] translations from PCS to MRI
% PCS.Comment % a char array describing the PCS 
%             % orientation (always 'Neuromag 122' here)
% PCS.PCSFiducial  % [3x4 double] of electrode fiducials
% PCS.CubeFiducial % [3x4 double] of MRI fiducials
% PCS.FiducialName = {'Nasion'  'LeftEar'  'RightEar'  'Origin'};
% 
% see also AVW_VIEW & ELEC_COREGISTER to obtain fiducials and
% calculate T, where PCS.R = T(1:3;1:3) and PCS.t = T(4,1:3).
% 
% See http://neuroimage.usc.edu/ for more information about 
% brainstorm, including a download PDF of the file formats.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% Author:   10/2002, Darren.Weber@flinders.edu.au
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('AVW2BRAINSTORM...\n'); tic

if ~exist('avw','var'),
    msg = sprintf('...no input avw - see help avw2brainstorm\n');
    error(msg);
end

fprintf('...converting avw to brainstorm variables\n');

Cube = uint8(avw.img);

Voxsize = double(avw.hdr.dime.pixdim(2:4));

if ~exist('Comment','var'),
    if isfield(avw,'fileprefix'),
        Comment = avw.fileprefix;
    else
        Comment = '';
    end
end

if ~exist('Segment','var'), Segment = []; end

if ~exist('Scalp','var'), Scalp = []; end

if ~exist('PCS','var'),
	PCS.R = eye(3);     % [3x3 double] rotations
	PCS.t = zeros(3,1); % [3x1 double] translations
	PCS.Comment = 'Neuromag 122';
	PCS.PCSFiducial  = []; % [3x4 double]
	PCS.CubeFiducial = []; % [3x4 double]
	PCS.FiducialName = {'Nasion'  'LeftEar'  'RightEar'  'Origin'};
else
    PCS.Comment = 'Neuromag 122';
end

savename = sprintf('%s_subjectimage.mat',avw.fileprefix);

fprintf('...saving brainstorm variables to %s\n', savename);

save(savename,'Cube','Voxsize','Comment','Segment','Scalp','PCS');

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
