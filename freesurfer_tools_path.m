function [FStoolsPath] = freesurfer_tools_path

% freesurfer_tools_path - locate the freesurfer*.m installation path
%
% [FStoolsPath] = freesurfer_tools_path
%
% This function looks for the installation path of
% freesurfer_read_surf.m in the matlab path.  If not 
% found, update the matlab path (try addpath).
%


% $Revision: 1.1 $ $Date: 2005/08/14 21:24:19 $

% Licence:  GNU GPL, no implied or express warranties
% History:  09/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FStoolsPath = fileparts(which('freesurfer_read_surf'));
if isempty(FStoolsPath),
    msg = sprintf('Cannot find freesurfer_read_surf on the matlab path.\nTry the addpath command.\n\n');
    error(msg);
end

return
