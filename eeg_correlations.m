function [COR] = eeg_correlations(volt)

% eeg_correlations - Pearsons correlations
%
% This utility simply implements the CORRCOEF command for a voltage
% data matrix with columns of electrodes and rows of sample points
%
% USEAGE: [COR] = eeg_correlations(volt)
%
% The utility returns a full matrix of correlation coefficients formed
% from voltage data array X whose each row is a sample point, and each
% column is an electrode.
%
% See also corrcoef
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:33 $

% Licence:  GNU GPL, no implied or express warranties
% Created:  07/00, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COR = corrcoef(volt);

return
