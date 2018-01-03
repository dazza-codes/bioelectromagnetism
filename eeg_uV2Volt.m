function [Volts] = eeg_uv2volt(uV);

% eeg_uv2volt - Convert microvolts to volts
%
%   [Volts] = eeg_uv2volt(uV)
%
%   Simply, Volts = uV ./ 10^6
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:33 $

% Licence:  GNU GPL, no implied or express warranties
% History:  11/2001, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Volts = uV ./ 10^6;
