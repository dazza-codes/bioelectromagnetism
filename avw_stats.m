function [stat] = avw_stats(stat)

% avw_stats - calculate region of interest stats for avw data
%
% [stat] = avw_stats(stat)
%
% stat is a struct with the following fields:
% 
% stat.roi  - a region of interest from avw.img, see avw_roi.
%
% stat.function - cell array of strings, which can be any matlab
% descriptive stats function.  The defaults are:
% {'mean','std','min','max','median','sum'}
% 
% For example,
%
% stat.roi.value = avw.img(1:10,1:10,1:10);
% stat.function = {'mean'};
% 
% stat.value - output array of stats values, with length(stat.function).
% 
% see also avw_roi
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/2003, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


version = '[$Revision: 1.1 $]';
fprintf('\nAVW_STATS [v%s]\n',version(12:16));  tic;

if ~exist('stat','var'),
    msg = sprintf('...no input stat struct.\n');
    error(msg);
end
if ~isfield(stat,'function'),
  fprintf('...no input stat.function, using defaults.\n');
  stat.function = {'mean','std','min','max','median','sum'};
elseif isempty(stat.function),
  fprintf('...empty input stat.function, using defaults.\n');
  stat.function = {'mean','std','min','max','median','sum'};
end


stat.value = zeros(1,length(stat.function));

for i = 1:length(stat.function),
  
  tmp = feval(stat.function{i},stat.roi.value);
  
  while length(tmp) > 1,
    tmp = feval(stat.function{i},tmp);
  end
  stat.value(i) = tmp;
  
end % stat roi

stat

return
