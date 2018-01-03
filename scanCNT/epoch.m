% epochedData = epoch(data,srate,events,epoch,sortCode)
% 
%   ------------- IN ------------------------------------------------------------
%   data                    ->  matrix of CNT data (elecs in columns)
%   events                  ->  matrix of event info (from loadCNTdata):
%                                   column 1 = event stim type
%                                   column 2 = offset in points              
%   epoch                   ->  2-element vector containing start and stop points
%   sortCode (optional)     ->  scalar event code to select epoch
%
%   ------------- OUT -----------------------------------------------------------
%   epochedData             <-  multidimensional array (points x elecs x epoch)
%
%   -- Note: Works only with Scan 4.1+ data files

function out = epoch(data,events,epoch,sortCode)

if (nargin == 4)
    ind = find(events(:,1) == sortCode);
    offset = events(ind,2);
elseif (nargin == 3)
    offset = events(:,2);
end

start = offset + epoch(1);
stop = offset + epoch(2);

for i = 1:length(offset)
    out(:,:,i) = data(start(i):stop(i),:);
end