% baselinedData = baselineData(epochedData,srate,epoch)
%                      
%   data                    ->  matrix or multidimensional matrix (elecs in columns)    
%   epoch (optional)        ->  2-element vector containing start and stop points
%
%   epochedData             <-  multidimensional array (points x elecs x epoch)
%
%   -- Note: Works only with Scan 4.1+ data files
function baselinedData = baselineData(epochedData,epoch)

if (nargin == 1)
    epoch = [1 size(epochedData,1)];
end

for i = 1:size(epochedData,3)
    meanData = mean(epochedData(epoch(1):epoch(2),:,i));
    meanData = meanData(ones(size(epochedData,1),1),:,ones(size(epochedData,3),1));
    baselinedData = epochedData - meanData;
end