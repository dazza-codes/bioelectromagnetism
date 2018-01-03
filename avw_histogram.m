function [bins,freq,freq_nozero] = avw_histogram(avw,binInterval,binNumber,plot)

% avw_histogram
%
% [bins,freq,freq_nozero] = avw_histogram(avw,binInterval,binNumber,plot)
%
% avw  - an Analyze 7.5 data struct, see avw_read
%
% binInterval - the intensity interval between bins (default = 5)
%
% binNumber - the total number of bins (default = 255)
%
% If binInterval is specified, it takes precedence over binNumber
%

% $Revision: 1.2 $ $Date: 2005/03/18 21:48:25 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2004, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '[$Revision: 1.2 $]';
fprintf('\nAVW_HISTOGRAM [v%s]\n',version(12:16));  tic;

if ~exist('binInterval','var'), binInterval = 5; end
if isempty(binInterval), binInterval = 5; end

if ~exist('binNumber','var'), binNumber = 255; end
if isempty(binNumber), binNumber = 255; end

if ~exist('plot','var'), plot = 1; end
if isempty(plot), plot = 1; end

% would be nice to use bins that reflect
% the bits per pixel, but it seems to take
% forever for a 16 bit image, otherwise
% use this code:
%     %check the bits per pixel of avw
%     bitpix = avw.hdr.dime.bitpix;
%     % set the bins according to the data type
%     if bitpix <= 8, bins = 0:255; end
%     if bitpix > 8,  bins = 0:65535; end

%bins = linspace(0,intensity_max,255);

intensity_min = min(min(min(avw.img)));
intensity_max = max(max(max(avw.img)));

% the determination of bins should incorporate binInterval and binNumber,
% for now, it only uses binInterval
bins = [intensity_min:binInterval:intensity_max]';

binNumber = length(bins);

fprintf('...calculating histogram for %d bins.\n',binNumber);

freq = zeros(1,binNumber);

for i=1:size(avw.img,3),
    % hist returns the counts for each bin for every column of the slice
    % input avw.img here.  Below, we sum across columns to get the
    % frequency count for the whole slice
    n = hist(avw.img(:,:,i), bins);
    if i == 1,
        freq = sum(n,2);
    else
        freq_slice = sum(n,2);
        freq = freq + freq_slice;
    end
end

if plot,
    figure('name','intensity histogram');
    bar(bins,freq);
end

freq_nozero = [freq(1) * 0; freq(2:end)];

return
