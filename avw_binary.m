function avw = avw_binary(avw,thresh)

% avw_binary - return 0 for avw.img <= thresh, 1 otherwise
%
% Usage: avw = avw_binary(avw,thresh)
%
% avw is the Analyze volume returned by avw_read
% thresh is the threshold intensity value (default = 0)
%


if ~exist('thresh','var'), thresh = 0; end
if isempty(thresh), thresh = 0; end

fprintf('...binarising volume at %g threshold...',thresh); tic;
Vindex = find(avw.img <= thresh);
avw.img(Vindex) = 0;
Vindex = find(avw.img);
avw.img(Vindex) = 1;
t = toc; fprintf('done (%5.2f sec)\n',t);

return
