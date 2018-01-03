function save_png(file, F, res)

% save_png - Save figure to PNG file
% 
% Usage: save_png(file,F,res)
% 
% 'file'    - the string filename
% 'F'       - a figure handle, when empty uses 'gcf'
% 'res'     - resolution, 600 dpi default
% 
% Outputs Portable Network Graphics (600 dpi)
% 

% $Revision: 1.3 $ $Date: 2006/08/15 00:11:50 $

% Licence:  GNU GPL, no implied or express warranties
% History:  02/2003, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('F','var'),
    if isempty(F), F = gcf; end
else
    F = gcf;
end
if exist('file','var'),
    if isempty(file), error('no filename'); end
else
    error('no filename');
end
if exist('res','var'),
    if isempty(res), res = 600; end
else
    res = 600;
end



driver = '-dpng';
option1 = '-noui';
option2 = sprintf('-r%d',res);
option3 = '-opengl';
fprintf('saving to:...%s\n',file);
print(F,driver,option1,option2,option3,file);

return
