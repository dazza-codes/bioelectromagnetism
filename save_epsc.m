function save_epsc(file,F,res,color)

% save_epsc - Save figure to EPS color file
% 
% Usage: save_epsc(file,F,res,color)
% 
% 'file'    - the string filename
% 'F'       - a figure handle, when empty uses 'gcf'
% 'res'     - resolution (dpi, 600 default)
% 'color'   - 0 = RGB, 1 = CMYK (default)
%
% This function uses the -painters option to render the output into a
% vector graphics format.  This may loose lighting/transparency effects in
% complex surfaces.
%
% Outputs Encapsulated Level 2 Color PostScript
% 

% $Revision: 1.3 $ $Date: 2008/05/04 18:23:43 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2004, Darren.Weber_at_radiology.ucsf.edu
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

if exist('color','var'),
    if isempty(color), color = 1; end
else
    color = 1;
end

fprintf('saving to:...%s\n',file);

driver = '-depsc2';
option1 = '-tiff'; % add tiff preview
option2 = '-noui';
option3 = '-painters'; % vector rendering, will loose lighting/transparency
option4 = '-loose'; % do not crop the image tightly
option5 = sprintf('-r%d',res);

if color,
    print(F,driver,option1,option2,option3,option4,option5,'-cmyk',file);
else
    print(F,driver,option1,option2,option3,option4,option5,file);
end

return
