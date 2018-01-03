function save_pdf(file,F)

% save_pdf - Save figure to PDF color file
% 
% Usage: save_pdf(file,F)
% 
% 'file'    - the string filename
% 'F'       - a figure handle, when empty uses 'gcf'
% 
% Outputs PDF
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

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

driver = '-dpdf';
option1 = '-noui';
%option3 = '-zbuffer';
fprintf('saving to:...%s\n',file);
print(F,driver,option1,file);

return
