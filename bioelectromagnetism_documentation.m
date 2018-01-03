
% bioelectromagnetism_documentation - script to use m2html for documentation

if exist('m2html'),
    location = which('eeg_toolbox');
    [path,file,ext] = fileparts(location);
    path = strrep(path,'eeg_toolbox','');
    cd(path);
    m2html('m','bioelectromagnetism','html','doc_m2html');
else
    error('cannot locate m2html, see http://www.madic.org/download/matlab/m2html/');
end

clear ext file location path
