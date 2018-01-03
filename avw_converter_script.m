% avw_converter_script - convert byte order of .img files

% This function will convert byte order of Analyze files
% This example uses SPM example data, so
% cd to that data directory before calling this script

stat = mkdir(pwd,'tmp');
if stat,
    files = {'T1' 'T2'}; % list of file prefixes
    Ibyte = 'ieee-be';
    Obyte = 'ieee-le';
    Iorient = []; % auto determine orientation
    Oorient = 0;  % output axial unflipped
    for f = 1:length(files),
        avw = avw_img_read(files{f},Iorient,Ibyte);
        % create new file in tmp subdirectory
        newfile = [pwd filesep 'tmp' filesep files{f}];
        % to write out header and image...
        avw_img_write(avw,newfile,Oorient,Obyte);
        % to write image only...
        % avw.hdr = [];
        % avw_img_write(avw,newfile,Oorient,Obyte);
        % to write header only, see avw_hdr_write
    end
end
