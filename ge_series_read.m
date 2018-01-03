function [ ge, lastfile ] = ge_series_read(examPath, series)

% ge_series_read - reads a volume of images from a GE series
% 
% [ ge, lastfile ] = ge_series_read(examPath, seriesPath)
% 
% examPath    - string path to exam directory
% seriesPath  - string path to series directory, 
%               relative to examPath
% 
% reads the volume of files in seriesPath, which is stored
% under examPath (it can return the name of the last file read).
% 
% The files are assumed to be arranged in series subdirectories
% below the examPath, so series files would be located as such:
% examPath/seriesN/*.MR or examPath/seriesN/I.*.  The function 
% will try to find all *.MR or I.* files and determine the 
% correct numerical order of the image files.
% 
% The function assumes the GE files are big endian.  This is
% an example of the returned struct, where ge.img contains
% the image volume data:
% 
% ge = 
% 
%      pix_hdr: [1x1 struct]
%    img.offset: 8432
%       su_hdr: [1x1 struct]
%       ex_hdr: [1x1 struct]
%       se_hdr: [1x1 struct]
%       hdr.image: [1x1 struct]
%          img: [256x256x11 double]
% 
% see also ge_hdr_read, ge_series2avw
% 


% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Licence: GPL, no express or implied warranties
% 
% Souheil J. Inati  <souheil.inati@nyu.edu> at 03/2003
% Dartmouth College, May 2000
% 
% Darren.Weber@flinders.edu.au, March 2003
% - Substantially redesigned file handling and function
%   call structures for integration with mri_toolbox at 
%   http://eeg.sf.net
% - Requested permission to distribute code under GPL licence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


version = '[$Revision: 1.1 $]';
fprintf('\nGE_SERIES_READ [v%s]\n',version(12:16));  tic;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Extract filenames and exam/series/image components

% check if series is int
if isnumeric(series),
    series = num2str(series);
end

seriesPath  = [examPath,filesep,series,filesep];

fprintf('...searching for files in: %s\n',seriesPath);

if exist(seriesPath) ~= 7,
    msg = sprintf('...cannot find path: %s!\n',seriesPath);
    error(msg);
end

dirFiles = dir(seriesPath);
n = 0;
for i = 1:length(dirFiles),
    % check whether this file is a directory
    if ~ dirFiles(i).isdir,
        % check whether it ends with .MR
        if findstr(dirFiles(i).name,'.MR'),
            %fprintf(sprintf('found *.MR file: %s\n',dirFiles(i).name));
            n = n + 1;
            seriesFileNames{n} = dirFiles(i).name;
            imageFormat = 'MR';
        elseif findstr(dirFiles(i).name,'I.'),
            %fprintf(sprintf('found I. file: %s\n',dirFiles(i).name));
            n = n + 1;
            seriesFileNames{n} = dirFiles(i).name;
            imageFormat = 'I';
        end
    end
end

if n == 0,
    msg = sprintf('...Found %d image files (I.* or *.MR)!\n',n);
    error(msg);
else
    fprintf('...Found %d image files\n',n);
end

if isequal(imageFormat,'MR'),
    
    fprintf('...sorting *.MR image files into numerical order\n');
    
    % extract exam/series/image numbers and sort into numerical order
    fileName.position.exam   = findstr(seriesFileNames{1},'E');
    fileName.position.series = findstr(seriesFileNames{1},'S');
    fileName.position.image  = findstr(seriesFileNames{1},'I');
    
    range = [ (fileName.position.exam + 1) : (fileName.position.series - 1) ];
    examN = seriesFileNames{1}(range);
    
    range   = [ (fileName.position.series + 1) : (fileName.position.image - 1) ];
    seriesN = seriesFileNames{1}(range);
    
    for i = 1:length(seriesFileNames),
        range = [ fileName.position.image + 1 ];
        imageNumbers{i} = seriesFileNames{i}(range:end);
        imageNumbers{i} = strrep(imageNumbers{i},'.MR','');
    end
    imageNumbers = str2double(imageNumbers);
    [n,i] = sort(imageNumbers);
    
    sortedFileNames = seriesFileNames(i);
    
elseif isequal(imageFormat, 'I'),
    
    fprintf('...sorting I.* image files into numerical order\n');
    
    % extract image numbers, which might already be in numerical order
    
    for i = 1:length(seriesFileNames),
        imageNumbers{i} = seriesFileNames{i}(3:end);
    end
    imageNumbers = str2double(imageNumbers);
    [n,i] = sort(imageNumbers);
    
    sortedFileNames = seriesFileNames(i);
    
end

% --- Create the name of the first file in examPath
firstfile = fullfile(seriesPath,sortedFileNames{1});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Read the header from the first file in the series
ge = ge_hdr_read(firstfile);



im_offset = ge.img.offset;  % this is not so good, as the ge.img is replaced below!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Read the image data from the GE series files

version = '[$Revision: 1.1 $]';
fprintf('\nGE_SERIES_READ [v%s]\n',version(12:16));
fprintf('...reading image data\n');

% initialize some variables
nX = ge.hdr.image.imatrix_X; % X Voxels
nY = ge.hdr.image.imatrix_Y; % Y Voxels
nZ = ge.hdr.image.slquant;   % Z Voxels (slice quantity)

sliceSize = nX*nY;
ge.img = zeros(nX, nY, nZ);

fprintf('...reading ');

for i = 1:nZ,
    
    imageFile = fullfile(seriesPath,sortedFileNames{i});
    
    % output filename to indicate progress
    if i == 1,
        backspaces = '';
    else
        filePreviousLength = length(fullfile(seriesPath,sortedFileNames{i-1}));
        backspaces = repmat('\b',1,filePreviousLength);
    end
    fprintf([backspaces,'%s'],imageFile);
    
    
    % Open the file
    [fid,message] = fopen(imageFile,'r','b'); % big endian
    if (fid == -1)
        fprintf('Cannot Open %s (big endian).\n',imageFile);
        
        % Can try to read little endian (shouldn't be necessary!)
        fprintf('Trying to read little endian\n');
        [fid,message] = fopen(imageFile,'r','l'); % little endian
        if (fid == -1),
            fprintf('Cannot Open %s (little endian).\n',imageFile);
            break
        end
    end
    
    % Skip the header, goto the data
    fseek(fid,im_offset,-1);
    % Read the slice data
    buffer = fread(fid,sliceSize,sprintf('int%d',ge.hdr.image.screenformat));
    % append the slice to the imageSet
    ge.img(:,:,i) = reshape(buffer, nX, nY);
    
    % Close the file
    status = fclose(fid);
    
    % Set the lastfile
    lastfile = imageFile;
end

t=toc; fprintf('\n...done (%5.2f sec).\n',t);

return
