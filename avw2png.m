function avw2png(avw,pngPath)

% avw2png - converts an avw struct to *.png files
% 
% avw2png(avw,pngPath)
%
% avw - Analyze volume, see avw_read
%
% pngPath - full path to FreeSurfer directory
%
% This function will translate an 8-bit
% uchar AVW file to a series of MGH-style 
% COR files that also have an 8-bit uchar 
% datatype.
%
% Example: avw2png(avw,'/data/subjects/bert/mri/orig')
% will use the workspace avw struct and output files
% named /data/subjects/bert/mri/orig/COR-[0-256]
%
% However, it is probably best to use the FreeSurfer
% command line tool, mri_convert, as it will handle
% various datatypes and it will reslice Analyze
% volumes that are not 256^3 mm FOV and 1^3 mm 
% voxels.
%

% An AFNI BRIK volume may be formed from the resulting 
% COR files by using the AFNI command:
% to3d -anat -prefix raw -session $SUBJECTS_DIR/$subjname/afni \
%      -view orig -datum byte -orient LSP \
%      $SUBJECTS_DIR/$subjname/mri/orig/COR-\[0-256]\*


% $Revision: 1.3 $ $Date: 2006/08/15 00:28:31 $

% Licence:  GNU GPL, no express or implied warranties
% History:  1/21/99: Timothy M. Ellmore, LBC/NIMH
%           Aug 2003, Darren.Weber_at_radiology.ucsf.edu
%                     Adapted to mri_toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '[$Revision: 1.3 $]';
fprintf('\nAVW2PNG [v%s]\n',version(12:16));
tic;

  if ~exist('avw','var'),
    doc avw2png;
    msg = sprintf('...no input avw\n');
    error(msg);
  elseif isempty(avw),
    doc avw2png;
    msg = sprintf('...empty input avw\n');
    error(msg);
  end



  % volume dimensions
  xdim = avw.hdr.dime.dim(2);
  ydim = avw.hdr.dime.dim(3);
  zdim = avw.hdr.dime.dim(4);
  
  % voxel dimensions
  xpixdim = avw.hdr.dime.pixdim(2);
  ypixdim = avw.hdr.dime.pixdim(3);
  zpixdim = avw.hdr.dime.pixdim(4);
  
  
  intensityMax = max(max(max(avw.img)));
  % try to adjust for large intensity extremes
  if intensityMax > 255,
    % got 16 bit data, at least
    adjust = 0.5;
  else
    adjust = 0.9;
  end
  clim = [0 (intensityMax * adjust)];


  % grayscale map
  map = gray(256);

  % use the zdim to loop over axial slices
  for i = 1:zdim,
    
    % create the .png file name
    [avwPath,avwFile,avwExt] = fileparts(avw.fileprefix);
    pngFile = sprintf('%s_%03d.png',avwFile,i);
    %fprintf('%s image.\n', pngFile);
    pngFile = fullfile(pngPath,pngFile);
    
    % the imwrite command gives a lot of trouble, so switching to
    % the alternate use of imagesc and saving the entire figure
    useImwrite = 0;
    if useImwrite,
      
      % extract the current sliceplane
      bitdepth = double(avw.hdr.dime.bitpix);
      if bitdepth == 8,
        pngImage = uint8(  avw.img(:,:,i) );
      elseif bitdepth == 16,
        pngImage = uint16( avw.img(:,:,i) );
      else
        pngImage = double( avw.img(:,:,i) );
      end
      
      imwrite(pngImage,map,pngFile,...
              'Author','UCSF DNL',...
              'Software','http://eeg.sf.net',...
              'InterlaceType','adam7',...
              'BitDepth', bitdepth);
      
      
    else,
      
      
      F = figure('Color','black','visible','off');
      
      pngImage = squeeze( avw.img(:,:,i) );
      
      %pcolor(pngImage'); colormap(gray); shading flat
      % OR
      %image('Cdata',pngImage','CDataMapping','scaled',...
      %    'XData',[0 xdim],'YData',[0 ydim]);
      % OR
      imagesc([0,xdim],[0,ydim],pngImage',clim);
      
      colormap(map)
      
      set(gca,'YDir','normal','XLimMode','manual','YLimMode','manual',...
              'ClimMode','manual','YColor',[1 1 1],'XColor',[1 1 1])

      titleString = sprintf('Axial (slice %03d)', i);
      title(titleString,'Color',[1 1 1])
      ylabel('Y')
      xlabel('(Right <<)  X  (>> Left)'); % default radiological orientation for Analyze
      
      axis image
      daspect([xpixdim,ypixdim,zpixdim])
      
      save_png(pngFile, F, 300)
      
      close(F)
      
    end
    
  end


  
  %    PNG-specific parameters
  %    -----------------------
  %    'Author'       A string
  % 
  %    'Description'  A string
  % 
  %    'Copyright'    A string
  % 
  %    'CreationTime' A string
  % 
  %    'Software'     A string
  % 
  %    'Disclaimer'   A string
  % 
  %    'Warning'      A string
  % 
  %    'Source'       A string
  % 
  %    'Comment'      A string
  % 
  %    'InterlaceType' Either 'none' or 'adam7'
  % 
  %    'BitDepth'     A scalar value indicating desired bitdepth;
  %                   for grayscale images this can be 1, 2, 4,
  %                   8, or 16; for grayscale images with an
  %                   alpha channel this can be 8 or 16; for
  %                   indexed images this can be 1, 2, 4, or 8;
  %                   for truecolor images with or without an
  %                   alpha channel this can be 8 or 16
  % 
  %    'Transparency' This value is used to indicate transparency
  %                   information when no alpha channel is used.
  %                   
  %                   For indexed images: a Q-element vector in
  %                     the range [0,1]; Q is no larger than the
  %                     colormap length; each value indicates the
  %                     transparency associated with the
  %                     corresponding colormap entry
  %                   For grayscale images: a scalar in the range
  %                     [0,1]; the value indicates the grayscale
  %                     color to be considered transparent
  %                   For truecolor images: a 3-element vector in
  %                     the range [0,1]; the value indicates the
  %                     truecolor color to be considered
  %                     transparent
  % 
  %                   You cannot specify 'Transparency' and
  %                   'Alpha' at the same time.
  % 
  %    'Background'   The value specifies background color to be
  %                   used when compositing transparent pixels.
  % 
  %                   For indexed images: an integer in the range
  %                     [1,P], where P is the colormap length
  %                   For grayscale images: a scalar in the range
  %                     [0,1]
  %                   For truecolor images: a 3-element vector in
  %                     the range [0,1]
  % 
  %    'Gamma'        A nonnegative scalar indicating the file
  %                   gamma
  % 
  %    'Chromaticities' An 8-element vector [wx wy rx ry gx gy bx
  %                   by] that specifies the reference white
  %                   point and the primary chromaticities 
  % 
  %    'XResolution'  A scalar indicating the number of
  %                   pixels/unit in the horizontal direction
  % 
  %    'YResolution'  A scalar indicating the number of
  %                   pixels/unit in the vertical direction
  % 
  %    'ResolutionUnit' Either 'unknown' or 'meter'
  % 
  %    'Alpha'        A matrix specifying the transparency of
  %                   each pixel individually; the row and column
  %                   dimensions must be the same as the data
  %                   array; may be uint8, uint16, or double, in
  %                   which case the values should be in the
  %                   range [0,1]
  % 
  %    'SignificantBits' A scalar or vector indicating how many
  %                   bits in the data array should be regarded
  %                   as significant; values must be in the range
  %                   [1,bitdepth]
  % 
  %                   For indexed images: a 3-element vector
  %                   For grayscale images: a scalar
  %                   For grayscale images with an alpha channel:
  %                     a 2-element vector
  %                   For truecolor images: a 3-element vector
  %                   For truecolor images with an alpha channel:
  %                     a 4-element vector
  % 
  %    In addition to these PNG parameters, you can use any
  %    parameter name that satisfies the PNG specification for
  %    keywords: only printable characters, 80 characters or
  %    fewer, and no leading or trailing spaces.  The value
  %    corresponding to these user-specified parameters must be a
  %    string that contains no control characters except for
  %    linefeed.
  %
