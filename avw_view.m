function [ varargout ] = avw_view(avw,parent,command),

% avw_view - create and navigate ortho views of Analyze 7.5 volume
%
% [ avw ] = avw_view([avw], [parent], [command])
%
% avw - a struct, created by avw_read; if omitted, a gui file locator will
% prompt for a .hdr file.
%
% parent  - an optional handle to the gui that calls this gui, useful for
% updating the UserData field of the parent. The avw structure may be
% returned to the parent, if possible.
%
% command - an optional string, such as 'init' or various callback
% commands, 'init' is the default.
%
% The navigation is by sliders, mouse clicks and arrow keys.  Right mouse
% clicks on any ortho-view will show a command menu, for simple block
% region of interest (ROI) calculations, image zoom, and save image.  The
% ROI calculations are returned into avw.stats.
%
% Fiducial points can be selected, which are returned into 'avw.fiducials'
% in the base workspace. These points are given in several coordinate
% frameworks (voxels, mm, and meters), with values given relative to an
% origin located at the "center" of the MRI volume (see avw_center, which
% returns abs and corner values, abs used here).
%
% The AC location can be selected and the values are returned into 'avw.ac'
% in the base workspace. These points are given in voxels, mm & meters; for
% the latter, the values are given as offsets from the "center" of the MRI
% volume (see avw_center).
%
% +X is left (L), +Y is anterior (A), +Z is superior (S), the default LAS
% orientation of the Analyze 7.5 format.  The coordinate system is left
% handed.  This is the radiological convention, as opposed to the
% neurological convention (RAS).  The latter can be emulated by using the
% 'Flip L/R' button.
%
% Example of loading and viewing the SPM T1 template:
% avw = avw_read('T1')
% avw = avw_view(avw);
%
% Similarly, just 'avw_view' can be typed at the command prompt and you can
% use the gui file locator to select any .hdr file.
%
% See also, avw_read, avw_img_read, avw_hdr_read
%

% $Revision: 1.5 $ $Date: 2007/05/09 00:06:58 $

% Licence:  GNU GPL, no express or implied warranties
% History:  06/2002, Darren.Weber_at_flinders.edu.au
%           10/2002, Darren.Weber_at_flinders.edu.au
%                    added fiducial point determination
%                    changed plots from surf to imagesc commands
%                    added handling of datatype for avw.img
%                    altered daspect to use avw.hdr.dime.pixdim
%                    altered color scheme
%           01/2003, Darren.Weber_at_flinders.edu.au
%                    added parent GUI handling
%           10/2003, Darren.Weber_at_radiology.ucsf.edu
%                    added right click options, including simple block ROI
%                    functions, zoom and save image
%           11/2003, Darren.Weber_at_radiology.ucsf.edu
%                    added arrow key navigation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gcbf,
  AVWVIEW = get(gcbf,'Userdata');
end

if ~exist('command','var'), command = 'init'; end

command = lower(command);

% Check for specific keys and assign command
if strcmp(command, 'keypress'),

  cc = get(AVWVIEW.gui,'CurrentCharacter');
  cc = double(cc);
  if cc,
    switch cc,
     case 27, command = 'quit';  % ESC
     case 28, command = 'left';  % left
     case 29, command = 'right'; % right
     case 30, command = 'up';    % up
     case 31, command = 'down';  % down
     otherwise, return;  % all other keys
    end
  end
end

switch command,

 case 'init',

  if ~exist('avw','var'),
    avw = avw_read;
  end

  if nargin == 0,
    AVWVIEW = init(avw);
  elseif isempty(inputname(1)),
    AVWVIEW = init(avw);
  else
    AVWVIEW = init(avw,inputname(1));
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'coordinates',

  AVWVIEW = set_coordinates(AVWVIEW);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'crosshairs',

  AVWVIEW = set_coordinates(AVWVIEW);
  %AVWVIEW = set_crosshairs(AVWVIEW);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case {'axial_image','coronal_image','sagittal_image'},

  switch command,
   case 'axial_image',
    AVWVIEW.view = 'axi'; axi_update = 0; cor_update = 1; sag_update = 1;
   case 'coronal_image',
    AVWVIEW.view = 'cor'; axi_update = 1; cor_update = 0; sag_update = 1;
   case 'sagittal_image',
    AVWVIEW.view = 'sag'; axi_update = 1; cor_update = 1; sag_update = 0;
  end

  AVWVIEW = get_current_position(AVWVIEW);

  if axi_update,
    axial_update(AVWVIEW);
  end
  if cor_update,
    coronal_update(AVWVIEW);
  end;
  if sag_update,
    sagittal_update(AVWVIEW);
  end;

  AVWVIEW = set_coordinates(AVWVIEW);
  %AVWVIEW = set_display_values(AVWVIEW);
  %AVWVIEW = set_crosshairs(AVWVIEW);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case {'axial_slider','coronal_slider','sagittal_slider'},

  switch command,

   case 'axial_slider',
    AVWVIEW.view = 'axi';
    AVWVIEW = get_slider_position(AVWVIEW);
    axial_update(AVWVIEW);

   case 'coronal_slider',
    AVWVIEW.view = 'cor';
    AVWVIEW = get_slider_position(AVWVIEW);
    coronal_update(AVWVIEW);

   case 'sagittal_slider',
    AVWVIEW.view = 'sag';
    AVWVIEW = get_slider_position(AVWVIEW);
    sagittal_update(AVWVIEW);
  end

  AVWVIEW = set_coordinates(AVWVIEW);
  %AVWVIEW = set_display_values(AVWVIEW);
  %AVWVIEW = set_crosshairs(AVWVIEW);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case {'left','right','up','down'},

  AVWVIEW = get_slider_position(AVWVIEW);

  % what axes are we in?
  if isequal(gca, AVWVIEW.handles.axial_axes),
    switch command,
     case 'left',
      % decrease sagittal slice
      if AVWVIEW.slices.sag > 1,
        AVWVIEW.slices.sag = AVWVIEW.slices.sag - 1;
      end
     case 'right',
      % increase sagittal slice
      if AVWVIEW.slices.sag < AVWVIEW.xdim,
        AVWVIEW.slices.sag = AVWVIEW.slices.sag + 1;
      end
     case 'up',
      % increase coronal slice
      if AVWVIEW.slices.cor < AVWVIEW.ydim,
        AVWVIEW.slices.cor = AVWVIEW.slices.cor + 1;
      end
     case 'down',
      % decrease coronal slice
      if AVWVIEW.slices.cor > 1,
        AVWVIEW.slices.cor = AVWVIEW.slices.cor - 1;
      end
    end
    switch command,
     case {'left','right'}
      set(AVWVIEW.handles.sagittal_slider,'Value',AVWVIEW.slices.sag);
      Ssag = squeeze(AVWVIEW.avw.img(AVWVIEW.slices.sag,:,:));
      set(AVWVIEW.handles.sagittal_image,'CData',Ssag');
      set(AVWVIEW.handles.sagittal_sliderN,'String',num2str(AVWVIEW.slices.sag));
      set(AVWVIEW.handles.sagittal_sliderN,'Value',AVWVIEW.slices.sag);
     case {'up','down'},
      set(AVWVIEW.handles.coronal_slider,'Value',AVWVIEW.slices.cor);
      Scor = squeeze(AVWVIEW.avw.img(:,AVWVIEW.slices.cor,:));
      set(AVWVIEW.handles.coronal_image,'CData',Scor');
      set(AVWVIEW.handles.coronal_sliderN,'String',num2str(AVWVIEW.slices.cor));
      set(AVWVIEW.handles.coronal_sliderN,'Value',AVWVIEW.slices.cor);
    end
  end

  if isequal(gca, AVWVIEW.handles.coronal_axes),
    switch command,
     case 'left',
      % decrease sagittal slice
      if AVWVIEW.slices.sag > 1,
        AVWVIEW.slices.sag = AVWVIEW.slices.sag - 1;
      end
     case 'right',
      % increase sagittal slice
      if AVWVIEW.slices.sag < AVWVIEW.xdim,
        AVWVIEW.slices.sag = AVWVIEW.slices.sag + 1;
      end
     case 'up',
      % increase axial slice
      if AVWVIEW.slices.axi < AVWVIEW.zdim,
        AVWVIEW.slices.axi = AVWVIEW.slices.axi + 1;
      end
     case 'down',
      % decrease axial slice
      if AVWVIEW.slices.axi > 1,
        AVWVIEW.slices.axi = AVWVIEW.slices.axi - 1;
      end
    end
    switch command,
     case {'left','right'}
      set(AVWVIEW.handles.sagittal_slider,'Value',AVWVIEW.slices.sag);
      Ssag = squeeze(AVWVIEW.avw.img(AVWVIEW.slices.sag,:,:));
      set(AVWVIEW.handles.sagittal_image,'CData',Ssag');
      set(AVWVIEW.handles.sagittal_sliderN,'String',num2str(AVWVIEW.slices.sag));
      set(AVWVIEW.handles.sagittal_sliderN,'Value',AVWVIEW.slices.sag);
     case {'up','down'},
      set(AVWVIEW.handles.axial_slider,'Value',AVWVIEW.slices.axi);
      Saxi = squeeze(AVWVIEW.avw.img(:,:,AVWVIEW.slices.axi));
      set(AVWVIEW.handles.axial_image,'CData',Saxi');
      set(AVWVIEW.handles.axial_sliderN,'String',num2str(AVWVIEW.slices.axi));
      set(AVWVIEW.handles.axial_sliderN,'Value',AVWVIEW.slices.axi);
    end
  end

  if isequal(gca, AVWVIEW.handles.sagittal_axes),
    switch command,
     case 'left',
      % decrease sagittal slice
      if AVWVIEW.slices.cor > 1,
        AVWVIEW.slices.cor = AVWVIEW.slices.cor - 1;
      end
     case 'right',
      % increase sagittal slice
      if AVWVIEW.slices.cor < AVWVIEW.ydim,
        AVWVIEW.slices.cor = AVWVIEW.slices.cor + 1;
      end
     case 'up',
      % increase axial slice
      if AVWVIEW.slices.axi < AVWVIEW.zdim,
        AVWVIEW.slices.axi = AVWVIEW.slices.axi + 1;
      end
     case 'down',
      % decrease axial slice
      if AVWVIEW.slices.axi > 1,
        AVWVIEW.slices.axi = AVWVIEW.slices.axi - 1;
      end
    end
    switch command,
     case {'left','right'}
      set(AVWVIEW.handles.coronal_slider,'Value',AVWVIEW.slices.cor);
      Scor = squeeze(AVWVIEW.avw.img(:,AVWVIEW.slices.cor,:));
      set(AVWVIEW.handles.coronal_image,'CData',Scor');
      set(AVWVIEW.handles.coronal_sliderN,'String',num2str(AVWVIEW.slices.cor));
      set(AVWVIEW.handles.coronal_sliderN,'Value',AVWVIEW.slices.cor);
     case {'up','down'},
      set(AVWVIEW.handles.axial_slider,'Value',AVWVIEW.slices.axi);
      Saxi = squeeze(AVWVIEW.avw.img(:,:,AVWVIEW.slices.axi));
      set(AVWVIEW.handles.axial_image,'CData',Saxi');
      set(AVWVIEW.handles.axial_sliderN,'String',num2str(AVWVIEW.slices.axi));
      set(AVWVIEW.handles.axial_sliderN,'Value',AVWVIEW.slices.axi);
    end
  end

  AVWVIEW = set_coordinates(AVWVIEW);
  %AVWVIEW = set_crosshairs(AVWVIEW);









  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case {'roi_9','roi_7','roi_5','roi_3'},

  position = [ AVWVIEW.slices.sag, AVWVIEW.slices.cor, AVWVIEW.slices.axi ];

  shape.type = 'block';
  if findstr(command,'9'), shape.size = [9,9,9]; end
  if findstr(command,'7'), shape.size = [7,7,7]; end
  if findstr(command,'5'), shape.size = [5,5,5]; end
  if findstr(command,'3'), shape.size = [3,3,3]; end

  stats.roi = avw_roi(AVWVIEW.avw,position,shape);

  AVWVIEW.avw.stats = avw_stats(stats);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case {'nasion','rpa','lpa','ac'},

  % return voxel coordinates into AVWVIEW.imgXYZ
  AVWVIEW = slices2metric(AVWVIEW);

  if get(AVWVIEW.handles.coord,'value') == 2,
    format = '%7.3f %7.3f %7.3f';
    imgXYZ = AVWVIEW.imgXYZ.mm;
    string = sprintf('%7.3f %7.3f %7.3f',imgXYZ);
  elseif get(AVWVIEW.handles.coord,'value') == 3,
    format = '%7.3f %7.3f %7.3f';
    imgXYZ = AVWVIEW.imgXYZ.meters;
    string = sprintf('%7.3f %7.3f %7.3f',imgXYZ);
  else
    imgXYZ = AVWVIEW.imgXYZ.voxels;
    string = sprintf('%7.0f %7.0f %7.0f',imgXYZ);
  end;

  switch command,
   case 'nasion',
    set(AVWVIEW.handles.nasion,'String',string);
    AVWVIEW.avw.fiducials.nasion.voxels(1,:) = AVWVIEW.imgXYZ.voxels;
    AVWVIEW.avw.fiducials.nasion.meters(1,:) = AVWVIEW.imgXYZ.meters;
    AVWVIEW.avw.fiducials.nasion.mm(1,:)     = AVWVIEW.imgXYZ.mm;
   case 'rpa',
    set(AVWVIEW.handles.rpa,'String',string);
    AVWVIEW.avw.fiducials.rpa.voxels(1,:) = AVWVIEW.imgXYZ.voxels;
    AVWVIEW.avw.fiducials.rpa.meters(1,:) = AVWVIEW.imgXYZ.meters;
    AVWVIEW.avw.fiducials.rpa.mm(1,:)     = AVWVIEW.imgXYZ.mm;
   case 'lpa',
    set(AVWVIEW.handles.lpa,'String',string);
    AVWVIEW.avw.fiducials.lpa.voxels(1,:) = AVWVIEW.imgXYZ.voxels;
    AVWVIEW.avw.fiducials.lpa.meters(1,:) = AVWVIEW.imgXYZ.meters;
    AVWVIEW.avw.fiducials.lpa.mm(1,:)     = AVWVIEW.imgXYZ.mm;
   case 'ac',
    set(AVWVIEW.handles.ac,'String',string);
    AVWVIEW.avw.ac.voxels(1,:) = AVWVIEW.imgXYZ.voxels;
    AVWVIEW.avw.ac.meters(1,:) = AVWVIEW.imgXYZ.meters;
    AVWVIEW.avw.ac.mm(1,:) = AVWVIEW.imgXYZ.mm;
  end







  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'flip',

  % flip X dim here...
  AVWVIEW.avw.img = flipdim(AVWVIEW.avw.img,1);

  if isfield(AVWVIEW.handles,'axial_image'),
    Saxi = squeeze(AVWVIEW.avw.img(:,:,AVWVIEW.slices.axi));
    set(AVWVIEW.handles.axial_image,'CData',Saxi');
  end;
  if isfield(AVWVIEW.handles,'coronal_image'),
    Scor = squeeze(AVWVIEW.avw.img(:,AVWVIEW.slices.cor,:));
    set(AVWVIEW.handles.coronal_image,'CData',Scor');
  end;
  if isfield(AVWVIEW.handles,'sagittal_image'),
    Ssag = squeeze(AVWVIEW.avw.img(AVWVIEW.slices.sag,:,:));
    set(AVWVIEW.handles.sagittal_image,'CData',Ssag');
  end;

  flipStatus = get(AVWVIEW.handles.flipStatus,'string');
  if strmatch(flipStatus,'R>>L (radiological)'),
    flipStatus = 'L>>R (neurological)';
  else,
    flipStatus = 'R>>L (radiological)';
  end;
  set(AVWVIEW.handles.flipStatus,'string',flipStatus);

  AVWVIEW = set_coordinates(AVWVIEW);




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case {'brighter','dimmer','setclimit'},

  switch command,
   case 'brighter',
    AVWVIEW.clim = AVWVIEW.clim .* 0.9;
   case 'dimmer',
    AVWVIEW.clim = AVWVIEW.clim .* 1.1;
   case 'setclimit',
    clim = get(AVWVIEW.handles.clim,'string');
    AVWVIEW.clim(2) = str2num(clim);
  end
  climString = sprintf('%05.2f',AVWVIEW.clim(2));
  set(AVWVIEW.handles.clim,'string',climString);

  if isfield(AVWVIEW.handles,'axial_image'),
    set(AVWVIEW.handles.axial_axes,'Clim',AVWVIEW.clim);
  end;
  if isfield(AVWVIEW.handles,'coronal_image'),
    set(AVWVIEW.handles.coronal_axes,'Clim',AVWVIEW.clim);
  end;
  if isfield(AVWVIEW.handles,'sagittal_image'),
    set(AVWVIEW.handles.sagittal_axes,'Clim',AVWVIEW.clim);
  end;




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'contrast'

  if isfield(AVWVIEW.handles,'axial_image'),
    X = squeeze(AVWVIEW.avw.img(:,:,AVWVIEW.slices.axi));
    AVWVIEW.cmap = contrast(X);
    colormap(AVWVIEW.handles.axial_axes,AVWVIEW.cmap);
  end;
  if isfield(AVWVIEW.handles,'coronal_image'),
    X = squeeze(AVWVIEW.avw.img(:,AVWVIEW.slices.cor,:));
    AVWVIEW.cmap = contrast(X);
    colormap(AVWVIEW.handles.axial_axes,AVWVIEW.cmap);
  end;
  if isfield(AVWVIEW.handles,'sagittal_image'),
    X = squeeze(AVWVIEW.avw.img(AVWVIEW.slices.sag,:,:));
    AVWVIEW.cmap = contrast(X);
    colormap(AVWVIEW.handles.axial_axes,AVWVIEW.cmap);
  end;




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'setcmap'

  cmapIndex = get(AVWVIEW.handles.cmap,'value');
  cmapString = get(AVWVIEW.handles.cmap,'string');

  AVWVIEW.cmapString = cmapString{cmapIndex};

  if isfield(AVWVIEW.handles,'axial_image'),
    AVWVIEW.cmap = colormap(AVWVIEW.handles.axial_axes,AVWVIEW.cmapString);
  end;
  if isfield(AVWVIEW.handles,'coronal_image'),
    AVWVIEW.cmap = colormap(AVWVIEW.handles.coronal_axes,AVWVIEW.cmapString);
  end;
  if isfield(AVWVIEW.handles,'sagittal_image'),
    AVWVIEW.cmap = colormap(AVWVIEW.handles.sagittal_axes,AVWVIEW.cmapString);
  end;




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'histogram',

  avw = AVWVIEW.avw;

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

  [bins,intensity_volume] = avw_histogram(avw);

  figure('name','intensity histogram');
  bar(bins,intensity_volume);




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'zoom'

  x = AVWVIEW.slices.sag;
  y = AVWVIEW.slices.cor;
  z = AVWVIEW.slices.axi;

  switch AVWVIEW.view,
   case 'axi',
    slice = rot90(AVWVIEW.avw.img(:,:,z));
   case 'cor',
    slice = rot90(squeeze(AVWVIEW.avw.img(:,y,:)));
   case 'sag',
    slice = rot90(squeeze(AVWVIEW.avw.img(x,:,:)));
  end

  figure;
  imagesc(slice,AVWVIEW.clim);
  colormap(gray); axis off; zoom on;
  daspect(AVWVIEW.daspect);




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'save_image',

  [filename, pathname] = uiputfile(...
      { '*.png','PNG Files (*.png)'; ...
        '*.jpg','JPG Files (*.jpg)'; ...
        '*.ppm;*.pgm;*.pbm','Portable Anymap (*.ppm,*.pgm,*.pbm)'; ...
        '*.*',  'All Files (*.*)'}, ...
      'IMWrite to file');

  if filename,

    x = AVWVIEW.slices.sag;
    y = AVWVIEW.slices.cor;
    z = AVWVIEW.slices.axi;

    pixelsPerMM = 1 ./ double(AVWVIEW.avw.hdr.dime.pixdim(2:4));
    pixelsPerMeter = pixelsPerMM .* 1000;

    switch AVWVIEW.view,
     case 'axi',
      slice = rot90(AVWVIEW.avw.img(:,:,z));
      xresolution = pixelsPerMeter(1);
      yresolution = pixelsPerMeter(2);
     case 'cor',
      slice = rot90(squeeze(AVWVIEW.avw.img(:,y,:)));
      xresolution = pixelsPerMeter(1);
      yresolution = pixelsPerMeter(3);
     case 'sag',
      slice = rot90(squeeze(AVWVIEW.avw.img(x,:,:)));
      xresolution = pixelsPerMeter(2);
      yresolution = pixelsPerMeter(3);
    end

    % scale the image values to between 0-1 for imwrite
    % (initially used max of slice, but actually want to
    % use the scaled intensity from the gui).
    %maxValue = max(max(slice));
    %scaledSlice = slice ./ maxValue;
    scaledSlice = slice ./ AVWVIEW.clim(2);

    % RGB = ind2rgb(X,map) converts the matrix X and corresponding
    % colormap map to RGB (truecolor)

    file = [pathname,filename];
    fprintf('saving to:...%s\n',file);

    %Most of the supported image file formats store uint8 data.
    %PNG and TIFF formats additionally support uint16 data. For
    %grayscale and RGB images, if the data array is double, the
    %assumed dynamic range is [0,1]. The data array is automatically
    %scaled by 255 before being written as uint8. If the data array
    %is uint8 or uint16, it is written without scaling as uint8 or
    %uint16, respectively.

    [pathname,filename,ext] = fileparts([pathname,filename]);

    switch ext,

     case '.png',
      %imwrite(uint8(image),colormap,file,format);
      %imwrite(slice,colormap(gray),file,format,...

      %     'XResolution'  A scalar indicating the number of
      %                    pixels/unit in the horizontal direction
      %
      %     'YResolution'  A scalar indicating the number of
      %                    pixels/unit in the vertical direction
      %
      %     'ResolutionUnit' Either 'unknown' or 'meter'
      format = 'png';
      imwrite(scaledSlice,file,format,...
              'BitDepth',16,'ResolutionUnit','meter',...
              'XResolution',xresolution,...
              'YResolution',yresolution);

     case '.jpg',
      format = 'jpg';
      imwrite(scaledSlice,file,format);
     case {'.ppm','.pgm','.pbm'},
      format = ext(2:end);
      imwrite(scaledSlice,file,format);
     otherwise
      fprintf('...cannot write %s image files\n',ext);
    end
  end





  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 case 'quit',

  flipStatus = get(AVWVIEW.handles.flipStatus,'string');
  if strmatch(flipStatus,'L>>R (neurological)'),
    fprintf('...returning avw.img to radiological orientation;\n');
    fprintf('...use avw_flip to permanently change the orientation.\n');
    AVWVIEW.avw.img = flipdim(AVWVIEW.avw.img,1);
  end;

  % try to assign back into the input struct
  if ~isempty(AVWVIEW.invarname),
    fprintf('...returning data to base workspace struct ''%s''\n',AVWVIEW.invarname);
    assignin('base',AVWVIEW.invarname,AVWVIEW.avw);
  elseif evalin('base','exist(''mri'',''var'')'),
    fprintf('...returning data to base workspace struct ''mri''\n');
    string = ['AVWVIEW = get(', num2str(AVWVIEW.gui),...
              ',''Userdata''); mri.data = AVWVIEW.avw; clear AVWVIEW;'];
    evalin('base',string);
  else
    fprintf('...returning data into base workspace struct ''avw''\n');
    assignin('base','avw',AVWVIEW.avw);
  end
  close gcbf;

 otherwise,

end


switch command,
 case 'quit',
 otherwise,
  set(AVWVIEW.gui,'UserData',AVWVIEW);
end

if nargout > 0,
  varargout{1} = AVWVIEW.avw;
end


return









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function axial_update(AVWVIEW)

if isfield(AVWVIEW.handles,'axial_image'),
  Saxi = squeeze(AVWVIEW.avw.img(:,:,AVWVIEW.slices.axi));
  set(AVWVIEW.handles.axial_image,'CData',Saxi');
end
if isfield(AVWVIEW.handles,'axial_sliderN'),
  set(AVWVIEW.handles.axial_sliderN,'String',num2str(AVWVIEW.slices.axi));
  set(AVWVIEW.handles.axial_sliderN,'Value',AVWVIEW.slices.axi);
end
if isfield(AVWVIEW.handles,'axial_slider'),
  set(AVWVIEW.handles.axial_slider,'Value',AVWVIEW.slices.axi);
end;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coronal_update(AVWVIEW)

if isfield(AVWVIEW.handles,'coronal_image'),
  Scor = squeeze(AVWVIEW.avw.img(:,AVWVIEW.slices.cor,:));
  set(AVWVIEW.handles.coronal_image,'CData',Scor');
end
if isfield(AVWVIEW.handles,'coronal_sliderN'),
  set(AVWVIEW.handles.coronal_sliderN,'String',num2str(AVWVIEW.slices.cor));
  set(AVWVIEW.handles.coronal_sliderN,'Value',AVWVIEW.slices.cor);
end
if isfield(AVWVIEW.handles,'coronal_slider'),
  set(AVWVIEW.handles.coronal_slider,'Value',AVWVIEW.slices.cor);
end;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sagittal_update(AVWVIEW)

if isfield(AVWVIEW.handles,'sagittal_image'),
  Ssag = squeeze(AVWVIEW.avw.img(AVWVIEW.slices.sag,:,:));
  set(AVWVIEW.handles.sagittal_image,'CData',Ssag');
end
if isfield(AVWVIEW.handles,'sagittal_sliderN'),
  set(AVWVIEW.handles.sagittal_sliderN,'String',num2str(AVWVIEW.slices.sag));
  set(AVWVIEW.handles.sagittal_sliderN,'Value',AVWVIEW.slices.sag);
end
if isfield(AVWVIEW.handles,'sagittal_slider'),
  set(AVWVIEW.handles.sagittal_slider,'Value',AVWVIEW.slices.sag);
end;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = set_crosshairs(AVWVIEW)

current_axes = gca;

[AVWVIEW, metric] = slices2metric(AVWVIEW);

if isfield(AVWVIEW.handles,'axial_axes'),

  axes(AVWVIEW.handles.axial_axes)

  axialXline = findobj('Tag','axialXline');
  if axialXline,
    delete(axialXline);
  end
  axialYline = findobj('Tag','axialYline');
  if axialYline,
    delete(axialYline);
  end
  
  AVWVIEW.handles.axial_xlim = get(AVWVIEW.handles.axial_axes,'Xlim');
  AVWVIEW.handles.axial_ylim = get(AVWVIEW.handles.axial_axes,'Ylim');
  
  AVWVIEW.handles.axial_xline = line('Xdata',[metric.sag metric.sag],'Ydata',AVWVIEW.handles.axial_ylim);
  AVWVIEW.handles.axial_yline = line('Ydata',[metric.cor metric.cor],'Xdata',AVWVIEW.handles.axial_xlim);
  
  set(AVWVIEW.handles.axial_xline,'Color','b','EraseMode','xor','Tag','axialXline');
  set(AVWVIEW.handles.axial_yline,'Color','b','EraseMode','xor','Tag','axialYline');

  if get(AVWVIEW.handles.crosshairs,'value'),
    set(AVWVIEW.handles.axial_xline,'visible','on');
    set(AVWVIEW.handles.axial_yline,'visible','on');
  else
    set(AVWVIEW.handles.axial_xline,'visible','off');
    set(AVWVIEW.handles.axial_yline,'visible','off');
  end
end

if isfield(AVWVIEW.handles,'coronal_axes'),

  axes(AVWVIEW.handles.coronal_axes);
  
  coronalXline = findobj('Tag','coronalXline');
  if coronalXline,
    delete(coronalXline);
  end
  coronalYline = findobj('Tag','coronalYline');
  if coronalYline,
    delete(coronalYline);
  end
  
  AVWVIEW.handles.coronal_xlim = get(AVWVIEW.handles.coronal_axes,'Xlim');
  AVWVIEW.handles.coronal_ylim = get(AVWVIEW.handles.coronal_axes,'Ylim');
  AVWVIEW.handles.coronal_xline = line('Xdata',[metric.sag metric.sag],'Ydata',AVWVIEW.handles.coronal_ylim);
  AVWVIEW.handles.coronal_yline = line('Ydata',[metric.axi metric.axi],'Xdata',AVWVIEW.handles.coronal_xlim);
  set(AVWVIEW.handles.coronal_xline,'Color','b','EraseMode','xor','Tag','coronalXline');
  set(AVWVIEW.handles.coronal_yline,'Color','b','EraseMode','xor','Tag','coronalYline');

  if get(AVWVIEW.handles.crosshairs,'value'),
    set(AVWVIEW.handles.coronal_xline,'visible','on');
    set(AVWVIEW.handles.coronal_yline,'visible','on');
  else
    set(AVWVIEW.handles.coronal_xline,'visible','off');
    set(AVWVIEW.handles.coronal_yline,'visible','off');
  end
end

if isfield(AVWVIEW.handles,'sagittal_axes'),
  
  axes(AVWVIEW.handles.sagittal_axes);
  
  sagittalXline = findobj('Tag','sagittalXline');
  if sagittalXline,
    delete(sagittalXline);
  end
  sagittalYline = findobj('Tag','sagittalYline');
  if sagittalYline,
    delete(sagittalYline);
  end
  
  AVWVIEW.handles.sagittal_xlim = get(AVWVIEW.handles.sagittal_axes,'Xlim');
  AVWVIEW.handles.sagittal_ylim = get(AVWVIEW.handles.sagittal_axes,'Ylim');
  AVWVIEW.handles.sagittal_xline = line('Xdata',[metric.cor metric.cor],'Ydata',AVWVIEW.handles.sagittal_ylim);
  AVWVIEW.handles.sagittal_yline = line('Ydata',[metric.axi metric.axi],'Xdata',AVWVIEW.handles.sagittal_xlim);
  set(AVWVIEW.handles.sagittal_xline,'Color','b','EraseMode','xor','Tag','sagittalXline');
  set(AVWVIEW.handles.sagittal_yline,'Color','b','EraseMode','xor','Tag','sagittalYline');

  if get(AVWVIEW.handles.crosshairs,'value'),
    set(AVWVIEW.handles.sagittal_xline,'visible','on');
    set(AVWVIEW.handles.sagittal_yline,'visible','on');
  else
    set(AVWVIEW.handles.sagittal_xline,'visible','off');
    set(AVWVIEW.handles.sagittal_yline,'visible','off');
  end

end

axes(current_axes);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = set_coordinates(AVWVIEW)

% set the axis coordinates to voxels, mm or meters

s = size(AVWVIEW.avw.img);
if length(s) > 0, xdim = s(1); else xdim = 1; end
if length(s) > 1, ydim = s(2); else ydim = 1; end
if length(s) > 2, zdim = s(3); else zdim = 1; end

% initialise for voxel coordinates
xpixdim = double(AVWVIEW.avw.hdr.dime.pixdim(2));
ypixdim = double(AVWVIEW.avw.hdr.dime.pixdim(3));
zpixdim = double(AVWVIEW.avw.hdr.dime.pixdim(4));
xdata = [0 xdim];
ydata = [0 ydim];
zdata = [0 zdim];

aspect = 1./AVWVIEW.daspect; %AVWVIEW.avw.hdr.dime.pixdim(2:4);


if get(AVWVIEW.handles.coord,'value') == 2,    % mm
  xdata = xdata .* xpixdim;
  ydata = ydata .* ypixdim;
  zdata = zdata .* zpixdim;
  aspect = [1 1 1];
end

if get(AVWVIEW.handles.coord,'value') == 3,    % meters
  xpixdim = xpixdim / 1000;
  ypixdim = ypixdim / 1000;
  zpixdim = zpixdim / 1000;
  xdata = xdata .* xpixdim;
  ydata = ydata .* ypixdim;
  zdata = zdata .* zpixdim;
  aspect = [1 1 1];
end


if isfield(AVWVIEW.handles,'axial_image'),
  set(AVWVIEW.handles.axial_axes,'Xlim',xdata);
  set(AVWVIEW.handles.axial_axes,'Ylim',ydata);
  set(AVWVIEW.handles.axial_axes,'Zlim',zdata);
  set(AVWVIEW.handles.axial_image,'Xdata',xdata);
  set(AVWVIEW.handles.axial_image,'Ydata',ydata);
  daspect(AVWVIEW.handles.axial_axes,aspect([1 2 3]));
end;
if isfield(AVWVIEW.handles,'coronal_image'),
  set(AVWVIEW.handles.coronal_axes,'Xlim',xdata);
  set(AVWVIEW.handles.coronal_axes,'Ylim',zdata);
  set(AVWVIEW.handles.coronal_axes,'Zlim',ydata);
  set(AVWVIEW.handles.coronal_image,'Xdata',xdata);
  set(AVWVIEW.handles.coronal_image,'Ydata',zdata);
  daspect(AVWVIEW.handles.coronal_axes,aspect([1 3 2]));
end;
if isfield(AVWVIEW.handles,'sagittal_image'),
  set(AVWVIEW.handles.sagittal_axes,'Xlim',ydata);
  set(AVWVIEW.handles.sagittal_axes,'Ylim',zdata);
  set(AVWVIEW.handles.sagittal_axes,'Zlim',xdata);
  set(AVWVIEW.handles.sagittal_image,'Xdata',ydata);
  set(AVWVIEW.handles.sagittal_image,'Ydata',zdata);
  daspect(AVWVIEW.handles.sagittal_axes,aspect([2 3 1]));
end;

AVWVIEW = set_display_values(AVWVIEW);
AVWVIEW = set_crosshairs(AVWVIEW);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = get_current_position(AVWVIEW),

AVWVIEW = get_slider_position(AVWVIEW);
[AVWVIEW, metric] = slices2metric(AVWVIEW);

switch AVWVIEW.view,
 case 'sag',
  currentpoint = get(get(AVWVIEW.handles.sagittal_image,'Parent'),'CurrentPoint');
  metric.cor = currentpoint(1,1);
  metric.axi = currentpoint(1,2);
 case 'cor',
  currentpoint = get(get(AVWVIEW.handles.coronal_image,'Parent'),'CurrentPoint');
  metric.sag = currentpoint(2,1);
  metric.axi = currentpoint(2,2);
 case 'axi',
  currentpoint = get(get(AVWVIEW.handles.axial_image,'Parent'),'CurrentPoint');
  metric.sag = currentpoint(2,1);
  metric.cor = currentpoint(2,2);
end

AVWVIEW = metric2slices(AVWVIEW,metric);
AVWVIEW = check_slices(AVWVIEW);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = get_slider_position(AVWVIEW),

[AVWVIEW.slices.sag,AVWVIEW.slices.cor,AVWVIEW.slices.axi] = deal(0);

if isfield(AVWVIEW.handles,'sagittal_slider'),
  if ishandle(AVWVIEW.handles.sagittal_slider),
    AVWVIEW.slices.sag = round(get(AVWVIEW.handles.sagittal_slider,'Value'));
  end
end
if AVWVIEW.slices.sag == 0,
  if isfield(AVWVIEW.handles,'sagittal_sliderN'),
    if ishandle(AVWVIEW.handles.sagittal_sliderN),
      AVWVIEW.slices.sag = round(get(AVWVIEW.handles.sagittal_sliderN,'Value'));
    end
  end
end

if isfield(AVWVIEW.handles,'coronal_slider'),
  if ishandle(AVWVIEW.handles.coronal_slider),
    AVWVIEW.slices.cor = round(get(AVWVIEW.handles.coronal_slider,'Value'));
  end
end
if AVWVIEW.slices.cor == 0,
  if isfield(AVWVIEW.handles,'coronal_sliderN'),
    if ishandle(AVWVIEW.handles.coronal_sliderN),
      AVWVIEW.slices.cor = round(get(AVWVIEW.handles.coronal_sliderN,'Value'));
    end
  end
end

if isfield(AVWVIEW.handles,'axial_slider'),
  if ishandle(AVWVIEW.handles.axial_slider),
    AVWVIEW.slices.axi = round(get(AVWVIEW.handles.axial_slider,'Value'));
  end
end
if AVWVIEW.slices.axi == 0,
  if isfield(AVWVIEW.handles,'axial_sliderN'),
    if ishandle(AVWVIEW.handles.axial_sliderN),
      AVWVIEW.slices.axi = round(get(AVWVIEW.handles.axial_sliderN,'Value'));
    end
  end
end

AVWVIEW = check_slices(AVWVIEW);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = check_slices(AVWVIEW),

adjust = 0;

[ SagSize, CorSize, AxiSize ] = size(AVWVIEW.avw.img);

if AVWVIEW.slices.sag > SagSize,
  AVWVIEW.slices.sag = SagSize;
  adjust = 1;
end;
if AVWVIEW.slices.sag < 1,
  AVWVIEW.slices.sag = 1;
  adjust = 1;
end;
if AVWVIEW.slices.cor > CorSize,
  AVWVIEW.slices.cor = CorSize;
  adjust = 1;
end;
if AVWVIEW.slices.cor < 1,
  AVWVIEW.slices.cor = 1;
  adjust = 1;
end;
if AVWVIEW.slices.axi > AxiSize,
  AVWVIEW.slices.axi = AxiSize;
  adjust = 1;
end;
if AVWVIEW.slices.axi < 1,
  AVWVIEW.slices.axi = 1;
  adjust = 1;
end;

if adjust,
  AVWVIEW = slices2metric(AVWVIEW);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = set_display_values(AVWVIEW),

% get coordinates of selected voxel and the image intensity there

sag = AVWVIEW.slices.sag;
cor = AVWVIEW.slices.cor;
axi = AVWVIEW.slices.axi;

imgvalue = AVWVIEW.avw.img(sag,cor,axi);

set(AVWVIEW.handles.imval,'String',sprintf('%7.2f',imgvalue));
set(AVWVIEW.handles.imval,'Value',imgvalue);

% Now update the image position text for the selected voxel

[AVWVIEW, metric] = slices2metric(AVWVIEW);
sag = metric.sag;
cor = metric.cor;
axi = metric.axi;

string = sprintf('%7.3f %7.3f %7.3f',sag,cor,axi);

set(AVWVIEW.handles.impos,'String',string);
set(AVWVIEW.handles.impos,'Value',[sag,cor,axi]);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AVWVIEW,metric] = slices2metric(AVWVIEW),

AVWVIEW.imgXYZ.voxels = [AVWVIEW.slices.sag,AVWVIEW.slices.cor,AVWVIEW.slices.axi];
AVWVIEW.imgXYZ.meters = AVWVIEW.imgXYZ.voxels .* AVWVIEW.scale2meters;
AVWVIEW.imgXYZ.mm     = AVWVIEW.imgXYZ.voxels .* AVWVIEW.scale2mm;

coord_value = get(AVWVIEW.handles.coord,'value');

if coord_value == 2,
  % using mm
  img_mm = AVWVIEW.imgXYZ.mm;
  metric.axi = img_mm(3);
  metric.cor = img_mm(2);
  metric.sag = img_mm(1);
elseif coord_value == 3,
  % using meters
  img_meters = AVWVIEW.imgXYZ.meters;
  metric.axi = img_meters(3);
  metric.cor = img_meters(2);
  metric.sag = img_meters(1);
else
  % voxels
  metric = AVWVIEW.slices;
end;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = metric2slices(AVWVIEW,metric),

coord_value = get(AVWVIEW.handles.coord,'value');

if coord_value == 2,
  % using mm
  xpix = double(AVWVIEW.avw.hdr.dime.pixdim(2));
  ypix = double(AVWVIEW.avw.hdr.dime.pixdim(3));
  zpix = double(AVWVIEW.avw.hdr.dime.pixdim(4));
  AVWVIEW.slices.axi = round(metric.axi / zpix);
  AVWVIEW.slices.cor = round(metric.cor / ypix);
  AVWVIEW.slices.sag = round(metric.sag / xpix);
elseif coord_value == 3,
  % using meters
  xpix = double(AVWVIEW.avw.hdr.dime.pixdim(2)) / 1000;
  ypix = double(AVWVIEW.avw.hdr.dime.pixdim(3)) / 1000;
  zpix = double(AVWVIEW.avw.hdr.dime.pixdim(4)) / 1000;
  AVWVIEW.slices.axi = round(metric.axi / zpix);
  AVWVIEW.slices.cor = round(metric.cor / ypix);
  AVWVIEW.slices.sag = round(metric.sag / xpix);
else
  % voxels
  AVWVIEW.slices.axi = round(metric.axi);
  AVWVIEW.slices.cor = round(metric.cor);
  AVWVIEW.slices.sag = round(metric.sag);
end;

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AVWVIEW = init(avw,invarname),

% try to keep track of the input struct
if exist('invarname','var'),
  if ~isempty(invarname),
    AVWVIEW.invarname = invarname;
  else
    AVWVIEW.invarname = '';
  end
else
  AVWVIEW.invarname = '';
end

% GUI General Parameters
GUIwidth  = 150;
GUIheight = 50;

version = '[$Revision: 1.5 $]';
name = sprintf('AVW View [v%s]',version(12:16));

if isfield(avw,'fileprefix'),
  if ~isempty(avw.fileprefix),
    format = strcat('%+',sprintf('%d',length(avw.fileprefix)+1),'s');
    name = strcat(name,' - ',sprintf(format,avw.fileprefix));
  end
end

% The Zbuffer provides smooth slice animations (OpenGL is no good)
GUI = figure('Name',name,'Tag','AVWVIEW','units','characters',...
             'BackingStore','off',...
             'NumberTitle','off','color',[0 0 0],...
             'MenuBar','figure','Position',[1 1 GUIwidth GUIheight],...
             'Renderer','zbuffer');

movegui(GUI,'center');

set(GUI,'KeyPressFcn','avw_view([],[],''keypress'');');

AVWVIEW.gui = GUI;

Font.FontName   = 'Helvetica';
Font.FontUnits  = 'Pixels';
Font.FontSize   = 9;
Font.FontWeight = 'normal';
Font.FontAngle  = 'normal';


AVWVIEW.shading = 'flat';
shading(AVWVIEW.shading)



% 	% determine the datatype of avw.img
% 	switch double(avw.hdr.dime.bitpix),
% 	case 1,
%         fprintf('...converting avw.img to uint8 for viewing only.\n\n');
%         avw.img = uint8(avw.img);
% 	case 8,
%         fprintf('...converting avw.img to uint8 for viewing only.\n\n');
%         avw.img = uint8(avw.img);
% 	case 16,
%         fprintf('...converting avw.img to uint16 for viewing only.\n\n');
%         avw.img = uint16(avw.img);
% 	case {32,64},
%         % make sure it is double, not single
%         avw.img = double(avw.img);
% 	otherwise,
%         % do nothing, leave it as is
% 	end


% calculate image stats
intensityMean = mean(mean(mean(avw.img)));
intensityMeanRobust = mean(mean(mean(avw.img(find(avw.img)))));
intensityStdev = std(std(std(avw.img)));

intensityMax = max(max(max(avw.img)));

% try to adjust for large intensity extremes
if intensityMax > 255,
  % got 16 bit data, at least
  adjust = 0.5;
else
  adjust = 0.9;
end
%AVWVIEW.clim = [0 (intensityMeanRobust + (10 * intensityStdev)) ];
AVWVIEW.clim = [0 (intensityMax * adjust)];



AVWVIEW.xdim = size(avw.img,1);
AVWVIEW.ydim = size(avw.img,2);
AVWVIEW.zdim = size(avw.img,3);

AVWVIEW.slices.sag = 1;
AVWVIEW.slices.cor = 1;
AVWVIEW.slices.axi = 1;
if AVWVIEW.xdim > 1, AVWVIEW.slices.sag = floor(AVWVIEW.xdim/2); end
if AVWVIEW.ydim > 1, AVWVIEW.slices.cor = floor(AVWVIEW.ydim/2); end
if AVWVIEW.zdim > 1, AVWVIEW.slices.axi = floor(AVWVIEW.zdim/2); end

% store the volume center for later reference when
% calculating fiducial locations
center = avw_center(avw);
AVWVIEW.center = center.abs.voxels;

% set the default origin at the center
AVWVIEW.origin  = AVWVIEW.center;

AVWVIEW.pixdim = double(avw.hdr.dime.pixdim(2:4));

AVWVIEW.scale2mm     = AVWVIEW.pixdim;          % vol scale in mm
AVWVIEW.scale2meters = AVWVIEW.pixdim ./ 1000;  % vol scale in meters
AVWVIEW.daspect      = AVWVIEW.pixdim ./ min(AVWVIEW.pixdim);

xPlotSize = 0.38;
yPlotSize = 0.38;

AVWVIEW.cmapString = 'gray';

AVWVIEW.cmap = colormap(AVWVIEW.cmapString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Axial Slice
if AVWVIEW.xdim > 1 & AVWVIEW.ydim > 1,

  handles.axial_subplot = subplot('position',[0.075 0.075 xPlotSize yPlotSize]);

  Saxial = squeeze(avw.img(:,:,AVWVIEW.slices.axi));

  %handles.axial_image = pcolor(double(Saxial)); colormap(gray); shading interp
  %pcolor(Saxial'); colormap(gray); shading flat

  %surf(avw.img(:,:,20)','edgecolor','none'); view(2); axis tight; colormap(flipdim(gray,1)); colorbar
  %surf(avw.img(:,:,20)','edgecolor','none'); view(2); axis tight; colormap(gray)
  %image(avw.img(:,:,20)'); axis image; colormap(gray)
  %imagesc(avw.img(:,:,20)'); axis image; colormap(gray)

  %handles.axial_image = image('Cdata',Saxial','CDataMapping','scaled',...
  %    'XData',[0 xdim],'YData',[0 ydim]);

  handles.axial_image = imagesc([0,AVWVIEW.xdim],[0,AVWVIEW.ydim],Saxial',AVWVIEW.clim);

  handles.axial_axes = gca;
  set(gca,'YDir','normal','XLimMode','manual','YLimMode','manual',...
          'ClimMode','manual','YColor',[1 1 1],'XColor',[1 1 1])

  title('Axial','Color',[1 1 1])
  ylabel('Y')
  xlabel('X')
  %xlabel('(Right <<) X (>> Left)'); % default radiological orientation for Analyze
  %xlabel('(Left <<) X (>> Right)')

  % This callback navigates with mouse button click
  set(handles.axial_image,'ButtonDownFcn','avw_view([],[],''axial_image'');');

  GUIheight = 0.46;

  if AVWVIEW.zdim > 1,
    slider_step(1) = 1/(AVWVIEW.zdim);
    slider_step(2) = 1/(AVWVIEW.zdim);
    handles.axial_slider = uicontrol('Parent',GUI,'Style','slider',...
                                     'Units','Normalized', Font, ...
                                     'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'center',...
                                     'BusyAction','queue',...
                                     'TooltipString','Axial slice navigation',...
                                     'Min',1,'Max',AVWVIEW.zdim,'SliderStep',slider_step,'Value',AVWVIEW.slices.axi,...
                                     'Callback','avw_view([],[],''axial_slider'');');
  end
  handles.axial_sliderN = uicontrol('Parent',GUI,'Style','text',...
                                    'Units','Normalized', Font, ...
                                    'Position',[.65 GUIheight .03 .03], 'HorizontalAlignment', 'center',...
                                    'BackgroundColor',[0 0 0],'ForegroundColor',[1 1 1],...
                                    'BusyAction','queue',...
                                    'TooltipString','Axial slice number',...
                                    'String',num2str(AVWVIEW.slices.axi),'Value',AVWVIEW.slices.axi);
  handles.axial_sliderT = uicontrol('Parent',GUI,'Style','text',...
                                    'Units','Normalized', Font, ...
                                    'Position',[.70 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                                    'BackgroundColor',[0 0 0],'ForegroundColor',[1 1 1],...
                                    'BusyAction','queue',...
                                    'TooltipString','Axial slice navigation',...
                                    'String','Axial');

  handles.axial_xlim = get(handles.axial_axes,'Xlim');
  handles.axial_ylim = get(handles.axial_axes,'Ylim');
  handles.axial_xline = line('Xdata',[AVWVIEW.slices.sag AVWVIEW.slices.sag],'Ydata',handles.axial_ylim);
  handles.axial_yline = line('Ydata',[AVWVIEW.slices.cor AVWVIEW.slices.cor],'Xdata',handles.axial_xlim);
  set(handles.axial_xline,'Color','b','EraseMode','xor','Tag','axialXline');
  set(handles.axial_yline,'Color','b','EraseMode','xor','Tag','axialYline');
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coronal Slice
if AVWVIEW.xdim > 1 & AVWVIEW.zdim > 1,

  handles.coronal_subplot = subplot('position',[0.075 0.575 xPlotSize yPlotSize]);

  Scor = squeeze(avw.img(:,AVWVIEW.slices.cor,:));
  handles.coronal_image = imagesc([0,AVWVIEW.xdim],[0,AVWVIEW.zdim],Scor',AVWVIEW.clim);

  handles.coronal_axes = gca;
  set(gca,'YDir','normal','XLimMode','manual','YLimMode','manual',...
          'ClimMode','manual','YColor',[1 1 1],'XColor',[1 1 1])

  %xlabel('(Left <<) X (>> Right)')
  %xlabel('(Right <<) X (>> Left)')
  xlabel('X')
  ylabel('Z')
  title('Coronal','Color',[1 1 1])

  % This callback navigates with left click
  set(handles.coronal_image,'ButtonDownFcn','avw_view([],[],''coronal_image'');');

  GUIheight = GUIheight - 0.04;

  if AVWVIEW.ydim > 1,
    slider_step(1) = 1/(AVWVIEW.ydim);
    slider_step(2) = 1/(AVWVIEW.ydim);
    handles.coronal_slider = uicontrol('Parent',GUI,'Style','slider',...
                                       'Units','Normalized', Font, ...
                                       'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'center',...
                                       'BusyAction','queue',...
                                       'TooltipString','Coronal slice navigation',...
                                       'Min',1,'Max',AVWVIEW.ydim,'SliderStep',slider_step,'Value',AVWVIEW.slices.cor,...
                                       'Callback','avw_view([],[],''coronal_slider'');');
  end
  handles.coronal_sliderN = uicontrol('Parent',GUI,'Style','text',...
                                      'Units','Normalized', Font, ...
                                      'Position',[.65 GUIheight .03 .03], 'HorizontalAlignment', 'center',...
                                      'BackgroundColor',[0 0 0],'ForegroundColor',[1 1 1],...
                                      'BusyAction','queue',...
                                      'TooltipString','Coronal slice number',...
                                      'String',num2str(AVWVIEW.slices.cor),'Value',AVWVIEW.slices.cor);
  handles.coronal_sliderT = uicontrol('Parent',GUI,'Style','text',...
                                      'Units','Normalized', Font, ...
                                      'Position',[.70 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                                      'BackgroundColor',[0 0 0],'ForegroundColor',[1 1 1],...
                                      'BusyAction','queue',...
                                      'TooltipString','Coronal slice navigation',...
                                      'String','Coronal');

  handles.coronal_xlim = get(handles.coronal_axes,'Xlim');
  handles.coronal_ylim = get(handles.coronal_axes,'Ylim');
  handles.coronal_xline = line('Xdata',[AVWVIEW.slices.sag AVWVIEW.slices.sag],'Ydata',handles.coronal_ylim);
  handles.coronal_yline = line('Ydata',[AVWVIEW.slices.axi AVWVIEW.slices.axi],'Xdata',handles.coronal_xlim);
  set(handles.coronal_xline,'Color','b','EraseMode','xor','Tag','coronalXline');
  set(handles.coronal_yline,'Color','b','EraseMode','xor','Tag','coronalYline');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sagittal Slice
if AVWVIEW.ydim > 1 & AVWVIEW.zdim > 1,

  handles.sagittal_subplot = subplot('position',[0.575 0.575 xPlotSize yPlotSize]);

  Ssag = squeeze(avw.img(AVWVIEW.slices.sag,:,:));
  handles.sagittal_image = imagesc([0,AVWVIEW.ydim],[0,AVWVIEW.zdim],Ssag',AVWVIEW.clim);

  handles.sagittal_axes = gca;
  set(gca,'YDir','normal','XLimMode','manual','YLimMode','manual',...
          'ClimMode','manual','YColor',[1 1 1],'XColor',[1 1 1])

  xlabel('Y')
  ylabel('Z')
  title('Sagittal','Color',[1 1 1])

  % This callback navigates with mouse click
  set(handles.sagittal_image,'ButtonDownFcn','avw_view([],[],''sagittal_image'');');

  GUIheight = GUIheight - 0.04;

  if AVWVIEW.xdim > 1,
    slider_step(1) = 1/(AVWVIEW.xdim);
    slider_step(2) = 1/(AVWVIEW.xdim);
    handles.sagittal_slider = uicontrol('Parent',GUI,'Style','slider',...
                                        'Units','Normalized', Font, ...
                                        'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'center',...
                                        'BusyAction','queue',...
                                        'TooltipString','Sagittal slice navigation',...
                                        'Min',1,'Max',AVWVIEW.xdim,'SliderStep',slider_step,'Value',AVWVIEW.slices.sag,...
                                        'Callback','avw_view([],[],''sagittal_slider'');');
  end
  handles.sagittal_sliderN = uicontrol('Parent',GUI,'Style','text',...
                                       'Units','Normalized', Font, ...
                                       'Position',[.65 GUIheight .03 .03], 'HorizontalAlignment', 'center',...
                                       'BackgroundColor',[0 0 0],'ForegroundColor',[1 1 1],...
                                       'BusyAction','queue',...
                                       'TooltipString','Sagittal slice number',...
                                       'String',num2str(AVWVIEW.slices.sag),'Value',AVWVIEW.slices.sag);
  handles.sagittal_sliderT = uicontrol('Parent',GUI,'Style','text',...
                                       'Units','Normalized', Font, ...
                                       'Position',[.70 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                                       'BackgroundColor',[0 0 0],'ForegroundColor',[1 1 1],...
                                       'BusyAction','queue',...
                                       'TooltipString','Sagittal slice navigation',...
                                       'String','Sagittal');

  handles.sagittal_xlim = get(handles.sagittal_axes,'Xlim');
  handles.sagittal_ylim = get(handles.sagittal_axes,'Ylim');
  handles.sagittal_xline = line('Xdata',[AVWVIEW.slices.cor AVWVIEW.slices.cor],'Ydata',handles.sagittal_ylim);
  handles.sagittal_yline = line('Ydata',[AVWVIEW.slices.axi AVWVIEW.slices.axi],'Xdata',handles.sagittal_xlim);
  set(handles.sagittal_xline,'Color','b','EraseMode','xor','Tag','sagittalXline');
  set(handles.sagittal_yline,'Color','b','EraseMode','xor','Tag','sagittalYline');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Contex Menus

menu=uicontextmenu;

% enable right click access to ROI tools
roi = uimenu(menu,'Label','ROI');
uimenu(roi,'Label','ROI (9x9x9 block)','Callback','avw_view([],[],''roi_9'');');
uimenu(roi,'Label','ROI (7x7x7 block)','Callback','avw_view([],[],''roi_7'');');
uimenu(roi,'Label','ROI (5x5x5 block)','Callback','avw_view([],[],''roi_5'');');
uimenu(roi,'Label','ROI (3x3x3 block)','Callback','avw_view([],[],''roi_3'');');

% save image to graphics file
uimenu(menu,'Label','Save Image','Callback','avw_view([],[],''save_image'');');

% zoom image to new figure
uimenu(menu,'Label','Zoom Image','Callback','avw_view([],[],''zoom'');');

if isfield(handles,'axial_image'),
  if isempty(get(handles.axial_image,'uicontextmenu')),
    set(handles.axial_image,'uicontextmenu',menu);
  end
end
if isfield(handles,'coronal_image'),
  if isempty(get(handles.coronal_image,'uicontextmenu')),
    set(handles.coronal_image,'uicontextmenu',menu);
  end
end
if isfield(handles,'sagittal_image'),
  if isempty(get(handles.sagittal_image,'uicontextmenu')),
    set(handles.sagittal_image,'uicontextmenu',menu);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image Intensity at Mouse Click

GUIheight = GUIheight - 0.04;

handles.Timval = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                           'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                           'BackgroundColor', [0 0 0],...
                           'ForegroundColor', [1 1 1],...
                           'BusyAction','queue',...
                           'String','Image Intensity');
handles.imval = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                          'Position',[.65 GUIheight .15 .03], 'HorizontalAlignment', 'right',...
                          'BackgroundColor', [0 0 0],...
                          'ForegroundColor', [1 1 1],...
                          'BusyAction','queue',...
                          'String','x','Value',0);

% Image Position at Mouse Click

GUIheight = GUIheight - 0.04;

handles.Timpos = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                           'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                           'BackgroundColor', [0 0 0],...
                           'ForegroundColor', [1 1 1],...
                           'BusyAction','queue',...
                           'String','Image Position');
handles.impos = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                          'Position',[.65 GUIheight .15 .03], 'HorizontalAlignment', 'right',...
                          'BackgroundColor', [0 0 0],...
                          'ForegroundColor', [1 1 1],...
                          'BusyAction','queue',...
                          'String','xyz','Value',[0 0 0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GUIheight = GUIheight - 0.04;

handles.flip = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                         'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                         'BackgroundColor', [.3 .3 .3],'ForegroundColor', [1 1 1],...
                         'BusyAction','queue',...
                         'TooltipString','Flip Left and Right (viewer only, see also avw_flip).',...
                         'String','Flip L/R',...
                         'Callback','avw_view([],[],''flip'');');
handles.flipStatus = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                               'Position',[.65 GUIheight .15 .03], 'HorizontalAlignment', 'right',...
                               'BackgroundColor', [0 0 0],...
                               'ForegroundColor', [1 1 1],...
                               'BusyAction','queue',...
                               'TooltipString','Flipped Status',...
                               'String','R>>L (radiological)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AC Location

GUIheight = GUIheight - 0.04;

handles.Tac = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                        'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                        'BackgroundColor', [.3 .3 .3],...
                        'ForegroundColor', [1 1 1],...
                        'BusyAction','queue',...
                        'TooltipString','AC point in (voxels) or (meter offset from center of volume)',...
                        'String','AC Point',...
                        'Callback','avw_view([],[],''ac'');');
handles.ac = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                       'Position',[.65 GUIheight .15 .03], 'HorizontalAlignment', 'right',...
                       'BackgroundColor', [0 0 0],...
                       'ForegroundColor', [1 1 1],...
                       'BusyAction','queue',...
                       'TooltipString','These values are offset from volume center.',...
                       'String','x,y,z');

% Nasion Location

GUIheight = GUIheight - 0.04;

handles.Tnasion = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                            'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                            'BackgroundColor', [.3 .3 .3],...
                            'ForegroundColor', [1 1 1],...
                            'BusyAction','queue',...
                            'TooltipString','Update Nasion - should be toward +Y',...
                            'String','Fiducial: Nas',...
                            'Callback','avw_view([],[],''nasion'');');
handles.nasion = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                           'Position',[.65 GUIheight .15 .03], 'HorizontalAlignment', 'right',...
                           'BackgroundColor', [0 0 0],...
                           'ForegroundColor', [1 1 1],...
                           'BusyAction','queue',...
                           'TooltipString','These values are offset from volume center, should be toward +Y',...
                           'String','x,y,z');

% Right Preauricular Location

GUIheight = GUIheight - 0.04;

handles.Trpa = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                         'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                         'BackgroundColor', [.3 .3 .3],...
                         'ForegroundColor', [1 1 1],...
                         'BusyAction','queue',...
                         'TooltipString','Update Right Preauricular - should be toward -X',...
                         'String','Fiducial: RPA',...
                         'Callback','avw_view([],[],''rpa'');');
handles.rpa = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                        'Position',[.65 GUIheight .15 .03], 'HorizontalAlignment', 'right',...
                        'BackgroundColor', [0 0 0],...
                        'ForegroundColor', [1 1 1],...
                        'BusyAction','queue',...
                        'TooltipString','These values are offset from volume center, should be toward -X',...
                        'String','x,y,z');

% Left Preauricular Location

GUIheight = GUIheight - 0.04;

handles.Tlpa = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                         'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                         'BackgroundColor', [.3 .3 .3],...
                         'ForegroundColor', [1 1 1],...
                         'BusyAction','queue',...
                         'TooltipString','Update Left Preauricular - should be toward +X',...
                         'String','Fiducial: LPA',...
                         'Callback','avw_view([],[],''lpa'');');
handles.lpa = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
                        'Position',[.65 GUIheight .15 .03], 'HorizontalAlignment', 'right',...
                        'BackgroundColor', [0 0 0],...
                        'ForegroundColor', [1 1 1],...
                        'BusyAction','queue',...
                        'TooltipString','These values are offset from volume center, should be toward +X',...
                        'String','x,y,z');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GUIheight = GUIheight - 0.04;

handles.contrast = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                             'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                             'BackgroundColor', [.3 .3 .3],'ForegroundColor', [1 1 1],...
                             'BusyAction','queue',...
                             'TooltipString','Auto contrast with gray colormap',...
                             'String','Auto Contrast',...
                             'Callback','avw_view([],[],''contrast'');');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GUIheight = GUIheight - 0.04;

handles.dimmer = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                           'Position',[.55 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                           'BackgroundColor', [.3 .3 .3],'ForegroundColor', [1 1 1],...
                           'BusyAction','queue',...
                           'TooltipString','Dim by 1%',...
                           'String','Dimmer',...
                           'Callback','avw_view([],[],''dimmer'');');

handles.brighter = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                             'Position',[.65 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                             'BackgroundColor', [.3 .3 .3],'ForegroundColor', [1 1 1],...
                             'BusyAction','queue',...
                             'TooltipString','Brighten by 1%',...
                             'String','Brighter',...
                             'Callback','avw_view([],[],''brighter'');');

handles.clim = uicontrol('Parent',GUI,'Style','edit','Units','Normalized', Font, ...
                         'Position',[.75 GUIheight .06 .03], 'HorizontalAlignment', 'right',...
                         'BackgroundColor', [0 0 0],...
                         'ForegroundColor', [1 1 1],...
                         'BusyAction','queue',...
                         'TooltipString','Image intensity Climit (inverse brightness)',...
                         'String',num2str(AVWVIEW.clim(2)),...
                         'Callback','avw_view([],[],''setClimit'');');

handles.cmap = uicontrol('Parent',GUI,'Style','popup','Units','Normalized', Font, ...
                         'Position',[.82 GUIheight .06 .03], 'HorizontalAlignment', 'left',...
                         'BackgroundColor', [0 0 0],...
                         'ForegroundColor', [1 1 1],...
                         'BusyAction','queue',...
                         'TooltipString','Color Map',...
                         'String',{'gray','bone','copper','hot','cool','spring','summer','autumn','winter','hsv','jet'},...
                         'Callback','avw_view([],[],''setCmap'');');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GUIheight = 0.46;

handles.crosshairs = uicontrol('Parent',GUI,'Style','checkbox','Units','Normalized', Font, ...
                               'Position',[.85 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                               'BackgroundColor', [.3 .3 .3],'ForegroundColor', [1 1 1],...
                               'BusyAction','queue',...
                               'TooltipString','Toggle Crosshairs on/off',...
                               'String','Crosshairs','Value',1,...
                               'Callback','avw_view([],[],''crosshairs'');');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GUIheight = GUIheight - 0.04;

handles.histogram = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                              'Position',[.85 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                              'BackgroundColor', [.3 .3 .3],'ForegroundColor', [1 1 1],...
                              'BusyAction','queue',...
                              'TooltipString','Histogram of Volume Intensity',...
                              'String','Histogram',...
                              'Callback','avw_view([],[],''histogram'');');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GUIheight = GUIheight - 0.04;

handles.coord = uicontrol('Parent',GUI,'Style','popupmenu','Units','Normalized', Font, ...
                          'Position',[.85 GUIheight .10 .03], 'HorizontalAlignment', 'left',...
                          'BackgroundColor', [.3 .3 .3],'ForegroundColor', [1 1 1],...
                          'BusyAction','queue',...
                          'TooltipString','Voxel or Mensurated Axis Coordinates',...
                          'String',{'Voxels','mm','meters'},...
                          'Callback','avw_view([],[],''coordinates'');');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Font.FontWeight = 'bold';

% View avw.hdr
handles.Bhdr = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                         'Position',[.92 .05 .07 .03],...
                         'String','HDR','BusyAction','queue',...
                         'TooltipString','View the .hdr parameters.',...
                         'BackgroundColor',[0.0 0.0 0.5],...
                         'ForegroundColor',[1 1 1], 'HorizontalAlignment', 'center',...
                         'Callback',strcat('AVWVIEW = get(gcbf,''Userdata''); ',...
                                           'avw_view_hdr(AVWVIEW.avw,AVWVIEW.gui);',...
                                           'clear AVWVIEW;'));

% OK: Return the avw!
handles.Bquit = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
                          'Position',[.92 .01 .07 .03],...
                          'String','RETURN','BusyAction','queue',...
                          'BackgroundColor',[0.0 0.5 0.0],...
                          'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
                          'Callback','avw_view([],[],''quit'');');

% Update the gui_struct handles for this gui
if exist('parent','var'), AVWVIEW.parent.gui = parent; end
AVWVIEW.avw = avw;
AVWVIEW.handles = handles;
set(AVWVIEW.gui,'Userdata',AVWVIEW);

AVWVIEW = set_coordinates(AVWVIEW);
set(AVWVIEW.gui,'HandleVisibility','callback');

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slice_img(avw),

figure
xslice = 128;
slice = squeeze( avw.img(xslice,:,:) );
imagesc(slice); axis image; colormap('gray')
figure
yslice = 128;
slice = squeeze( avw.img(:,yslice,:) );
imagesc(slice); axis image; colormap('gray')
figure
zslice = 128;
slice = squeeze( avw.img(:,:,zslice) );
imagesc(slice); axis image; colormap('gray')

return
