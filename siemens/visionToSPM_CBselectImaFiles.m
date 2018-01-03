function visionToSPM_CBselectImaFiles( figHandle);

loadTitel = [...
      'select an ima-file of the volume       '; ...
      'select a mosaic ima of fMRI session    '; ...
      'select the FIRST !!! ima of volume-scan'];
objData.imaFileList = [];

imaType       = get( findobj( figHandle, 'Tag', 'imaType'), 'Value');
actionType    = get( findobj( figHandle, 'Tag', 'actionType'), 'Value');
verschachtelt = get( findobj( figHandle, 'Tag', 'verschachtelt'), 'Value');
imgFilename   = get( findobj( figHandle, 'Tag', 'imgFile'), 'String');
paramsMan     = get( findobj( figHandle, 'Tag', 'paramsManualy'), 'value');
delImas       = get( findobj( figHandle, 'Tag', 'delImas'), 'value');
%
% clear ima parameter / imagePreview
set( findobj( figHandle, 'Tag', 'imaHeaderInfo'), 'String', '');
set( findobj( figHandle, 'Tag', 'Axes1Image1'), 'Visible', 'off');

if visionToSPM_SHOWimgFile( figHandle)
   if isempty( imgFilename) 
      errordlg('please enter a filename for analyze file first !');
      return
   end
else
   imgFilename = [];
end
%
objData.imaType = imaType;
%   
% file open message depends on action and imaType
if actionType == 3 & imaType == 1 
   [objData.filename objData.path] = uigetfile('*.ima', loadTitel( 3,:));
else
   [objData.filename objData.path] = uigetfile('*.ima', loadTitel(imaType,:));
end

if ~objData.filename
   set( findobj( figHandle, 'Tag', 'selectImaFiles'), 'UserData', objData);
   set( findobj( figHandle, 'Tag', 'action'), 'visible', 'off');
   return;
end

% collect imas in objData.imaFileList
if actionType == 1 | actionType == 2
   objData.imaFileList = findFileList( objData.filename, objData.path, imaType);
   % let the user confirm the filelist
   if size( objData.imaFileList,1)==0
      errordlg('No continous ima list found ! This may due to a selection of imas numbered in mosaic style, while you selected NORMAL from ima type selection (and vice versa) !');
      objData.filename=0;
      set( findobj( figHandle, 'Tag', 'action'), 'visible', 'off');
      set( findobj( figHandle, 'Tag', 'selectImaFiles'), 'UserData', objData);
      return;      
   end
   msg = char( sprintf('continous list of %d ima''s found.', size( objData.imaFileList,1)),...
      sprintf('first ima: %s', char( objData.imaFileList(  1,:))),...
      sprintf('last  ima: %s', char( objData.imaFileList(end,:))),...
      'is this OK ?');
   if strcmp( questdlg( msg, 'visionToSPM', 'Yes','No','Yes'), 'No')
      objData.filename=0;
      set( findobj( figHandle, 'Tag', 'action'), 'visible', 'off');
      set( findobj( figHandle, 'Tag', 'selectImaFiles'), 'UserData', objData);
      return;
   end
end

showPreView( figHandle, objData.path, objData.filename);

if ~isempty( objData.imaFileList)
   params.file            = detImaParams( [objData.path char(objData.imaFileList(1,:)) ]);
else
   params.file            = detImaParams( [objData.path objData.filename]);
end
if paramsMan
   params.file = replaceParams( figHandle, params.file);
end
params.user.verschachtelt = verschachtelt;
params.user.imaType       = imaType;
params.user.actionType    = imaType;
params.user.delImas       = delImas;

if actionType==3
   [ imgFilename objData.path] = uigetfile('*.img', 'select img-file, transformation should applied to (cancel for none)');
   if imgFilename == 0 % file open menue canceled   
      imgFilename = 'M.mat';
      objData.path = '';
      fprintf('transform written to: %s\n', [pwd '\' imgFilename]);
   end
end;

showInfo( figHandle, params.file, objData.filename, objData.path);
set( findobj( figHandle, 'Tag', 'selectImaFiles'), 'UserData', objData);
set( findobj( figHandle, 'Tag', 'action'), 'visible', 'on');
return;