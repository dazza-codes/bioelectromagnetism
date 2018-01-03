function visionToSPM_CBgo( figHandle);
%
% (c) Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%
% init session
file=get( findobj( figHandle, 'Tag', 'selectImaFiles'), 'UserData');

if ~file.filename
   errordlg('please select a list of ima files first!');
   return;
end
%
% read user selection
imaType       = get( findobj( figHandle, 'Tag', 'imaType'), 'Value');
actionType    = get( findobj( figHandle, 'Tag', 'actionType'), 'Value');
verschachtelt = get( findobj( figHandle, 'Tag', 'verschachtelt'), 'Value');
imgFilename   = get( findobj( figHandle, 'Tag', 'imgFile'), 'String');
paramsMan     = get( findobj( figHandle, 'Tag', 'paramsManualy'), 'value');
delImas       = get( findobj( figHandle, 'Tag', 'delImas'), 'value');
matFilePrefix = get( findobj( figHandle, 'Tag', 'matFilePrefix'), 'String');

if visionToSPM_SHOWimgFile( figHandle)
   if isempty( imgFilename) 
      errordlg('please enter a filename for analyze file first !');
      return
   end
else
   imgFilename = [];
end

if ~isempty( file.imaFileList)
   params.file            = detImaParams( [file.path char(file.imaFileList(1,:)) ]);
else
   params.file            = detImaParams( [file.path file.filename]);
end
if paramsMan
   params.file = replaceParams( figHandle, params.file);
end

params.user.verschachtelt = verschachtelt;
params.user.imaType       = imaType;
params.user.actionType    = imaType;
params.user.delImas       = delImas;

switch actionType
   case 1
      err = 0;
      h = waitbar(0,'starting bitmap generation ...');
      err = imaToImg_bitmap( imgFilename, file.imaFileList, file.path, params, h);
      close( h);
      if imaType == 1
         err = imaToImg_geometrie( imgFilename, file.path, params);
         fprintf('processing: geometry\n');
      else
         h = waitbar(0,'calculating geometrical transforms (*.mat files) ...');
         fprintf('processing: geometry\n');
         cDir = pwd; cd( file.path);         
         fn1 = sprintf( '%s', char( file.imaFileList(1,:)));
         fn1 = [ matFilePrefix fn1(1:end-4) ];
         % determine the transformation from first ima
         err = imaToImg_geometrie( fn1, file.path, params);
         % and copy the rest ...
         fid = fopen( [fn1 '.mat'],'r','ieee-le'); F = fread( fid); fclose(fid); 
         for i=2:size( file.imaFileList,1)
            fn2 = sprintf( '%s', char( file.imaFileList(i,:)) );
            fn2 = [ matFilePrefix strrep( fn2, '.ima', '.mat') ]; 
            try
               % much more faster than copyfile
               fid = fopen(fn2,'w','ieee-le'); fwrite(fid,F); fclose( fid);
            catch
               err=1;
            end            
            waitbar(i/size( file.imaFileList,1),h);   
         end
         cd( cDir);
         close( h);
      end
            
      if err == 1
         fprintf('\nan error occured during img file creation, imas not deleted !\n\n');
      else
         if delImas
            cDir = pwd; cd( file.path);
            for i = 1 : size(file.imaFileList,1);
               try
                  delete( char( file.imaFileList(i,:)));
               catch
               end
            end
            cd(cDir);
         end
      end
      
   case 2
      h = waitbar(0,'starting bitmap eneration ...');
      err = imaToImg_bitmap( imgFilename, file.imaFileList, file.path, params, h);
      close( h);
      
      if err == 1
         fprintf('\nan error occured during img file creation, imas not deleted !\n\n');
      else
         if delImas
            for i = 1 : size(file.imaFileList,1);
               try
                  delete( char( file.imaFileList(i,:)));
               catch
               end
            end
         end
      end      
      
   case 3
      [ imgFilename file.path] = uigetfile('*.img', 'select img-file, transformation should applied to (cancel for none)');
      if imgFilename == 0 % file open menue canceled   
         imgFilename = 'M.mat';
         file.path = '';
         fprintf('transform written to: %s\n', [pwd '\' imgFilename]);
      end
      imaToImg_geometrie( [ matFilePrefix imgFilename(1:end-4) ], file.path, params); 
   case 4
end
cd( file.path);
fprintf('finished !\n\n');
return;