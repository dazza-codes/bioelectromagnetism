function showPreView( figHandle, path, filename);
%
% (c) Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%
% init session
imaFid = fopen( [ path filename], 'r', 'ieee-be');  
fseek( imaFid, 6144, 'bof');
[img n]= fread( imaFid, inf, 'uint16');

%if params.user.imaType == 1
%   if n ~= 256*256 
%      fprintf('unable to read image completly. %2.4f %% readed => terminating !', ...
%         n/256*256);
%      return;
%   end
%end

cMap = colormap;
set( findobj( figHandle, 'Tag', 'Axes1'), 'CLim', [ min(img) max(img)] );

set( figHandle, 'Colormap', gray);
s = sqrt(size(img,1));
set( findobj( figHandle, 'Tag', 'Axes1'), 'DataAspectRatio', [ s s 1 ]);
set( findobj( figHandle, 'Tag', 'Axes1'), 'XLim', [0 s]);
set( findobj( figHandle, 'Tag', 'Axes1'), 'YLIM', [0 s]);

set( findobj( figHandle, 'Tag', 'Axes1Image1'), 'CData', uint16( ( reshape(img, s, s)'  )) );
set( findobj( figHandle, 'Tag', 'Axes1Image1'), 'Visible', 'on');

set( figHandle, 'Colormap', cMap);

