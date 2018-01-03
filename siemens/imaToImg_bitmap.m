function err = imaToImg_bitmap( imgFilename, imaFilelist, path, params, hProgBar);
%
% Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%

% the hProgBar is an optional parameter, a call like
% imaToImg_bitmap( imgFilename, imaFilelist, path, params) is also possible !

switch params.user.imaType
   case 1 % normal
      if nargin == 5
         err = simpleIma2img( imgFilename, imaFilelist, path, params, hProgBar);
      else
         err = simpleIma2img( imgFilename, imaFilelist, path, params);
      end
      
   case 2 % mosaic
      if nargin == 5
         err = mosaicIma2img( imaFilelist, path, params, hProgBar);
      else
         err = mosaicIma2img( imaFilelist, path, params);
      end
      
end
