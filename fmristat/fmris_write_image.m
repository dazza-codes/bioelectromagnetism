function [d]=fmris_write_image(d,Z,T)

% (c) John Aston & Roger Gunn

switch lower(d.file_name(min(findstr('.',d.file_name))+(1:3)))
case 'mnc'
   if nargin == 3
      fmris_write_minc(d,Z,T);
   else
      fmris_write_minc(d);
   end
case 'img'
   if nargin == 3
      fmris_write_analyze(d,Z,T);
   else
      fmris_write_analyze(d);
   end
end

return


