function [d] = fmris_read_image(file,Z,T)

% (c) John Aston & Roger Gunn

file = deblank(file);

if strcmp(file((length(file)-2):length(file)),'.gz')
   unix(['gunzip ' file]);
   file=file(1:(length(file)-3));
end

switch lower(file((length(file)-3):length(file)))
case '.mnc'
   if nargin == 3
      d = fmris_read_minc(file,Z,T);
   else
      d = fmris_read_minc(file);
   end
case '.img'
   if nargin == 3
      d = fmris_read_analyze(file,Z,T);
   else
      d = fmris_read_analyze(file);
   end
end

if ~isfield(d,'parent_file')
   d.parent_file=file;
end

return