function	[d]=fmris_write_minc(d,Z,T);

if ~isfield(d,'dim') & isfield(d,'parent_file')
   d2=fmris_read_image(d.parent_file,0,0);
   d.dim=d2.dim;
end
if d.dim(4)==1
   d.dim(4)=0;
end

d.dim=d.dim(length(d.dim):-1:1);

if nargin<2
   Z=1:d.dim(2);
end   
if nargin<3
   T=1:(d.dim(1)+(d.dim(1)==0));
end

if isfield(d,'precision')
   precision=d.precision;
else
   precision='float';
end

if isfield(d,'file_path')
   file=[d.file_path d.file_name];
else
   file=d.file_name;
end

if T(1)==1 & Z(1)==1
   newh=newimage(file,d.dim,d.parent_file,precision);
else
   newh=openimage(file,'w');
end

data=squeeze(reshape(d.data,d.dim(4)*d.dim(3),length(Z),length(T)));

if d.dim(1)==0
   putimages(newh,data,Z);
elseif length(T)==1|length(Z)==1
   putimages(newh,data,Z,T);
else
   for i=1:length(T)
      putimages(newh,squeeze(data(:,:,i)),Z,T(i));
   end
end

closeimage(newh);

return
