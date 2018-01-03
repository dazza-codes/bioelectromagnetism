function [d]=fmris_read_minc(file,Z,T)

%file_l    		= length(file);
sep = '/';
sep_list 		= findstr(file,sep);

if ~isempty(sep_list)
   %d.file_path 	= file(1:sep_list(length(sep_list)));
   %d.file_name 	= file((sep_list(length(sep_list))+1):length(file));
   d.file_path 	= '';
   d.file_name 	= file;
else
   d.file_path 	= [pwd sep];
   d.file_name 	= file;
end

d.dim 	= fliplr(miinquire(file,'imagesize')');
d.dim 	= d.dim + (d.dim == 0);

if nargin == 1 
   Z = 1:d.dim(3);
   T = 1:d.dim(4);
end

if (T~=0)|(Z~=0)
   if d.dim(4)==1
      images 	= mireadimages(file,Z-1);
   else
      n=length(T);
      images=[];
      for i=0:floor((n-1)/256)
         frames=T((i*256+1):min((i+1)*256,n));
         images 	= [images mireadimages(file,Z-1,frames-1)];
      end
   end
   d.data 		= reshape(images,[d.dim(1:2) length(Z) length(T)]);
   d.calib		= [min(min(min(min(d.data)))) max(max(max(max(d.data))))];
end

[x_vox, y_vox, z_vox] 	=  miinquire(file,'attvalue','xspace','step','attvalue','yspace','step','attvalue','zspace','step');
d.vox 					= [x_vox y_vox z_vox];
d.vox_units				= '';
d.vox_offset			= 0;
d.precision				= '';
d.calib_units			= '';
[x_origin, y_origin, z_origin] = miinquire(file,'attvalue','xspace','start','attvalue','yspace','start','attvalue','zspace','start');
d.origin 				= [x_origin y_origin z_origin];
d.descrip				= '';

return;

