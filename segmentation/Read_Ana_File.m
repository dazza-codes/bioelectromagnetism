function [data,head,format]=Read_Ana_File(filename)
% [data,header]=Read_Ana_File(filename)
% READs ANAlize FILE
% filename without extension. It is assumed that there are two files on disk:
% filename.hdr and filename.img.

% rgpm 19-7-2000

if exist([filename,'.hdr'])~=2 ,
   myerror(' Missing Header file / in Read_Ana_File ');
end

if exist([filename,'.img'])~=2 ,
   myerror(' Missing Image file / in Read_Ana_File ');
end

[header,formato]=Read_Ana_Header([filename,'.hdr']);

Ndim=header.dime.dim(1);  % number of dimensions
Nx=header.dime.dim(2);    % pixels in X
Ny=header.dime.dim(3);    % pixels in Y
Nz=header.dime.dim(4);    % pixels in Z slices
Nt=header.dime.dim(5);    % number of time frames

if Ndim >4, myerror(' Option not available in Read_Ana_File ' ); end

filetype=Ana_data_file_type(header.dime.datatype,[]);

disp(' Read_Ana_File -->  dimensions :'); [Nx,Ny,Nz,Nt]

data=Read_Binary_Vol4([filename,'.img'],Nx,Ny,Nz,Nt,filetype,formato);

if nargout>=2,
   head=header;
end
if nargout==3,
   format=formato;
end

%over   

return

% data type !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
switch header.dime.bitpix,
case 32, filetype='uint32',   
case 16, filetype='uint16', 
case  8, filetype='uint8',
otherwise
   myerror('File type option not available in Read_Ana_File ');
end
