function vol4=Read_Binary_Vol4(filename,Nx,Ny,Nz,Nt,filetype,formato)
% vol4=Read_Binary_Vol4(filename,Nx,Ny,Nz,Nt,filetype,formato)
% Reads Binary data corresponding with a 4-D array.
% To read Analyze format or just binary saved data
% Nt can be one as it is for anatomical data (MRI),
% or bigger than one as for fMRI data.
% filetype correspond to precision, see fread Matlab command.
% Ouput matrix vol4 is uint16 (2 bytes unsigned integer) 

% rgpm 19-7-2000
% note 13-2-2003 According to Ivan X and Y are changed with respect to the
% way matlab stores arrays, however I find no problem at all. 


if nargin<7, formato='native'; end 

[fid,message] = fopen(filename,'r',formato);
if fid==-1, myerror([message, ' / Read_Binary_Vol4 ']), end

vol4=repmat(uint16(1),[Nx,Ny,Nz,Nt]);

Nxy=Nx*Ny;
h=waitbar(0,['Reading data from --> ',filename]);
for t=1:Nt,
   if Nt>1,waitbar(t/Nt,h); end 
   for z=1:Nz,
      if Nt==1, waitbar(z/Nz); end
      [temp,cont]=fread(fid,[Nx,Ny],filetype);
      %temp, pause
      if cont~=Nxy, myerror(' Read_Binary_Vol4 '); end
      %vol4(:,:,z,t)=convert(temp,'uint16');
      vol4(:,:,z,t)=temp;
   end
end
close(h)
fclose(fid);


%over

