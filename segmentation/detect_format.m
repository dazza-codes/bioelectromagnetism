function format=detect_format(filename)
% format=detect_format(filename)
% tries to detect the format of filename based on the header size 348
% the very first  number of the header.

% rgpm 19-7-2000

format='native';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end

format='cray';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end

format='ieee-be';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end

format='ieee-le';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end

format='ieee-be.164';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end

format='ieee-le.164';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end

format='vaxd';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end

format='vaxg';
fid=fopen(filename,'r',format);
nada= fread(fid,1,'uint32');
fclose(fid);
if nada==348, return, end
