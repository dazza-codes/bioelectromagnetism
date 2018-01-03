function [ana,formato]=Read_Ana_Header(filename)
% [ana,formato]=Read_Ana_Header(filename)
% Reads the header of analize file filename.hdr
% Note that filename can contain float (4 bytes real) data that are not allowed 
% in MatLab.
% References:  anz_hdr_struct.pro (IDL program from Ivan Zimini) and
% www.mayo.edu/bir/analyze/analyze_ov.html

% rgpm 19-7-2000

format=detect_format(filename), % detectando formato

if nargout==2, formato=format; end

[fid,message] = fopen(filename,'r',format);

if fid==-1, myerror(message), end

%            header_key

nada= fread(fid,1,'uint32');   
header_key.sizeof_hdr =nada;

nada= fread(fid,[1,10],'uchar');   
header_key.data_type=char(nada);

nada= fread(fid,[1,18],'uchar');   
header_key.db_name=char(nada);

nada= fread(fid,1,'uint32');   
header_key.extents=nada;

nada= fread(fid,1,'uint16');   
header_key.session_error=nada;

nada= fread(fid,1,'uchar');   
header_key.regular=char(nada);

nada= fread(fid,1,'uchar');   
header_key.hkey_un0=char(nada);

%               image_dimensions

nada= fread(fid,[1,8],'uint16');   
image_dimensions.dim=nada;

nada= fread(fid,[1,4],'uchar');   
image_dimensions.vox_units=char(nada);

nada= fread(fid,[1,8],'uchar');   
image_dimensions.cal_units=char(nada);

nada= fread(fid,1,'uint16');   
image_dimensions.unused1=nada;

nada= fread(fid,1,'uint16');   
image_dimensions.datatype=nada;

nada= fread(fid,1,'uint16');   
image_dimensions.bitpix=nada;

nada= fread(fid,1,'uint16');   
image_dimensions.dim_un0=nada;

nada= fread(fid,[1,8],'float32');   
image_dimensions.pixdim=nada; % in disk it is a float !!!!!!

nada= fread(fid,1,'float32');   
image_dimensions.vox_offset=nada; % in disk it is a float !!!!!!

nada= fread(fid,1,'float32');   
image_dimensions.funused1=nada;  % in disk it is a float !!!!!!

nada= fread(fid,1,'float32');   
image_dimensions.funused2=nada; % in disk it is a float !!!!!!

nada= fread(fid,1,'float32');   
image_dimensions.funused3=nada; % in disk it is a float !!!!!!

nada= fread(fid,1,'float32');   
image_dimensions.cal_max=nada; % in disk it is a float !!!!!!

nada= fread(fid,1,'float32');   
image_dimensions.cal_min=nada; % in disk it is a float !!!!!!

nada= fread(fid,1,'uint32');   
image_dimensions.compressed=nada;

nada= fread(fid,1,'uint32');   
image_dimensions.verified=nada;

nada= fread(fid,1,'uint32');   
image_dimensions.glmax=nada;

nada= fread(fid,1,'uint32');   
image_dimensions.glmin=nada;

%                          data_history

nada= fread(fid,[1,80],'uchar');   
data_history.descrip=char(nada);

nada= fread(fid,[1,24],'uchar');   
data_history.aux_file=char(nada);

nada= fread(fid,1,'uchar'); % segun SPM   
data_history.orient=nada;

nada= fread(fid,[1,5],'int16');  % segun SPM 
data_history.originator=nada;

nada= fread(fid,[1,10],'uchar');   
data_history.generated=char(nada);

nada= fread(fid,[1,10],'uchar');   
data_history.scannum=char(nada);

nada= fread(fid,[1,10],'uchar');   
data_history.patient_id=char(nada);

nada= fread(fid,[1,10],'uchar');   
data_history.exp_date=char(nada);

nada= fread(fid,[1,10],'uchar');   
data_history.exp_time=char(nada);

nada= fread(fid,[1,3],'uchar');   
data_history.hist_un0=char(nada);

nada= fread(fid,1,'uint32');   
data_history.views=nada;

nada= fread(fid,1,'uint32');   
data_history.vols_added=nada;

nada= fread(fid,1,'uint32');   
data_history.start_field=nada;

nada= fread(fid,1,'uint32');   
data_history.field_skip=nada;

nada= fread(fid,1,'uint32');   
data_history.omax=nada;

nada= fread(fid,1,'uint32');   
data_history.omin=nada;

nada= fread(fid,1,'uint32');   
data_history.smax=nada;

nada= fread(fid,1,'uint32');   
data_history.smin=nada;

fclose(fid);

ana.hk = header_key;
ana.dime = image_dimensions;
ana.hist = data_history;

