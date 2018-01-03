function status = GE_writeSPMHeader(fname,header)
%
% status = GE_writeSPMHeader(fname,header)
%
% Write SPM header into fname.hdr
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           %
% Write Analyze header file %
%                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fid,message] = fopen(fname,'w');

if (fid == -1)
	disp(message);
	return;
end

% write ANALYZE Header Key
%---------------------------------------------------------------------------
fseek(fid,0,'bof');
fwrite(fid,header.sizeOfHeader,		'uint32');
fwrite(fid,header.dataType,	        'char' );
fwrite(fid,header.dbName,		'char');
fwrite(fid,header.extents,		'int32');
fwrite(fid,header.sessionError,		'int16' );
fwrite(fid,header.regular,		'char' );
fwrite(fid,header.hkeyUnused,		'char' );

% write ANALYZE Image Dimension
%---------------------------------------------------------------------------
fseek(fid,40,'bof');
fwrite(fid,header.dim,			'int16');
fwrite(fid,header.voxUnits,		'char' );
fwrite(fid,header.calUnits,		'char' );
fwrite(fid,header.iUnused1,		'int16' );
fwrite(fid,header.datatype,		'int16' );
fwrite(fid,header.bitpix,		'int16');
fwrite(fid,header.dimUnused,		'int16');
fwrite(fid,header.pixDim,		'float32');
fwrite(fid,header.voxOffset,		'float32');
fwrite(fid,header.fUnused1,		'float32');
fwrite(fid,header.fUnused2,		'float32');
fwrite(fid,header.fUnused3,		'float32');
fwrite(fid,header.calMax,		'float32');
fwrite(fid,header.calMin,		'float32');
fwrite(fid,header.compressed,		'float32');
fwrite(fid,header.verified,		'float32');
fwrite(fid,header.glMax,		'int32');
fwrite(fid,header.glMin,		'int32');

% write ANALYZE data history
%---------------------------------------------------------------------------
fwrite(fid,header.description,	'char');
fwrite(fid,header.auxFile,	'char');
fwrite(fid,header.orient,	'int8');
fwrite(fid,header.originator,	'int16');
fwrite(fid,header.generated,	'char');
fwrite(fid,header.scanNum,	'char');
fwrite(fid,header.patientID,	'char');
fwrite(fid,header.expDate,	'char');
fwrite(fid,header.expTime,	'char');
fwrite(fid,header.histUN0,	'char');
fwrite(fid,header.views,	'int32');
fwrite(fid,header.volsAdded,	'int32');
fwrite(fid,header.startField,	'int32');
fwrite(fid,header.fieldSkip,	'int32');
fwrite(fid,header.oMax,		'int32');
fwrite(fid,header.oMin,		'int32');
fwrite(fid,header.sMax,		'int32');

if fwrite(fid,header.sMin,	'int32') ~= 1
	fclose(fid);
	error(['Error writing ' fname '. Check your disk space.']);
end

if ftell(fid) ~= -1
  status = 1;
else
  status = 0;
end

fclose(fid);

return






