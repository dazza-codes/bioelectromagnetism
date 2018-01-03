function avw = LoadAVW ( Filename )
  fid = fopen ( Filename, 'r' );
  
  line = fgetl ( fid );
  if ~strmatch ( line, 'AVW_ImageFile' )
    error ( sprintf ( '%s is not in AVW format', Filename ) );
  end
  
  Info = sscanf ( line, 'AVW_ImageFile %f %d' );
  Offset = Info(2);
  
  avw = struct ( 'Image', 0 );
  
  % Read until 'EndInformation'
  while isempty ( strmatch ( line, 'EndInformation' ) )
    line = trim ( fgetl ( fid ) );
    if ~isempty ( strfind ( line, '=' ) )
      idx = strfind ( line, '=' );
      name = trim ( line(1:idx-1) );
      value = line(idx+1:end);
      
      if strfind ( name, '|' )
	name = ['DICOMTag_' strrep(name,'|','_')];
      end
      
      if ~isempty ( str2num ( value ) )
        value = str2num ( value );
      else
        % Strip out '"'
        value = strrep ( value, '"', '' );
      end
      
      name = trim ( validname ( name ) );
      try
        avw = setfield ( avw, name, value );
      catch
        % disp ( lasterr )
        warning ( sprintf ( 'Failed to set %s to %s', name, num2str ( value ) ) );
      end
    end
  end
  
  % How much do we read
  DataSize = avw.Width * avw.Height * avw.Depth * avw.NumVols;
  fclose ( fid );
  if strcmp ( avw.Endian, 'Little' )
    % Not exactly sure what's up with this being switched...
    fid = fopen ( Filename, 'r', 'l' );
  else
    fid = fopen ( Filename, 'r', 'b' );
  end
  fseek ( fid, Offset, 'bof' );
  % Figure out the datatype
  switch upper(avw.DataType)
    case 'AVW_UNSIGNED_SHORT'
     DataType = 'uint16';
    case 'AVW_SIGNED_SHORT'
     DataType = 'int16';
    case 'AVW_UNSIGNED_CHAR'
     DataType = 'uchar';
    case 'AVW_SIGNED_CHAR'
     DataType = 'schar';
    case 'AVW_FLOAT'
     DataType = 'float';
   otherwise
     error ( sprintf ( 'Unknown datatype: %s', avw.DataType ) );
  end


  Data = reshape ( fread ( fid, DataSize, DataType ), avw.Width, avw.Height, avw.Depth, avw.NumVols );
  % The image needs to be rotated
  avw.Image = zeros ( avw.Height, avw.Width, avw.Depth, avw.NumVols );
  yy = avw.Height:-1:1;
  for vol = 1:avw.NumVols
    for dd = 1:avw.Depth
      t = Data(:,:,dd,vol)';
      avw.Image(:,:,dd,vol) = t(yy,:);
    end
  end
  clear Data;
  fclose ( fid );
  avw.Image = squeeze ( avw.Image );
  
  
function sout = trim ( s )
  [r,c] = find ( ~isspace ( s ) );
  if isempty ( c )
    sout = s([]);
  else
    sout = s(:,min(c):end);
  end
  sout = deblank ( sout );
  
  
function s = validname ( s )
  ns = '';
  for ii = 1:length ( s )
    if ~isnumeric ( s(ii) )
      ns = [ns s(ii)];
    end
  end
  s = deblank ( ns );
  
