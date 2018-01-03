function SaveAVW ( Filename, avw )
% SaveAVW ( Filename, avw )
% Save an AVW file from the given avw structure
  
  % Order is picky!
  HeaderNames = { 'DataType', 'Width', 'Height', 'Depth', 'NumVols', 'Endian', 'ColormapSize' };
  avw.MinimumDataValue = floor ( min ( avw.Image(:)) );
  avw.MaximumDataValue = floor ( max ( avw.Image(:)) );
  
  NonHeaderNames = GetNonHeaderNames ( avw, HeaderNames );
  
  HeaderInfo = FormatValues ( avw, HeaderNames, 0 );
  Info = FormatValues ( avw, NonHeaderNames, 1 );
  
  Header = HeaderInfo;
  Header = [Header sprintf('BeginInformation\n') Info sprintf('EndInformation\n') sprintf('MoreInformation=-1\n')];
  Header = [Header sprintf('Vol Slc  Offset    Length       Cmp Format\n')];
  Header = [Header sprintf('.CONTIG\n')];
  Header = [Header sprintf('EndSliceTable\n')];

  HeaderSize = 4096 * ceil ( (100+length(Header)) / 4096 );
  
  Header = [sprintf('AVW_ImageFile   1.00     %d\n',HeaderSize) Header];
  l = length ( Header );
  b = strrep ( blanks ( HeaderSize - l ), ' ', '#' );
  Header = [Header b];
  
  Endian = 'b';
  if strcmp ( avw.Endian, 'Little' )
    Endian = 'l';
  end
  fid = fopen ( Filename, 'w', Endian );
  % Watch out, if we don't use %s, escape sequences get mucked up
  fprintf ( fid, '%s', Header );
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
  
  % Need to undo the flip flop
  yy = size(avw.Image,1):-1:1;
  for vol = 1:avw.NumVols
    for dd = 1:size(avw.Image,3)
      t = squeeze(avw.Image(:,:,dd,vol))';
      t = t(:,yy);
      fwrite ( fid, t(:), DataType );
    end
  end
  fclose ( fid );
  
  
function Header = FormatValues ( avw, Names, Quote )
  Header = '';
  for hidx = 1:length(Names)
    Name = Names{hidx};
    Value = getfield ( avw, Name );
    
    if ~isempty ( findstr ( Name, 'DICOMTag_' ) )
      idx = strfind ( Name, '_' );
      idx = idx + 1;
      Name = [Name(idx(1):idx(1)+3) '|' Name(idx(2):idx(2)+3)];
    end
    
    if isstr ( Value ) & Quote
      Value = sprintf ( '"%s"', Value );
    end
    if isnumeric ( Value )
      Value = num2str ( Value );
    end
    if Quote
      t = sprintf ( '  %s=%s\n', Name, Value );
    else
      t = sprintf ( ' %s=%s\n', Name, Value );
    end      
    Header = [Header t];
  end
  
  
function n = GetNonHeaderNames ( avw, HeaderNames )
  names = fieldnames ( avw );
  n = {};
  for nn = 1:length(names)
    if strcmp ( names{nn}, 'Image' )
      continue
    end
    
    IsHeader = 0;
    for hidx = 1:length(HeaderNames)
      if strcmp ( names{nn}, HeaderNames{hidx} )
	IsHeader = 1;
      end
    end
    
    if ~IsHeader
      n{length(n)+1} = names{nn};
    end
  end
  
      
  