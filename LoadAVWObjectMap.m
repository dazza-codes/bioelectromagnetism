function objectmap = LoadAVWObjectMap ( Filename )
% objectmap = LoadAVWObjectMap ( Filename )
% Load the objectmap given by Filename.  Returns a structure
% containing the map (as field Image), and a list of objects (as
% field Objects)
%
% NOTE: Matlab versions have an offset of 1, i.e. the 1-based
% indexing of Matlab applies to the objectmap.Image as well as the
% objectmap.Objects.  That is objectmap.Objects(1) is defined as
% all the voxels in objectmap.Image == 1

  % Filename = 'BWLabelsHippocampus.obj';

objectmap = struct ( 'Version', -1 );  
  
fid = fopen ( Filename, 'r', 'b');
if fid < 0,
  error ( ['Can not open ' Filename ' for reading'] );
end

objectmap.Version = fread ( fid, 1, 'int' );
if ((objectmap.Version ~= 910926) && (objectmap.Version ~=20050829))
  error ( ['Unknown ObjectMap version: ' int2str(objectmap.Version) ' expected 910926 or 20050829'] );
end

objectmap.Width = fread ( fid, 1, 'int' );
objectmap.Height = fread ( fid, 1, 'int' );
objectmap.Depth = fread ( fid, 1, 'int' );
objectmap.NumberOfObjects = fread ( fid, 1, 'int' );
% 20050829 versions added one int to the stream to handle 4D data
% Should be the number of volumes.
if isequal(objectmap.Version,20050829)
    objectmap.NumberOfVolumes=fread(fid,1,'int');
else
    objectmap.NumberOfVolumes=1;
end
for o = 1:objectmap.NumberOfObjects
  Object = struct ( 'Name', 'Foo' );
  Object.Name = char(fread ( fid, 32, 'uchar' ));
  Object.Name = sscanf ( Object.Name, '%s' );
  Object.Display = fread ( fid, 1, 'int' );
  Object.Copy = fread ( fid, 1, 'uchar' );
  Object.Mirror = fread ( fid, 1, 'uchar' );
  Object.Status = fread ( fid, 1, 'uchar' );
  Object.NUsed = fread ( fid, 1, 'uchar' );
  Object.Shades = fread ( fid, 1, 'int' );
  Object.StartRed = fread ( fid, 1, 'int' );
  Object.StartGreen = fread ( fid, 1, 'int' );
  Object.StartBlue = fread ( fid, 1, 'int' );
  Object.EndRed = fread ( fid, 1, 'int' );
  Object.EndGreen = fread ( fid, 1, 'int' );
  Object.EndBlue = fread ( fid, 1, 'int' );
  Object.XRotation = fread ( fid, 1, 'int' );
  Object.YRotation = fread ( fid, 1, 'int' );
  Object.ZRotation = fread ( fid, 1, 'int' );
  Object.XShift = fread ( fid, 1, 'int' );
  Object.YShift = fread ( fid, 1, 'int' );
  Object.ZShift = fread ( fid, 1, 'int' );
  Object.XCenter = fread ( fid, 1, 'int' );
  Object.YCenter = fread ( fid, 1, 'int' );
  Object.ZCenter = fread ( fid, 1, 'int' );
  Object.XRotationIncrement = fread ( fid, 1, 'int' );
  Object.YRotationIncrement = fread ( fid, 1, 'int' );
  Object.ZRotationIncrement = fread ( fid, 1, 'int' );
  Object.XShiftIncrement = fread ( fid, 1, 'int' );
  Object.YShiftIncrement = fread ( fid, 1, 'int' );
  Object.ZShiftIncrement = fread ( fid, 1, 'int' );
  Object.XMinimum = fread ( fid, 1, 'short' );
  Object.YMinimum = fread ( fid, 1, 'short' );
  Object.ZMinimum = fread ( fid, 1, 'short' );
  Object.XMaximum = fread ( fid, 1, 'short' );
  Object.YMaximum = fread ( fid, 1, 'short' );
  Object.ZMaximum = fread ( fid, 1, 'short' );
  Object.Opacity = fread ( fid, 1, 'float' );
  Object.OpacityThickness = fread ( fid, 1, 'int' );
  Dummy = fread ( fid, 1, 'int' );
  objectmap.Objects(o) = Object;
end

% Read the rest of the file
Map = fread ( fid, inf, 'uchar' );
% objectmap.RLE = Map;
% Undo RLE
T = zeros ( objectmap.Width, objectmap.Height, objectmap.Depth );

ImageIndex = 1;
Counts = Map(1:2:end);
Values = Map(2:2:end);

for idx = 1:length(Counts)
    try
        T(ImageIndex:ImageIndex+Counts(idx)-1) = Values(idx) + 1;
        ImageIndex = ImageIndex + Counts(idx);
    catch
            idx
    ImageIndex
    end
end

objectmap.Image = zeros ( objectmap.Height, objectmap.Width, objectmap.Depth );
yy = objectmap.Height:-1:1;
for dd = 1:objectmap.Depth
  t = T(:,:,dd)';
  objectmap.Image(:,:,dd) = t(yy,:);
end

fclose ( fid );

