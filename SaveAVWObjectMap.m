function SaveAVWObjectMap ( Filename, objectmap )
% SaveAVWObjectMap ( Filename, objectmap )
% Save the objectmap in Filename

fid = fopen ( Filename, 'w', 'b' );
if fid < 0,
  error ( ['Can not open ' Filename ' for writing'] );
end

fwrite ( fid, objectmap.Version, 'int' );
if objectmap.Version ~= 910926
  error ( ['Unknown ObjectMap version: ' int2str(objectmap.Version) ' expected 910926'] );
end

fwrite ( fid, objectmap.Width, 'int' );
fwrite ( fid, objectmap.Height, 'int' );
fwrite ( fid, objectmap.Depth, 'int' );
fwrite ( fid, objectmap.NumberOfObjects, 'int' );

for o = 1:objectmap.NumberOfObjects
  Object = objectmap.Objects(o);
  N = double ( Object.Name );
  N = [N zeros(1,32-length(N))];
  fwrite ( fid, N, 'uchar' );
  fwrite ( fid, Object.Display, 'int' );
  fwrite ( fid, Object.Copy, 'uchar' );
  fwrite ( fid, Object.Mirror, 'uchar' );
  fwrite ( fid, Object.Status, 'uchar' );
  fwrite ( fid, Object.NUsed, 'uchar' );
  fwrite ( fid, Object.Shades, 'int' );
  fwrite ( fid, Object.StartRed, 'int' );
  fwrite ( fid, Object.StartGreen, 'int' );
  fwrite ( fid, Object.StartBlue, 'int' );
  fwrite ( fid, Object.EndRed, 'int' );
  fwrite ( fid, Object.EndGreen, 'int' );
  fwrite ( fid, Object.EndBlue, 'int' );
  fwrite ( fid, Object.XRotation, 'int' );
  fwrite ( fid, Object.YRotation, 'int' );
  fwrite ( fid, Object.ZRotation, 'int' );
  fwrite ( fid, Object.XShift, 'int' );
  fwrite ( fid, Object.YShift, 'int' );
  fwrite ( fid, Object.ZShift, 'int' );
  fwrite ( fid, Object.XCenter, 'int' );
  fwrite ( fid, Object.YCenter, 'int' );
  fwrite ( fid, Object.ZCenter, 'int' );
  fwrite ( fid, Object.XRotationIncrement, 'int' );
  fwrite ( fid, Object.YRotationIncrement, 'int' );
  fwrite ( fid, Object.ZRotationIncrement, 'int' );
  fwrite ( fid, Object.XShiftIncrement, 'int' );
  fwrite ( fid, Object.YShiftIncrement, 'int' );
  fwrite ( fid, Object.ZShiftIncrement, 'int' );
  fwrite ( fid, Object.XMinimum, 'short' );
  fwrite ( fid, Object.YMinimum, 'short' );
  fwrite ( fid, Object.ZMinimum, 'short' );
  fwrite ( fid, Object.XMaximum, 'short' );
  fwrite ( fid, Object.YMaximum, 'short' );
  fwrite ( fid, Object.ZMaximum, 'short' );
  fwrite ( fid, Object.Opacity, 'float' );
  fwrite ( fid, Object.OpacityThickness, 'int' );
  % Dummy
  fwrite ( fid, 0, 'int' );
end


Output = 0;
NOut = 0;

% RLE is done image by image
for slice = 1:size(objectmap.Image,3)
  X = squeeze(objectmap.Image(:,:,slice));
  % Write the volume run length encoded.
  % Rotate back
  % X = rot90 ( X, 2 );
  X = X(:)';
  Counts = diff([ 0 find(X(1:end-1) ~= X(2:end)) length(X) ]);
  Values = X([ find(X(1:end-1) ~= X(2:end)) length(X) ]);

  for idx = 1:length(Counts)
    Count = Counts(idx);
    Value = Values(idx) - 1;
    while Count >= 255
      fwrite ( fid, [255 Value], 'uchar' );
      Count = Count - 255;
      % Output = Output + 255;
      % NOut = NOut + 1;
    end
    fwrite ( fid, [Count Value], 'uchar' );
    % Output = Output + Count;
    % NOut = NOut + 1;
  end
end
fclose ( fid );
