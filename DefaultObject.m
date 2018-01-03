function O = DefaultObject ( Name )
% O = DefaultObject ( Name )
% Create an object for an AVWObjectMap.  Creates a color
% pseudo-randomly using the Name as a seed, so that the same name
% will always give the same color.
  
  S = rand ( 'state' );
  rand ( 'seed', sum ( double ( Name ) ) );
% Return an empty Object
  Object = struct ( 'Name', Name );
  Object.Display = 1;
  Object.Copy = 0;
  Object.Mirror = 0;
  Object.Status = 0;
  Object.NUsed = 1;
  Object.Shades = 64;
  Object.StartRed = 0;
  Object.StartGreen = 0;
  Object.StartBlue = 0;
  Object.EndRed = 128 + round ( 127 * rand ( 1 ) );
  Object.EndGreen = 128 + round ( 127 * rand ( 1 ) );
  Object.EndBlue = 128 + round ( 127 * rand ( 1 ) );
  Object.XRotation = 0;
  Object.YRotation = 0;
  Object.ZRotation = 0;
  Object.XShift = 0;
  Object.YShift = 0;
  Object.ZShift = 0;
  Object.XCenter = 0;
  Object.YCenter = 0;
  Object.ZCenter = 0;
  Object.XRotationIncrement = 0;
  Object.YRotationIncrement = 0;
  Object.ZRotationIncrement = 0;
  Object.XShiftIncrement = 0;
  Object.YShiftIncrement = 0;
  Object.ZShiftIncrement = 0;
  Object.XMinimum = 0;
  Object.YMinimum = 0;
  Object.ZMinimum = 0;
  Object.XMaximum = 0;
  Object.YMaximum = 0;
  Object.ZMaximum = 0;
  Object.Opacity = 1.0;
  Object.OpacityThickness = 0;
  O = Object;
  rand ( 'state', S );