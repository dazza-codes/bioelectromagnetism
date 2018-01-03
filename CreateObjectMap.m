function objectmap = CreateObjectMap ( height, width, depth )
% objectmap = CreateObjectMap ( height, width, depth )
% Create an empty objectmap.

objectmap = struct ( 'Version', 910926 );
objectmap.Width = width;
objectmap.Height = height;
objectmap.Depth = depth;
objectmap.NumberOfObjects = 1;

objectmap.Objects = DefaultObject ( 'Original' );
objectmap.Image = ones ( objectmap.Height, objectmap.Width, objectmap.Depth );
