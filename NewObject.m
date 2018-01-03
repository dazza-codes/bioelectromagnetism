function objectmap = NewObject ( objectmap, Name )
% Index = NewObject ( objectmap, Name )
% Add a new object to the objectmap called Name and return it's index.

  Index = length ( objectmap.Objects ) + 1;
  if length(Name) > 31
    Name = Name(1:31);
  end
  objectmap.Objects(Index) = DefaultObject ( Name );
  objectmap.NumberOfObjects = length ( objectmap.Objects );
  return

  
