function visionToSPM_CBdelImas( figHandle);

if visionToSPM_SHOWdelImas( figHandle)
   set( findobj( figHandle, 'Tag', 'delImas'), 'Visible', 'off');
else
   set( findobj( figHandle, 'Tag', 'delImas'), 'Visible', 'on');
end
