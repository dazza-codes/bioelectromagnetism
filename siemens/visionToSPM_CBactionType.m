function visionToSPM_CBactionType( figHandle)

if visionToSPM_SHOWimgFile( figHandle)
   set( findobj( figHandle, 'Tag', 'imgFile'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'imgFileDesc'), 'Visible', 'on');
else
   set( findobj( figHandle, 'Tag', 'imgFile'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'imgFileDesc'), 'Visible', 'off');
end   

if visionToSPM_SHOWdelImas( figHandle)
   set( findobj( figHandle, 'Tag', 'delImas'), 'Visible', 'on');
else
   set( findobj( figHandle, 'Tag', 'delImas'), 'Visible', 'off');
end
