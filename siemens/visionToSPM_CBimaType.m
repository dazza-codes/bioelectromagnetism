function visionToSPM_CBimaType( figHandle)

file=get( findobj( figHandle, 'Tag', 'selectImaFiles'), 'UserData');
iT = get( findobj( figHandle, 'Tag', 'imaType'), 'Value');

if visionToSPM_SHOWimgFile( figHandle)
   set( findobj( figHandle, 'Tag', 'imgFile'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'imgFileDesc'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'matFilePrefix'), 'String', '');
else
   set( findobj( figHandle, 'Tag', 'imgFile'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'imgFileDesc'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'matFilePrefix'), 'String', 'ra');
end   

if file.imaType ~=iT
   set( findobj( figHandle, 'Tag', 'imaHeaderInfo'), 'String', '');
   set( findobj( figHandle, 'Tag', 'Axes1Image1'), 'Visible', 'off');
   file.filename=0;
   set( findobj( figHandle, 'Tag', 'action'), 'visible', 'off');
   set( findobj( figHandle, 'Tag', 'selectImaFiles'), 'UserData', file);
end
