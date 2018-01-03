function ans = visionToSPM_SHOWimgFile( figHandle);

aT = get( findobj( figHandle, 'Tag', 'actionType'), 'Value');
iT = get( findobj( figHandle, 'Tag', 'imaType'), 'Value');
if (aT == 3 | aT == 4) | iT==2
   ans = 0;
else
   ans = 1;
end

