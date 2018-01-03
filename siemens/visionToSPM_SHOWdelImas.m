function ans = visionToSPM_SHOWdelImas( figHandle);

aT = get( findobj( figHandle, 'Tag', 'actionType'), 'Value');

if (aT == 3 | aT == 4)
   ans = 0;
else
   ans = 1;
end
