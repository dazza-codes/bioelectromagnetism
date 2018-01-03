function where=maxima(f)
% where=maxima(f)
% find local maxima of vector f
% for minima use maxima(-f)

% rgpm 10-2-2000

dy = diff(f);
where = find(([dy;0]<0) & ([0;dy]>0)); %indices of maxima
% over