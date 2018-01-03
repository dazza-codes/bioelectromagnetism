function cmap = grayish(cmap,minval);

% grayish - gray out the lower end of a spectrum for less severe palette
% function cmap = grayish(cmap,minval);
% assumes standard two-dimensional color map of three columns

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $


[m,n] = size(cmap);

if(length(minval) == 1),
  minval = minval + zeros(1,3);
end

for i = 1:3,
 cmap(:,i) = max(cmap(:,i),minval(i)); % set threshold
end

startndx = 0; % looking for highest index
for i = 1:3,
  startndx = max(startndx,min(find(cmap(:,1) > minval(i))));
end

% so color above this are above the minval

cmap = cmap(startndx:end,:); % only those values

% now rescale to full range by interpolating
mtrunc = size(cmap,1); % new length
ndx = [0:(mtrunc-1)]/(mtrunc-1)*(m-1); % scales to new range

cmap = interp1(ndx,cmap,[0:(m-1)]);

return
