function X_zero=assignPerimeterPositions(perimeterDists); 
% X_zero=assignPerimeterPositions(unfoldMesh.orderedUniquePerimeterPoints,perimeterDists); 
% Set distances around circle to match actual distances from the center.
% We should use something like this instead of just forcing the perimeter points to lie on a circle because the floodfill 
% perimeter finding routine only guarantees that all the points will be >at least< 'radius' distance away. Some may be significantly more..
numPerimPoints=length(perimeterDists);

radius=1000;% What does this mean really?

angles=linspace(0,2*pi,numPerimPoints)';
X_zero=zeros(numPerimPoints,2);

% Perimeter must be convex for Tule's algorithm  to work properly - 
% Removing the next line makes a more interesting-looking flat map but
% it can introduce crossing errors in the perimeter.
%perimeterDists=ones(size(perimeterDists))*mean(abs(perimeterDists(:)));

X_zero(:,1)=perimeterDists'.*cos(angles);
X_zero(:,2)=perimeterDists'.*sin(angles);
