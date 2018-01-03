% Generating random meshes to test mrFlatMesh
busyHandle=0;
showFigures=1;
statusHandle=0;
perimDist=9;
rand('state',0);

nPoints=1151;
x = (rand(1,nPoints)*1); 
y = (rand(1,nPoints)*1); 
z=rand(1,nPoints)*4;

TRI = delaunay(x,y);
testMesh.uniqueVertices=[x;y;z]';
testMesh.uniqueFaceIndexList=TRI;

[testMesh.connectionMatrix]=findConnectionMatrix(testMesh);

% testMesh.connectionMatrix(10,:)=0;
% testMesh.connectionMatrix(11,:)=0;
% 
% testMesh.connectionMatrix(:,10)=0;
% testMesh.connectionMatrix(:,11)=0;
% 

% Try mark and Josh's version of DIJKSTRA's algorithm
%function D = dijkstra( G , S )

% --------------------------------------------------------------------

%      Mark Steyvers, Stanford University, 12/19/00

% --------------------------------------------------------------------

%

% DIJKSTRA  Find shortest paths in graphs

% 	D = dijkstra( G , S ) use the full or sparse matrix G in which 

% 	an entry (i,j) represents the arc length between nodes i and j in a 

%	graph. In a full matrix, the value INF represents the absence of an arc; 

%	in a sparse matrix, no entry at (i,j) naturally represents no arc.

%		

%	S is the one-dimensional matrix of source nodes for which the shortest

%	to ALL other nodes in the graphs will be calculated. The output matrices

%	D and P contain the shortest distances and predecessor indices respectively.

%	An infinite distance is represented by INF. The predecessor indices contain

%	the node indices of the node along the shortest path before the destination

%	is reached. These indices are useful to construct the shortest path with the

%	function pred2path (by Michael G. Kay).
% Need to generate the distance 
[G]=sqrt(find3DNeighbourDists(testMesh));



% Find a point in the middle to use as a start node
medX=median(x(:));
startNode=582; %find(x==medX)

testMesh.dist=dijkstra(G,startNode);
% 
% [nodeList, edgeList]=generateManDistNodes(testMesh);
% dimdist =[1 1 1];
% testMesh.dist = mrManDist(nodeList,edgeList,startNode,dimdist,-1,0); 
testMesh.distMap=makeDistanceImage([x;y]',testMesh.dist,128);
	
colors=repmat(testMesh.dist,3,1)';
colors=colors/max(colors(:));

figure(1);
subplot(1,2,1);
hold off;
trisurf(TRI,x,y,z,testMesh.dist); 
colormap hot;

view(2),...
hold on; 
plot(x,y,'o');

set(gca,'box','on');
hold off;

figure(2);
subplot(2,1,1);
imagesc(testMesh.distMap);
axis image;
		
colormap hot;
title('Manifold distance map');
colorbar;
	
testMesh.heightMap=makeDistanceImage([x;y]',z,128);
subplot(2,1,2);
imagesc(testMesh.heightMap);
colorbar;
title('Z');
axis image;

mesh=testMesh;
mesh.vertices=mesh.uniqueVertices;
mesh.faceIndexList=TRI;

nVerts=length(mesh.vertices);

mesh.uniqueVertices=mesh.vertices;
mesh.vertsToUnique=[1:nVerts]';
mesh.UniqueToVerts=[1:nVerts]';

mesh.uniqueFaceIndexList=findUniqueFaceIndexList(mesh); % this gets rid of faces with multiple duplicate vertices and different permutations of the same face
mesh.uniqueFaceIndexList=TRI;


statusStringAdd(statusHandle,'Using threshold to find perimeter(s).');

% Find perims with simple thold
insideNodes=find(mesh.dist<=perimDist);
nVerts=length(insideNodes);

% Separate this new mesh straight away
newMesh.uniqueVertices=mesh.uniqueVertices(insideNodes);
newMesh.connectionMatrix=mesh.connectionMatrix(insideNodes,insideNodes);
newMesh.edgeList=findEdgesInGroup(newMesh,(1:nVerts));


niceVerts=mesh.uniqueVertices(insideNodes,:);
figure(15);
plot(niceVerts(:,1),niceVerts(:,2),'.');

insideNodes=insideNodes(:);
numBadNodes=9999999;
[perimeterEdges,eulerCondition]=findGroupPerimeter(mesh,insideNodes);

coordsInGroup=mesh.uniqueVertices(insideNodes);
facesInGroup=mesh.uniqueFaceIndexList(findFacesInGroup(mesh,insideNodes),:);
 
 
 figure(20);
 Z=ones(length(mesh.uniqueVertices),1);
 f=trimesh(facesInGroup,mesh.uniqueVertices(:,1),mesh.uniqueVertices(:,2),Z);
 
 insideNodes=removeHangingNodes(mesh,insideNodes);
 badPerimNodes=findBadPerimNodes(mesh,perimeterEdges);
 Z2=Z;
 
 Z2(badPerimNodes)=10;
 
  f2=trimesh(facesInGroup,mesh.uniqueVertices(:,1),mesh.uniqueVertices(:,2),Z2);
 

   
while (numBadNodes>0)
    [perimeterEdges,eulerCondition]=findGroupPerimeter(mesh,insideNodes);
	length(perimeterEdges)
	length(unique(perimeterEdges,'rows'))
    fprintf('\nEuler number=%d',eulerCondition);
    
    badPerimNodes=findBadPerimNodes(mesh,perimeterEdges);
    numBadNodes=length(badPerimNodes);
    fprintf('\nThere are %d bad perim nodes.',numBadNodes);
	
    
    if(numBadNodes)
		[insideNodes]=correctBadNodes(mesh,insideNodes,badPerimNodes);
        %[insideNodes]=setdiff(insideNodes,badPerimNodes);
		%insideNodes=removeHangingNodes(mesh,insideNodes);
    end
end

% The euler number will tell you whether you've got a mesh with no holes in it.

% The routine above can generate islands - get round it by zeroing the connection matrix for the largest perimeter
% and then doing a flood fill with no limits to generate the inside group

messageString=sprintf('Euler number for this set=%d',eulerCondition);
statusStringAdd(statusHandle,messageString);
uniquePerimPoints=unique(perimeterEdges(:));
messageString=sprintf('%d unique perimeter points found.',length(uniquePerimPoints));
statusStringAdd(statusHandle,messageString);
[orderedUniquePerimeterPoints,biggest]=orderMeshPerimeterPointsAll(mesh,perimeterEdges);
nPerims=size(orderedUniquePerimeterPoints);
fprintf('\n%d different perimeters found',nPerims);
orderedUniquePerimeterPoints=orderedUniquePerimeterPoints{biggest}.points;

tempConMat=mesh.connectionMatrix; % save it for later
mesh.connectionMatrix(orderedUniquePerimeterPoints,:)=0;
mesh.connectionMatrix(:,orderedUniquePerimeterPoints)=0;
[insideNodes,insideNodeStruct]=floodFillFindPerim(mesh,Inf,startNode,busyHandle);
insideNodes=[insideNodes(:);orderedUniquePerimeterPoints(:)];
mesh.connectionMatrix=tempConMat;
[perimeterEdges,eulerCondition]=findGroupPerimeter(mesh,insideNodes);
fprintf('\nEuler condition=%d',eulerCondition);


% internal points (the ones we want)
% ************************************************************************


 figure(21);
 [perimeterEdges,eulerCondition]=findGroupPerimeter(mesh,insideNodes);

coordsInGroup=mesh.uniqueVertices(insideNodes);
facesInGroup=mesh.uniqueFaceIndexList(findFacesInGroup(mesh,insideNodes),:);
 
 Z=ones(length(mesh.uniqueVertices),1);
 Z2=Z;
 Z2(badPerimNodes)=10;
f2=trimesh(facesInGroup,mesh.uniqueVertices(:,1),mesh.uniqueVertices(:,2),Z2);
unfoldMesh.connectionMatrix=mesh.connectionMatrix(insideNodes,insideNodes);
unfoldMesh.uniqueVertices=mesh.uniqueVertices(insideNodes,:);
unfoldMesh.dist=mesh.dist(insideNodes);


% Convert the edges to feed into orderMeshPerimeterPoints
fullEdgePointList=perimeterEdges(:);

[numEdges,x]=size(perimeterEdges);

newEdges=zeros((numEdges*2),1);
statusStringAdd(statusHandle,'Finding sub-mesh edges.');

for t=1:(numEdges*2)  
	if ((~mod(t,100)) & busyHandle)
		updateBusybar(busyHandle,t);
	end
	newEdges(t)=find(insideNodes==fullEdgePointList(t));
end

newEdges=reshape(newEdges,numEdges,2);
statusStringAdd(statusHandle,'Finding sub-mesh perim.');

% Find the perimeter points.
unfoldMesh.orderedUniquePerimeterPoints=zeros(length(orderedUniquePerimeterPoints),1);

for t=1:length(orderedUniquePerimeterPoints)
	f1=find(insideNodes==orderedUniquePerimeterPoints(t));%

	unfoldMesh.orderedUniquePerimeterPoints(t)=f1;%orderMeshPerimeterPoints(newEdges);

end


% Unfolding bit...
% Now we'd like to unfold this.
% Need the following things...
% Connection matrix N (nxn) - almost the same as the connection matrix except that rows=rows/sum(rows)
% and all points on the perimeter have been removed.
% X0 - a px2 matrix containing the 2D locations of the perimeter points.
% P - The perimeter connection matrix: (nxp) 'whose (i,j)th entry is 1/mi when perimeter node j is connected
% to sample node i. 

% Find the N and P connection matrices
statusStringAdd(statusHandle,'Finding sub-mesh con. mat.');
[N,P,unfoldMesh.internalNodes]=findNPConnection(unfoldMesh);
[unfoldMesh.distSQ]=find3DNeighbourDists(unfoldMesh);   % Here we find the 3D distance from each point to its neighbours.


% Assign the initial perimeter points - they're going to go in a circle for now...
numPerimPoints=length(unfoldMesh.orderedUniquePerimeterPoints);

perimeterDists=mesh.dist(orderedUniquePerimeterPoints); % Distance of each perimeter point from the start node

statusStringAdd(statusHandle,'Assigning perimeter points');

% We'd like to place the perimeter points down in an intelligent manner. We can place them at the correct distance from the start node
% and at the correct distance from each other. I think.
% We already know their distance from the start node, now we'd like to get their distances
% from each other.

% Should be able to extract this from unfoldMesh.distSQ

%  start at the first perimeter point in unfoldMesh.orderedUniquePerimeterPoints
% Find the distance between this and the next point, etc etc...
nPerimPoints=length(unfoldMesh.orderedUniquePerimeterPoints);
interPerimDists=zeros(nPerimPoints,1);

for thisPerimPoint=1:nPerimPoints
	nextIndex=mod((thisPerimPoint+1),nPerimPoints);
	if(nextIndex==0)
		nextIndex=1;
	end
	
	interPerimDists(thisPerimPoint)=unfoldMesh.distSQ(unfoldMesh.orderedUniquePerimeterPoints(thisPerimPoint),unfoldMesh.orderedUniquePerimeterPoints(nextIndex));
end

interPerimDists=sqrt(interPerimDists);


unfoldMesh.X_zero=assignPerimeterPositions(perimeterDists); % Can set distances around circle to match actual distances from the center. 
% Angular positions will come next.

statusStringAdd(statusHandle,'Solving position equation (slow)');
X=(speye(size(N)) - N) \ (sparse(P * unfoldMesh.X_zero));


unfoldMesh.N=N;
unfoldMesh.P=P;
unfoldMesh.X=X;

% Find out the differences between  
dist2DSQ=find2DNeighbourDists(unfoldMesh);
d1=sparse(unfoldMesh.distSQ)-sparse(dist2DSQ);
d1=abs(d1);


goodness=sum(sum(d1.^2));

messageString=sprintf('Current error per node: %d',full(sqrt(goodness))/length(insideNodes));
statusStringAdd(statusHandle,messageString);

% Show the mesh
if (showFigures)

	statusStringAdd(statusHandle,'Displaying unfold');
	figure(50);
	
	hold off;
	gplot(unfoldMesh.N,unfoldMesh.X);
	
	axis equal;
	axis off;
	zoom on

end

