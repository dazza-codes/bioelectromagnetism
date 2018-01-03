function groupedNodeStruct=findGroups(mesh,nodeList)
% Takes a list of nodes and tries to group it into unconnected unique clusters
% For each cluster, it returns the number of faces and edges contained within that cluster

numVertices=length(mesh.connectionMatrix);
groupsFound=0;
currentGroup=1;

rowsToSearch=[1:numVertices];
foundAll=0;

searchConnectionMatrix=mesh.connectionMatrix(nodeList,nodeList);

counter=0;

while (~foundAll)
   nodesFoundThisTime=1;
   allNodesFound=[];
   [startNode x]=find(searchConnectionMatrix);
   currentNodes=startNode(1); % First non-zero row is the initial node
  
   
   while(nodesFoundThisTime)
      [y connectedNodes]=find(searchConnectionMatrix(currentNodes,:));
      connectedNodes=unique(connectedNodes(:));  % All the things connected to current node set
      nodesFoundThisTime=length(connectedNodes);
      
      if (nodesFoundThisTime)
         searchConnectionMatrix(:,currentNodes)=0;
         allNodesFound=[allNodesFound;connectedNodes];
         currentNodes=connectedNodes;
         
      end % if nothing found
     
      
   end % Loop while you're finding something
   groupedNodeStruct{currentGroup}.nodeList=unique(nodeList(allNodesFound));
   groupedNodeStruct{currentGroup}.tempList=unique(allNodesFound);
   currentGroup=currentGroup+1;
   
   foundAll=~(sum(searchConnectionMatrix(:)));
   
   counter=counter+1;
   disp(counter);
   

end

fprintf('%d groups found',counter);

% Re-generate this...
searchConnectionMatrix=mesh.connectionMatrix(nodeList,:);
searchConnectionMatrix=searchConnectionMatrix(:,nodeList);
searchConnectionMatrix=searchConnectionMatrix;

%for t=1:(counter)
%   groupNodes=groupedNodeStruct{t}.nodeList;
 %  tempNodes=groupedNodeStruct{t}.tempList;
   
   % Need to know V, E, F - number of vertices, edges and faces in this group
   
   % Edges are easy...
 %  groupedNodeStruct{t}.E=sum(sum(searchConnectionMatrix(:,tempNodes)))/2;
   
   % How many individual triangles are entirely composed of these nodes?
   % Have identified good vertices. Now make a list of good faces - faces that contain at least 2 good vertices
   % This bit taken from locateStartNode
 %  ufl=mesh.uniqueFaceIndexList(:,1);
 %  ufi(:,1)=ismember(ufl,groupNodes);
 %  ufl=mesh.uniqueFaceIndexList(:,2);
 %  ufi(:,2)=ismember(ufl,groupNodes);
 %  ufl=mesh.uniqueFaceIndexList(:,3);
 %  ufi(:,3)=ismember(ufl,groupNodes);


  %  goodFaces=squeeze(find(sum(ufi')==3)); %All three vertices of the face must be in the outside set

 %   groupedNodeStruct{t}.F=length(unique(goodFaces));
 %   groupedNodeStruct{t}.V=length(unique(tempNodes));
  %  V=groupedNodeStruct{t}.V;
  %  E=groupedNodeStruct{t}.E;
  %  F=groupedNodeStruct{t}.F;
   
   
  %  eulerNum=(F-E+V);
    
  %  groupedNodeStruct{t}.eulerNum=eulerNum;
    
    % Think this is still wrong somehow. How?
   
    
% end
 



      