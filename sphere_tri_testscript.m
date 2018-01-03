% script to test the sphere triangulation
% order of vertices in faces is illustrated with RGB

% radius = 1, clockwise vertices in faces
FV = sphere_tri('ico',3,[],0);
% radius = 1, counterclockwise vertices in faces
%FV = sphere_tri('ico',3,[],1);

facecolor = [.7 .7 .7];
figure
Hp = patch('faces',FV.faces,'vertices',FV.vertices,...
    'facecolor',facecolor,'facealpha',0.8,'edgecolor',[.8 .8 .8]); 
camlight('headlight','infinite'); daspect([1 1 1]); axis vis3d; axis off
material dull; rotate3d
hold on

f = 10;
vertex_index1 = FV.faces(f,1);
vertex_index2 = FV.faces(f,2);
vertex_index3 = FV.faces(f,3);
vertex1 = FV.vertices(vertex_index1,:);
vertex2 = FV.vertices(vertex_index2,:);
vertex3 = FV.vertices(vertex_index3,:);

plot3(vertex1(1),vertex1(2),vertex1(3),'ro')
plot3(vertex2(1),vertex2(2),vertex2(3),'go')
plot3(vertex3(1),vertex3(2),vertex3(3),'bo')

vertNormals = get(Hp,'vertexnormals');

quiver3(vertex1(1),vertex1(2),vertex1(3),...
    vertNormals(vertex_index1,1),vertNormals(vertex_index1,2),vertNormals(vertex_index1,3));

view(3)
rotate3d
