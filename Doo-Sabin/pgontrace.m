function varargout=pgontrace(pt,pgon,colorfile)
% PGONTRACE Raytraces a polyhedral mesh
%
% [LI,PAT,IVN]=PGONTRACE(PT,PGON[,COLORFILE]) raytraces the
% polyhedral mesh defined by PT and PGON using Vertexnormal-
% Interpolation and returns an array of light handles LI, an
% array of patch handles PAT and the interpolated 
% vertexnormals IVN.
%
% PT should be an m-by-3 matrix containing the coordinates
% of the points, PGON should be a cell array of polygons, each
% polygon being a pointlist (of indices to points).
% COLORFILE is the name of an m-file that defines the color of
% a point on the surface in dependency of its coordinates 
% and its normalvector. See COLORFILE for details.
%
% See also COLORFILE DOO PGONDISP


%************************* INPUT ARGUMENT CHECK ***********************
%**********************************************************************

error(nargchk(2,3,nargin));

if size(pt,2)~=3
   pt=pt';
end
   
if size(pt,2)~=3
   error(' pt should be a m-by-3 matrix');
end

if nargin==2
   colorfile='defaultcolormethod';
end


%**************************** GRAPHICS ********************************
%**********************************************************************

map=get(gcf,'colormap');      % set finer colormap
newmap=map;
if size(map,1)<=64
   for i=1:size(map,1)-1
      newmap(4*i-3:4*i,:)=[map(i,:); 0.75*map(i,:)+0.25*map(i+1,:);...
                                     0.50*map(i,:)+0.50*map(i+1,:);...
                                     0.25*map(i,:)+0.75*map(i+1,:)];
   end
end
set(gcf,'colormap',newmap);            


cla;
axis off;                                % some graphic settings
set(gca,'Projection','perspective');
view(17,14);
rotate3d on;
li=light;


pat=zeros(length(pgon),1);               % patch handles
vn =cell(1,length(pgon));                % vertex normals

for i=1:length(pgon)
   
   pgoni=pgon{i};                        % speed
   ptpgoni=pt(pgoni,:);                  % speed
   pat(i)=patch(ptpgoni(:,1) ,  ptpgoni(:,2) , ptpgoni(:,3) ,[0 0 0],...%sum(ptpgoni')',...
               'FaceLighting','Phong','EdgeColor','None',...
               'EdgeLighting','None','Backfacelighting','unlit');
            
   vni=get(pat(i),'VertexNormal');       % get Vertexnormals
            
   for j=1:size(pgoni,2)
      vnij=vni(j,:);                 ,   % speed
      vn{i}(j,:)=vnij/norm(vnij);        % normalize Vertexnormals
   end
end

ptdata=fullpointinfo(pt,pgon);

ivn=zeros(length(pt),3);                 % interpolated Vertexnormal 
colsize=length(feval(colorfile,[0 0 0],[1 1 1]));
col=zeros(length(pt),colsize);           % color of pt

for i=1:length(pt)
   ptdatai=ptdata{i};                            % speed
   ivni=vn{ptdatai(1,1)}(ptdatai(2,1),:);        % first vertexnormal
   for j=2:size(ptdatai,2)
      ivni=ivni+vn{ptdatai(1,j)}(ptdatai(2,j),:);    % add Vertexnormals belonging
   end                                               % to the same point
   ivni=ivni/norm(ivni);                             % normalize,
   ivn(i,:)=ivni;                                    % and store them
   
   col(i,:)=feval(colorfile,pt(i,:),ivni);           % get color of pt
end

for i=1:length(pgon)
   set(pat(i),'VertexNormal',ivn(pgon{i}(:),:),...
              'FaceVertexCdata',col(pgon{i}(:),:),...
              'FaceColor','interp');
end

   

%************************* OUTPUT ARGUMENT CHECK ***********************
%***********************************************************************


switch nargout
case 1
   varargout={li};
case 2
   varargout={li pat};
case 3
   varargout={li pat ivn};
end


%**************************** SUBROUTINES ******************************
%***********************************************************************

   
function newpt=fullpointinfo(pt,pgon)
%***** To each point, this function calculates the polygons
%***** containing it, and the indices of the point in these polygons.

newpt=cell(length(pt),1);

for i=1:length(pgon)
   pgoni=pgon{i};            % speed
   for j=1:size(pgoni,2)      
      pgonij=pgoni(j);       % speed
      newpt{pgonij}=[newpt{pgonij} [i;j]];
   end                       % concatenate pgon-data([i,j]) to
end                          % existing pt-data


function color=defaultcolormethod(pt,nv)
% ***** If no colorfile is specified, defaultcolormethod is used.
color=sum(nv);




