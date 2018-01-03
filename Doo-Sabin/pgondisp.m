function varargout=pgondisp(pt,pgon)
% PGONDISP Display polyhedral mesh.
%
% [LI,PAT]=PGONDISP(PT,PGON) displays the polyhedral mesh
% defined by PT and PGON and returns a light handle LI and
% an array of patch handles PAT.
% 
% PT should be an m-by-3 matrix containing the coordinates
% of the points, PGON should be a cell array of polygons, each
% polygon being a pointlist(of indices to points).
%
% Example(Dice):
%
% pt=[0 0 0;1 0 0;1 1 0;0 1 0;...
%        0 0 1;1 0 1;1 1 1;0 1 1];
% pgon={[1 2 3 4] [5 8 7 6] [1 5 6 2]...
%             [3 7 8 4] [2 6 7 3] [4 8 5 1]};
% pgondisp(pt,pgon);
%
% See also DOO PGONTRACE PGONORIENT CHSURF


%************************* INPUT ARGUMENT CHECK ***********************
%**********************************************************************

error(nargchk(2,2,nargin));

if size(pt,2)~=3
   pt=pt';
end
   
if size(pt,2)~=3
   error(' pt should be an m-by-3 matrix');
end

%**************************** GRAPHICS ********************************
%**********************************************************************

axis off;
view(17,14);
rotate3d on;

pat=zeros(1,length(pgon));
for i=1:length(pgon)
   ptpgoni=pt(pgon{i},:);    % speed
   pat(i)=patch(ptpgoni(:,1) ,  ptpgoni(:,2) , ptpgoni(:,3) ,[1.0000    0.7812    0.4975],...
      'facelighting','flat','edgecolor',[0 0 0],...
      'edgelighting','flat','backfacelighting','lit');
end


li=light;
material dull;

%************************* OUTPUT ARGUMENT CHECK ***********************
%***********************************************************************

switch nargout
case 1
   varargout{1}=li;
case 2
   varargout{1}=li;
   varargout{2}=pat;
end

 
   
