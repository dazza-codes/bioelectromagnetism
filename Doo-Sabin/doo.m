function varargout=qdoo(argpt,argpgon,n,varargin)
% DOO Doo-Sabin subdivision algorithm
%
% Given a polyhedral mesh defined by points PT and polygons PGON,
%
% [PTOUT,PGONOUT]=DOO(PT,PGON,N,PARAM1,VALUE1,PARAM2,...)
%
% applies the Doo-Sabin subdivision algorithm N times.
% MESH MUST BE CLOSED AND ORIENTABLE.
% For almost all inputdata (wrt to Lebesgue-measure), the sequence
% of meshes converges to a C^1 surface. (If for example all points
% of a polygon have the same coordinates (a pathological polygon),
% you may end up with a cusp.)
%
% PT should be an m-by-3 matrix containing the coordinates
% of the points, PGON should be a cell array of polygons, each
% polygon being a pointlist(of indices to points).
%
% Parameters:
%    ORIENTATION [Positive|Negative|{None}]
%    OUTPUT            [Full |{Simple}]
%    STATS             [{On} | Off]
% 
% Set ORIENTATION to 'Positive' if each edge of the
% mesh is passed once in each direction, and polygonpoints
% rotate counterclockwise with respect to the normalvector.
% Set ORIENTATION to 'Negative' if points rotate clockwise.
%
% If OUTPUT is set to 'Simple' (the default) the output
% has the same form as the inputs PT and PGON.
% If OUTPUT is set to 'Full' PTOUT is an m-by-5 matrix, one
% row per point, the first three entries being the coordinates,
% the fourth the number of a polygon that the point belongs to,
% and the fifth the index of the point in that polygon.
% OUTPGON{i} is a 2-by-nrpts+1 matrix, the first row 
% containing the points of which the polygon consists with the
% first point listed again at the end, the second the indices 
% of the adjactent polygons.
%
% If STATS is 'On' some calculating statistics are displayed.
% 
%
% Example(Dice):
%
% pt=[0 0 0;1 0 0;1 1 0;0 1 0;...
%        0 0 1;1 0 1;1 1 1;0 1 1];
% pgon={[1 2 3 4] [5 8 7 6] [1 5 6 2]...
%             [3 7 8 4] [2 6 7 3] [4 8 5 1]};
% doo(pt,pgon,3);
%
% See also PGONDISP PGONTRACE ORIENT CHSURF




global  pt ptlen pgon pgonlen nonedgepgonlen nrpts sumnrpts newpt newpgon w;


%******************** argument check ******************************************
%******************************************************************************

error(nargchk(3,7,nargin));
if mod(nargin,2)==0
   error('wrong number of input arguments');
end


if size(argpt,2)~=3
   pt=argpt';
else
   pt=argpt;
end

if size(pt,2)~=3
   error(' pt should be a m-by-3 matrix ');
end

pgon=reshape(argpgon,length(argpgon),1);

orientation='none';
output='simple';
checks='on';
stats='on';

orientation=0;
or=zeros(length(pgon),1);
or(1)=1;

for i=1:2:length(varargin)
   switch lower(varargin{i}(1:2))    
   case 'or'  % orientation
      switch lower(varargin{i+1}(1:2))
      case 'po' % positive
         orientation=1;
      case 'ne' % negative
         for j=1:length(pgon)
            pgoni=pgon{i};
            pgon{i}=pgoni(length(pgoni):-1:1);
         end
         orientation=1;
      case 'no' % none
      otherwise
         error('Bad value for property ''orientation''.');
      end     
      
   case 'ou'  % output
      switch lower(varargin{i+1}(1))
      case 'f' % full
         output='full';
      case 's' % simple
      otherwise
         error('Bad value for property ''output''.');
      end     
      
   case 'st' % stats
      switch lower(varargin{i+1}(1:2))
      case 'on'
      case 'of'
         stats='off';
      otherwise
         error('Bad value for property ''stats''.');
      end  
      
   otherwise
      error(['Invalid property : ''' varargin{i} '''']);
   end
end

      






%******************* get additional information about inputdata ***************
%******************************************************************************

tic;

if orientation==0
   t0=clock;
   
   pgon=pgonorient(pgon); % orients the polygons   
   
   if stats(1:2)=='on'
      ortime=etime(clock,t0);
      fprintf('\n Orientation  :  %.2fs ',ortime);
   end   
end

initpgonglobals(pgon);

pgon=order(pgon);  
%***** orders the polygonlist such that it starts with the extraordinary 
%***** pgons (=pgons with n~=4 vertices) 

pt=fullpointinfo(pt,pgon);
%***** To each point, this function calculates the polygons
%***** containing it, and the indices of the point in these polygons.

initpointglobals(pt);


pgon=addfirstpoint(pgon);
%***** Adds the first point of the pgon to the
%***** end of the pgons's pointlist (convenient when referencing)

pgon=fulledgeinfo(pt,pgon);
%***** To each polygon's edge, this function calculates the index of the
%***** adjactent polygon, and the index of the edge in this polygon.
%***** POLYGONS MUST BE ORIENTED NEGATIVELY


w=costable(pgon);     
%***** calculates the weights for convex-combination (doo-weights)





%************************ first doo step *************************************
%***** It is useful to treat this step separately, because the polygon-point
%***** data looks much simpler after its completion.(e.g. all vertices 
%***** connect exactly 4 edges)
%*****************************************************************************



[planeidx,vertidx,edgeidx,ptidx]=linearindexing(pt,pgon);
%***** Defines linear indices of both points and pgons of the next doo-step.
%***** planeidx contains indices of planepolygons, edgeidx->edgepolygons
%***** vertidx->vertexpolygons

[newpt,newpgon]=planepolygons(pt,pgon,planeidx,edgeidx);
%***** Calculates the new points of a doo-sabin step, and the polygons, that
%***** derive from the polygons of the previous step by shrinking.

[newpt,newpgon]=vertexpolygons(pt,pgon,newpt,newpgon,vertidx,edgeidx,ptidx);
%***** Calculates the polygons that derive from cutting off the vertices.

newpgon=edgepolygons(newpt,newpgon,pgon,edgeidx,vertidx,ptidx);
%***** Calculates the polygons that derive from cutting off the edges.



%****************** get some information about new data ********************
%***************************************************************************





redefineglobals;

nrpts=zeros(pgonlen,1);               
for i=1:pgonlen
   nrpts(i)=size(pgon{i},2)-1;
end
sumnrpts=sum(nrpts);

%***** this overrides some settings of redefineglobals.
%***** we merely need it after the first doo-step.






%************************ rest of the doo steps ****************************
%***** algorithm is the same in principal, but takes advantage of the
%***** improved knowledge about the data and the improved data structure.
%***** For greater speed we use global variables, but the data dependence
%***** is exactly the same as in the non-quick-routines. (the q stands for
%***** quick)
%***************************************************************************



if stats(1:2)=='on'
   tnew=toc;
   fprintf('\n DooStep 1  :  %.2fs \n',tnew);
end


for i=1:n-1
   qlinearindexing;  
   qplanepolygons;   
   qvertexpolygons;
   qedgepolygons;
   redefineglobals;
   if stats(1:2)=='on'
      told=tnew;
      tnew=toc;
      fprintf(' DooStep %i  :  %.2fs \n',i+1,tnew-told);
   end   
end

if stats(1:2)=='on'
   tnew=toc;
   fprintf('\n Total time    :  %.2fs \n',tnew);
   fprintf('\n Polygons : %i \n',pgonlen);   
   fprintf(' Points      : %i \n\n',ptlen);
end



%************************ output argument check *****************************
%****************************************************************************

switch nargout
case 0
   pt=pt(:,1:3);
   for i=1:pgonlen
     pgon{i}=pgon{i}(1,1:nrpts(i));
   end
   pgondisp(pt,pgon);
case 2
   if output(1)=='s'  % simple
      pt=pt(:,1:3);
      for i=1:pgonlen
         pgon{i}=pgon{i}(1,1:nrpts(i));
      end
      varargout={pt pgon};
   else               % full
      for i=1:nonedgepgonlen
         pgon{i}=pgon{i}(1:2,:);
      end
      varargout={pt pgon};
   end   
otherwise
   error(' Wrong number of output arguments');
end






%*************************** subroutines ************************************
%****************************************************************************



function initpgonglobals(pgon);
%***** we initialize our polygon-globals
global pgonlen nrpts sumnrpts

pgonlen=length(pgon);                 
nrpts=zeros(pgonlen,1);               
for i=1:pgonlen
   nrpts(i)=size(pgon{i},2);
end
sumnrpts=sum(nrpts);








function redefineglobals
%***** This converts the output of a doo-step to input for the next
global  pt ptlen pgon pgonlen nrpts newpt newpgon nonedgepgonlen sumnrpts  

pt=newpt;pgon=newpgon;

nonedgepgonlen=pgonlen+ptlen;
%***** polygons that do not derive from edges 

pgonlen=pgonlen+ptlen+sumnrpts/2;
%***** #planepgons + #vertpgons + #edgepgons 

oldptlen=ptlen;
ptlen=sumnrpts;
nrpts=[nrpts ; 4*ones(oldptlen+sumnrpts/2,1)];
%***** new polygons have exactly 4 vertices

sumnrpts=sumnrpts+4*(oldptlen+sumnrpts/2);








function newpgon=addfirstpoint(pgon);
%***** Adds the first point of the pgon to the
%***** end of the pgons's pointlist (convenient when referencing)

for i=1:length(pgon)
   pgoni=pgon{i};
   newpgon{i}=[pgoni pgoni(1)];
end






function newpgon=order(pgon)
%***** This function takes the extraordinary pgons to the
%***** beginning of the pgonlist.

global pgonlen nrpts 

expgonidx=find(nrpts~=4);
expgon=pgon(expgonidx);

expgonlen=length(expgon);
orpgonlen=pgonlen-expgonlen;

pgon(expgonidx)=[];
newpgon=[expgon ; pgon];

nrpts =[nrpts(expgonidx); 4*ones(orpgonlen,1)];







function newpt=fullpointinfo(pt,pgon)
%***** To each point, this function calculates the polygons
%***** containing it, and the indices of the point in these polygons.

%***** Data format: coos = 1-by-3 vector of point's coordinates
%*****              pgon = a 2-by-nrpts matrix, the first row is an
%*****                     indexvector to the polygons containing this point,
%*****                     the second indicates the position of the point in 
%*****                     these polygons

global pgonlen nrpts ptlen

newpt=cell(length(pt),1);
for i=1:length(pt)
   newpt{i}.coos=pt(i,:);  % convert points coordinates to new
   newpt{i}.pgon=[];       % data format  
end

for i=1:pgonlen
   pgoni=pgon{i};
   for j=1:nrpts(i)      
      newpt{pgoni(j)}.pgon=[newpt{pgoni(j)}.pgon [i;j]];
   end                       % concatenate pgon-data([i,j]) to
end                          % existing pt.data







function initpointglobals(pt)
%***** we initialize our point-globals

global ptlen nrpgons expgonlen

ptlen=length(pt);           
nrpgons=zeros(ptlen,1);       
for i=1:ptlen
   nrpgons(i)=size(pt{i}.pgon,2);
end   

troubles=find(nrpgons<3);
if ~isempty(troubles)
   fprintf('\n An error has occured. The points : [%s]\n',...
            int2str(troubles));               
   fprintf(' are contained in less than 3 polygons \n\n');
   error('Input Data Error');
end






function newpgon=fulledgeinfo(pt,pgon,orientation,or)
%***** To each polygon's edge, this function calculates the index of the
%***** adjactent polygon, and the index of the edge in this polygon.
%***** POLYGONS MUST BE ORIENTED NEGATIVELY

%***** datastructure: newpgon{i} is a 3-by-nrpts+1 matrix, the first row 
%***** containing the points of which the polygon consists, the second the 
%***** indices of the adjactent polygons, the third the indices of the edges
%***** in the adjactent polygons. Edges are 1-2,2-3,...,end-1.

global pgonlen ptlen nrpts nrpgons

for i=1:pgonlen
   newpgon{i}=zeros(3,nrpts(i)+1);   % reserve space
   newpgon{i}(1,:)=pgon{i};          % and convert the polygondata       

end


for i=1:pgonlen
   pgoni=pgon{i};             % speed
   
   for j=1:nrpts(i)
      if newpgon{i}(2,j)==0       % this pgon has not already 
                                  % been taken into account
         
         thispt=pgoni(j);         % clearness
         nextpt=pgoni(j+1);       % clearness

         adjpgon = intersection(pt{thispt}.pgon(1,:),pt{nextpt}.pgon(1,:),i); 
               
                     % the index of the adjactent polygon 
                     % can be determinded by intersecting the polygons containg
                     % the first edgepoint with those containig the second.
                     
                     
         if isempty(adjpgon)
            fprintf('\n An error has occured at edge %i \n',j);               
            fprintf(' of polygon %i. I did not find an adjactent polygon. \n',i);
            error('Input Data Error');
         end
  
         adjedge=findfirst(pgon{adjpgon},nextpt);  
                                                                        
                     % now get the index of the edge in the calculated pgon
                     % here we need the orintation of the polygons
                     
                     
         newpgon{i}(2:3,j)=[adjpgon adjedge];     % and set the values
         newpgon{adjpgon}(2,adjedge)=i;   
         
                     % we do not need the other indices of the adjactent pgon
      end            % any more, so we leave them zero
   end 
end

function val=intersection(a,b,n)    
%***** if a and b contain exactly two equal
%***** values, and one of them is given as n,
%***** this function calculates the other.

val=[];
for i=1:length(a)
   if a(i)~=n
      f=b(find(b==a(i)));     
      if ~isempty(f)
         val=f(1);
         break;
      end
   end   
end

function first=findfirst(array,number)
%***** If array is known to contain number, this function returns
%***** the index of the first appearance of number in array.

i=1;
while 1
   if array(i)==number
      first=i;
      break;
   end
   i=i+1;
end










function out=costable(pgons)
%***** calculates the weights for convex-combination of the polygonpoints

global nrpts nrpgons

for n=3:max([nrpts ; nrpgons ; 4]) % this number does not increase along iterations
   out{n}=zeros(n,n);
   w(1)=1/4*(1+5/n);                      
   for j=2:n
      w(j)=(3+2*cos(2*pi*(j-1)/n))/(4*n); 
   end                         % calculate the doo-weights for first point
   out{n}(1,:)=w;
   
   for k=2:n                   % and now rotate them to get the doo-weights
      w=[w(n) w(1:n-1)];       % for the other points
      out{n}(k,:)=w;
   end
end







function [planeidx,vertidx,edgeidx,ptidx]=linearindexing(pt,pgon);
%***** Defines linear indices of both points and pgons of the next doo-step.

global pgonlen nrpgons nrpts 

pgons4=(nrpgons==4);          % we determine which points are ordinary

orptidx=find(pgons4);         % and get their indices
orptlen=length(orptidx);      

exptidx=find(~pgons4);        % the same for extraordinary points
exptlen=length(exptidx);      


expgonlen=length(find(nrpts~=4));  % number of extraordinary polygons
orpgonlen=pgonlen-expgonlen;

lxlow=1;                  % linear index
lxhigh=expgonlen+1;       % the first pgons derive from extraordinary planepgons.
planeidx=lxlow:lxhigh-1;  % they are extraordinary themselves


lxlow=lxhigh;             % then we reserve space for the extraordinary
lxhigh=lxlow+exptlen;     % pgons that derive from extraordinary vertices.
vertidx(exptidx)=lxlow:lxhigh-1;

lxlow=lxhigh;             % all extraordinary pgons have been taken into account,
lxhigh=lxlow+orpgonlen;   % next polygons are those from ordinary pgons
planeidx=[planeidx lxlow:lxhigh-1];


lxlow=lxhigh;             % polygons from ordinary vertices
lxhigh=lxlow+orptlen;
vertidx(orptidx)=lxlow:lxhigh-1;

lx=lxhigh;                % polygons that derive from edges are last

edgeidx=cell(pgonlen,1);
for i=1:pgonlen
   edgeidx{i}=zeros(2,nrpts(i)+1);
end

for i=1:pgonlen
   pgoni=pgon{i};         % speed
   pgoni2=pgoni(2,:);
   pgoni3=pgoni(3,:);

   for j=1:nrpts(i)
      if edgeidx{i}(1,j)==0  % this pgon has not already been taken into account
         
         edgeidx{i}(:,j)=[lx ; 4];  % this edgepgon derives from the                           
                                    % j-th edge of the pgon i
         
                                    % the edgepolygon starts at the current
                                    % point, and therefore its fourth edge
                                    % will be adjactent to the shrunken pgon i

         edgeidx{pgoni2(j)}(:,pgoni3(j))=[lx 2];
                          % the appropriate settings for the adjactent pgon
         lx=lx+1;                            
      end
   end
end

ptlx=1;                 % while the linearindexing of pgons looks rather
for i=1:pgonlen         % complicated, linearindexing of points is very simple
   ptlxhigh=ptlx+nrpts(i)-1;
   ptidx{i}=[ptlx:ptlxhigh ptlx];    
   ptlx=ptlxhigh+1;                  
end








function [newpt,newpgon]=planepolygons(pt,pgon,planeidx,edgeidx)
%***** Calculates the new points of a doo-sabin step, and the polygons, that
%***** derive from the polygons of the previous step by shrinking.

global pgonlen ptlen nrpts sumnrpts nrpgons w

newpgon=cell(pgonlen+ptlen+sumnrpts/2,1);% #planepgons + #vertpgons + #edgepgons 

newpt=ones(sumnrpts,5); % note that we change the data-structure of the points.
                        % this is doo to the fact, that all new points are
                        % ordinary (connect 4 edges). The first three entries
                        % are the points coordinates, the fourth is a pgon
                        % that contains the point (no matter which one) and
                        % the fifth is the index of the point in that pgon
                        
lx=1;                   % linear index
for i=1:pgonlen
   nrptsi=nrpts(i);     % speed
   lxhigh=lx+nrptsi-1;
   
   newpgon{planeidx(i)}=[[lx:lxhigh lx] ; edgeidx{i}];     
                        % the shrunken polygon is found easily
                        
                                                
   p=zeros(nrptsi,3);    % now for the new points:
   for j=1:nrptsi        % we gather the pgon's point's coordinates in p                                
      p(j,:)=pt{pgon{i}(1,j)}.coos;   
   end
   
   newpt(lx:lxhigh,:)=[w{nrptsi}*p   planeidx(i)*ones(nrptsi,1)   [1:nrptsi]'];
                         % and convex-combine them using the weights w
   lx=lxhigh+1;   
end








function [newpt,newpgon]=vertexpolygons(pt,pgon,newpt,newpgon,vertidx,edgeidx,ptidx)
%***** Calculates the polygons that derive from cutting off the vertices.
%***** This algorithm may look a bit clumsy, but consider that we do not
%***** know the order in which to connect the points that form our new
%***** vertexpolygon. The points themselves are easily found.

global ptlen nrpts nrpgons

for i=1:ptlen                             
   nrpgonsi=nrpgons(i);                  % speed
   vertpgoni=zeros(3,nrpgonsi+1);        % speed                                      
   
   pgonnr  =pt{i}.pgon(1,1);    % pgonnr stands for a the index of a polygon 
   pgonptnr=pt{i}.pgon(2,1);    % containing the point i. We start arbitrarily 
   if pgonptnr==1               % with the first pgon in the pgonlist of i.
      edgenr=nrpts(pgonnr);     % pgonptnr is the index of the current point 
   else                         % in this pgon, and edgenr is the number of       
      edgenr=pgonptnr-1;        % the edge leading to this point.
   end
   
   
   for j=1:nrpgonsi-1            
      vertpgoni(1:3,j)=[ptidx{pgonnr}(pgonptnr); edgeidx{pgonnr}(1,edgenr);...
                        edgeidx{pgonnr}(2,edgenr)-1];   
                  
         % stores the information about the j-th vertex of the
         % vertexpolygon in vertpgoni. First are the points the
         % vertexpolygon consists of, second the adjactent pgons
         % and third the edgeindices.
         
      pgonnr=pgon{pgonnr}(2,edgenr);           
      pgonptnr=findfirst(pgon{pgonnr}(1,:),i); 
         
         % The new pgonnr is the number of the polygon adjactent 
         % to the old pgonnr along the edge edgenr. This polygon 
         % must contain our point i,so we search for its index. 
         % Doing this nrpgonsi times we encircle our vertex i and
         % find our vertexpolygon.
      
      if pgonptnr==1                               
         edgenr=nrpts(pgonnr);             
      else                                         
         edgenr=pgonptnr-1;
      end

   end    
  
   vertpgoni(1:3,nrpgonsi)=...
      [ptidx{pgonnr}(pgonptnr); edgeidx{pgonnr}(1,edgenr);...
       edgeidx{pgonnr}(2,edgenr)-1];   

   
   newpgon{vertidx(i)}=vertpgoni;   
         % now take all the info about the vertexpgon to newpgon
   newpgon{vertidx(i)}(1,nrpgonsi+1)=vertpgoni(1,1);
         % for conveniece's sake we add the first point as last again
end       
                                                      



                                                      

                                                      
                                                      
function newpgon=edgepolygons(newpt,newpgon,pgon,edgeidx,vertidx,ptidx)
%***** Calculates the polygons that derive from cutting off the edges.

global pgonlen nrpts sumnrpts

edgepgon=cell(sumnrpts/2,1);     % We know how many this will be

lx=edgeidx{1}(1,1);
for i=1:pgonlen
   pgoni=pgon{i};
   pgoni1=pgoni(1,:);
   pgoni2=pgoni(2,:);
   pgoni3=pgoni(3,:);
   ptidxi=ptidx{i};
   
   for j=1:nrpts(i)
      pgoni3j=pgoni3(j);       % speed
      
      if pgoni3j~=0            % an edge belongs to 2 pgons => no double counts,     
         pgoni2j=pgoni2(j);    % (we left pgon{i}(3,j) zero for pgons that
         ptidxpgoni2j=ptidx{pgoni2j};
                               % share edge j with a pgon of lower lx
                               % (in fulledgeinfo))
         newpgon{lx}...                                        
         =[ptidxi(j)                       ptidxpgoni2j(pgoni3j+1)...
           ptidxpgoni2j(pgoni3j)         ptidxi(j+1)             ptidxi(j);
           vertidx(pgoni1(j))              0 ...
           vertidx(pgoni1(j+1))            0               0];
     
                  % the last zero is dummy, the other ones replace 
                  % pgondata that will not be used in subsequent calculations.
         lx=lx+1;   
     end
   end 
end


function qlinearindexing;
%***** Defines linear indices of both points and pgons of the next doo-step.
%***** Quick version, taking advantage of improved datastructure

global  pt ptlen nrpts pgon pgonlen nonedgepgonlen vertidx edgeidx ptidx


lx=pgonlen+ptlen;      % we do not need planeidx anymore, since planeidx(i)=i
vertidx=pgonlen+1:lx;  % so we start with vertidx

lx=lx+1;      

edgeidx=cell(pgonlen,1);

for i=nonedgepgonlen+1:pgonlen      % then edgeidx
   edgeidx{i}=2*ones(2,nrpts(i)+1); % the last coordinate is dummy (concat fit)
end                                 % the edgeindices are known to
                                    % be equal to 2 (sketch).
                                    

for i=1:nonedgepgonlen
   nrptsi=nrpts(i);                 % speed
   pgoni=pgon{i};

   edgeidx{i}=[lx:lx+nrptsi ; 4*ones(1,nrptsi+1)];  
    
                                    % this edgepgon derives from the                           
                                    % j-th edge of the pgon i
                                    % the edgepolygon starts at the current
                                    % point, and therefore its fourth edge
                                    % coincides with the current.  
   for j=1:nrptsi               
      edgeidx{pgoni(2,j)}(1,pgoni(3,j))=lx;                                        
      lx=lx+1;          % the appropriate settings for the adjactent pgon                          
   end
end



ptlx=1;
for i=1:pgonlen                     % linearindexing of points stays the
   ptlxhigh=ptlx+nrpts(i)-1;        % same.
   ptidx{i}=[ptlx:ptlxhigh ptlx];      
   ptlx=ptlxhigh+1;                  
end






function qplanepolygons
%***** Calculates the new points of a doo-sabin step, and the polygons, that
%***** derive from the polygons of the previous step by shrinking.
%***** Quick version, taking advantage of improved datastructure

global pt ptlen pgon pgonlen nrpts newpt newpgon edgeidx sumnrpts w

newpgon=cell(pgonlen+ptlen+sumnrpts/2,1);%#planepgons + #vertpgons + #edgepgons 
newpt=ones(sumnrpts,5);

lx=1;                    % linear index
for i=1:pgonlen
   nrptsi=nrpts(i);  
   
   lxhigh=lx+nrptsi-1;
   
   newpgon{i}=[[lx:lxhigh lx] ; edgeidx{i}];      % planepolygon easily found
   
   p=pt(pgon{i}(1,1:nrptsi),1:3); % We gather the pgon's point's coordinates in p
   
   newpt(lx:lxhigh,:) = [w{nrptsi}*p i*ones(nrptsi,1) [1:nrptsi]'];
                                  % and convex-combine them to the new points
   lx=lxhigh+1;
   
end

function qvertexpolygons
%***** Calculates the polygons that derive from cutting off the vertices.
%***** This algorithm may look a bit clumsy, but consider that we do not
%***** know the order in which to connect the points that form our new
%***** vertexpolygon. The points themselves are easily found.
%***** Quick version, taking advantage of improved datastructure


global pt ptlen pgon nrpts newpgon edgeidx vertidx ptidx


for i=1:ptlen                          
     
   vertpgoni=zeros(3,5);                                             
   
   pgonnr  =pt(i,4);  % Due to the changed data format 
   pgonptnr=pt(i,5);  
   if pgonptnr==1                       
      edgenr=nrpts(pgonnr);           
   else                                  
      edgenr=pgonptnr-1;
   end

   
   for j=1:3                                      
      
      vertpgoni(1,j)=ptidx{pgonnr}(pgonptnr);        
      vertpgoni(2,j)=edgeidx{pgonnr}(1,edgenr);       
                                                    
                                                
                                                   
      pgonnr=pgon{pgonnr}(2,edgenr);                
      pgonptnr=findfirst(pgon{pgonnr}(1,:),i);      
      if pgonptnr==1                                     
         edgenr=nrpts(pgonnr);             
      else                                         
         edgenr=pgonptnr-1;
      end
   end   
         
   vertpgoni(1,4)=ptidx{pgonnr}(pgonptnr);  
   vertpgoni(2,4)=edgeidx{pgonnr}(1,edgenr);    

   
   vertpgoni(3,1:4)=[3 1 3 1];  % due to the ordering of the polygons,
                                % we know this in advance (sketch)
                                             
   newpgon{vertidx(i)}=vertpgoni;
   newpgon{vertidx(i)}(1,5)=vertpgoni(1,1);
end



function qedgepolygons
%***** Calculates the polygons that derive from cutting off the edges.
%***** Quick version, taking advantage of improved datastructure

global  pgon nonedgepgonlen nrpts newpgon vertidx edgeidx ptidx

lx=edgeidx{1}(1);
for i=1:nonedgepgonlen
   pgoni=pgon{i};
   pgoni1=pgoni(1,:);
   pgoni2=pgoni(2,:);
   pgoni3=pgoni(3,:);
   ptidxi=ptidx{i};

   for j=1:nrpts(i)
      pgoni2j=pgoni2(j);      
      pgoni3j=pgoni3(j);    
      ptidxpgoni2j=ptidx{pgoni2j};
      
      newpgon{lx}...                                        
         =[ptidxi(j)                    ptidxpgoni2j(pgoni3j+1)...
           ptidxpgoni2j(pgoni3j)        ptidxi(j+1)       ptidxi(j);...
           vertidx(pgoni1(j))           0 ... 
           vertidx(pgoni1(j+1))         0                 0];
     lx=lx+1;
   end 
end

