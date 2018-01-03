function opgon=pgonorient(argpgon)
% PGONORIENT orients a closed surface of polygons
%
% Let pgon be a cell array of polygons, each polygon
% being a pointlist(of indices to points), then
% orpgon=PGONORIENT(pgon) orients the polygons defined by
% pgon with respect to each other.
% More accurately, it reverses the order of some polygons
% such that in the outgoing surface each edge is run once
% in each direction. SURFACE MUST BE CLOSED AND ORIENTABLE.
% The first polygon is assumed to have the right orientation.
%
% Example(Tetraeder):
%
%   >>pgon={[1 2 3] [1 2 4] [2 3 4] [1 3 4]};
%   >>newpgon=pgonorient(pgon);  
%
% See also DOO  PGONDISP  PGONTRACE  CHSURF

global opgon opgonlen onrpts or nrorpgons orientation
%***** To keep the recursive function orient simple, globals are used


opgon=argpgon;                  % pass the argument to global
opgonlen=length(opgon);     


maxpt=zeros(opgonlen,1);        % speed
onrpts=zeros(opgonlen,1);       % speed
for i=1:opgonlen
   maxpt(i)=max(opgon{i});
   onrpts(i)=length(opgon{i});  % number of points in opgon{i}
end
ptlen=max(maxpt);               % number of points alltogether


pt=fullpointinfo(ptlen,opgon,opgonlen,onrpts);
%***** To each point, this function calculates the polygons
%***** containing it.

opgon=addfirstpoint(opgon);
%***** Adds the first point of the pgon to the
%***** end of the pgons's pointlist (convenient when referencing)

opgon=fulledgeinfo(pt,opgon,opgonlen,onrpts);
%***** To each polygon's edge, this function calculates the index 
%***** of the adjacent polygon.

or=zeros(1,opgonlen);  % or=1 for oriented pgons, zero for not yet oriented
or(1)=1;               % the first pgon is assumed oriented correctly
nrorpgons=1;           % hence the number of oriented pgons is one
orientation='pos';     % If not, orient will set it to 'none'

orient(1);             % this recursive function does all the work


for i=1:opgonlen
   opgon{i}=opgon{i}(1,1:onrpts(i));  % brings pgon into input form
end

if orientation(1)=='p'
   fprintf('\n Polyhedral mesh was already oriented. \n');
end



%**************************************************************************
%******************************* SUBROUTINES ******************************



function pt=fullpointinfo(ptlen,opgon,opgonlen,onrpts)
%***** To each point, this function calculates the polygons
%***** containing it.

pt=cell(ptlen,1);
for i=1:ptlen
   pt{i}=[];                 % reserve space 
end

for i=1:opgonlen
   for j=1:onrpts(i)      
      pgonij=opgon{i}(j);     % speed
      pt{pgonij}=[pt{pgonij} i];
   end                       % concatenate pgon-data i to
end                          % existing pt-data



%**************************************************************************



function newpgon=addfirstpoint(opgon);
%***** Adds the first point of the pgon to the
%***** end of the pgons's pointlist (convenient when referencing)

for i=1:length(opgon)
   pgoni=opgon{i};
   newpgon{i}=[pgoni pgoni(1)];
end



%**************************************************************************



function newpgon=fulledgeinfo(pt,opgon,opgonlen,onrpts)
%***** To each polygon's edge, this function calculates the index of the
%***** adjactent polygon.

%***** datastructure: newpgon{i} is a 2-by-onrpts+1 matrix, the first row   
%***** containing the points of which the polygon consists, the second the  
%***** indices of the adjactent polygons.

newpgon=cell(opgonlen,1);

for i=1:opgonlen
   pgoni=opgon{i};             % speed
   nrptsi=onrpts(i);           % speed
   row2=zeros(1,nrptsi);  % speed

   
   for j=1:nrptsi
         
         thispt=pgoni(j);         % speed,clearness
         nextpt=pgoni(j+1);       % speed,clearness

         adjpgon = intersection(pt{thispt},pt{nextpt},i); 
               
                     % the index of the adjactent polygon 
                     % can be determinded by intersecting the polygons containg
                     % the first edgepoint with those containig the second.
                     
         if isempty(adjpgon)
            fprintf('\n An error has occured at edge %i \n',j);               
            fprintf(' of polygon %i. Did not find an adjacent polygon. \n\n',i);
            error('Input Data Error');
         end
                     
         row2(j)=adjpgon;     % and set the value
         
   end 
   newpgon{i}=[pgoni ; [row2 row2(1)] ];
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



%**************************************************************************



function orient(pgonnr)
%***** Actual orientation. Algorithm is as follows: Starting with the first
%***** pgon we orient its neighbouring pgons, then the neighbouring
%***** pgons of these if not already oriented and so on.

global opgon opgonlen onrpts or nrorpgons orientation

notoriented=[];           % neighbouring pgons that are not yet oriented
pgonpgonnr=opgon{pgonnr};  % speed

for j=1:onrpts(pgonnr)
      
   adjpgon=pgonpgonnr(2,j);     % adjactent pgon
     
   if or(adjpgon)==0            % if adjactent pgon is not oriented
      notoriented=[notoriented adjpgon];  % put it in the list
      
      pgonadjpgon=opgon{adjpgon};      % speed
      thispt=pgonpgonnr(1,j);         % speed,clearness
      nextpt=pgonpgonnr(1,j+1);       % speed,clearness
      
      idxnextpt=findfirst(opgon{adjpgon}(1,:),nextpt);
      %***** the index of the next point in the adjactent polygon
            
      if pgonadjpgon(1,idxnextpt + 1) ~= thispt     % wrong orientation
         pgonadjpgon=fliplr(pgonadjpgon);           % reverse order
         pgonadjpgon(2,:)=[pgonadjpgon(2,2:onrpts(pgonnr)+1) 0];
                                                    % correct edge information
         opgon{adjpgon}=pgonadjpgon;                 % and store the data
         orientation='none';                                           
      end
      
      
      or(adjpgon)=1;                  % now the adjactent pgon is oriented
      nrorpgons=nrorpgons+1;          % and the number of oriented pgons 
                                      % increases by one.
      
   end
end

if nrorpgons ~=opgonlen                % if not already finished
   for j=notoriented                  % orient the neighbouring pgons
      orient(j);                      
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



