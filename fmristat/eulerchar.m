function ec=eulerchar(input_file, input_thresh, mask_file, mask_thresh)

%EULERCHAR finds Euler characteristics for excursion sets.
%
% EC = EULERCHAR(INPUT_FILE, INPUT_THRESH [, MASK_FILE [, MASK_THRESH]])
%
% INPUT_THRESH: a vector of thresholds that defines the excursion sets 
% as INPUT_FILE >= INPUT_THRESH.
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
% If MASK_FILE is empty (default) the search volume is the whole volume.
%
% EC: vector of Euler characteristics for each threshold in INPUT_THRESH.
% The EC is defined as = #points - #edges + #triangles - #tetrahedra 
% inside the intersection of the excursion set and search volume.

ismask=(nargin >= 3);
if ismask
   if nargin < 4
      mask_thresh=[];
   end
   if isempty(mask_thresh)
      mask_thresh=fmri_mask_thresh(mask_file);
   end
   base=input_file(1:size(input_file,2)-4); 
   ext=input_file(size(input_file,2)+(-3:0)); 
   X_file=[base '_mesh' ext];
   m_file=[base '_mask' ext];
   delete(X_file);
   delete(m_file);
   mask_mesh(input_file,base,mask_file,mask_thresh,0);
   m_thresh=0.5;
else
   X_file=input_file;
end

d=fmris_read_image(X_file,0,0);
d.dim
numslices=d.dim(3);
J=d.dim(2);
I=d.dim(1);
IJ=I*J;

% Set up:

i=kron(ones(1,J),1:I);
j=kron(1:J,ones(1,I));

ex=find(i<I);
ex1=[ex; ex+IJ]';
ex2=[find(i>1); find(i>1)+IJ]';

ey=find(j<J);
ey1=[ey; ey+IJ]';
ey2=[find(j>1); find(j>1)+IJ]';

ez=1:(IJ);
ez1=ez;
ez2=ez+IJ;

exye=find((rem(i+j,2)==0)&(i<I)&(j<J));
exyo=find((rem(i+j,2)==1)&(i<I)&(j<J));
exy=[exye exyo];
exy1=[exye     exyo+1; exye+1+IJ exyo+IJ]';
exy2=[exye+1+I exyo+I; exye+I+IJ exyo+1+I+IJ]';

exze=find((rem(i+j,2)==0)&(i<I));
exzo=find((rem(i+j,2)==1)&(i<I));
exz =[exze exzo];
exz1=[exze exzo+1];
exz2=[exze+1+IJ exzo+IJ];

eyze=find((rem(i+j,2)==0)&(j<J));
eyzo=find((rem(i+j,2)==1)&(j<J));
eyz =[eyze eyzo];
eyz1=[eyze eyzo+I];
eyz2=[eyze+I+IJ eyzo+IJ];

% edges:

edge1=[ex1(:,1)' ey1(:,1)' exy1(:,1)' ...
       ez1 exz1 eyz1 ...
       ex1(:,2)' ey1(:,2)' exy1(:,2)'];
edge2=[ex2(:,1)' ey2(:,1)' exy2(:,1)' ...
       ez2 exz2 eyz2 ...
       ex2(:,2)' ey2(:,2)' exy2(:,2)'];
 
% triangles:

tri1=[exy1(:,1)' exy1(:,1)' exz1 exz1 eyz1 eyz1 exy1(:)' exy1(:)' exy1(:,2)' exy1(:,2)'];
tri2=[exy2(:,1)' exy2(:,1)' exz2 exz2 eyz2 eyz2 exy2(:)' exy2(:)' exy2(:,2)' exy2(:,2)'];
tri3=[exye+1 exyo exye+I exyo+1+I ...
    exze+1 exzo exze+IJ exzo+1+IJ ...
    eyze+I eyzo eyze+IJ eyzo+I+IJ ...
    exye+1+IJ exyo+IJ exye exyo+1 ...
    exye+I+IJ exyo+1+I+IJ exye+1+I exyo+I ...
    exye+IJ exyo+1+IJ exye+1+I+IJ exyo+I+IJ];

% tetrahedra:

tet1=[exye exyo         exye exyo+1+I+IJ      exye exyo+1+I+IJ ...
      exye exyo+1+I+IJ  exye+1+I exyo+1+I+IJ];
tet2=[exye+1 exyo+1     exye+I exyo+IJ        exye+1+IJ exyo+IJ ...
      exye+1+IJ exyo+IJ exye+1+IJ exyo+1 ];
tet3=[exye+1+I exyo+I   exye+1+I exyo+1+IJ    exye+I+IJ exyo+I+IJ ...
      exye+I+IJ exyo+1  exye+I+IJ exyo+I ];
tet4=[exye+1+IJ exyo+IJ exye+I+IJ exyo+1      exye+1+I exyo+I ...
      exye+IJ exyo+I    exye+1+I+IJ exyo+1+I ];

% START:

mask=zeros(IJ,2);
X=repmat(-Inf,IJ,2);
flip=2;
if ismask
   m=fmris_read_image(m_file,1,1);
   mask(:,flip)=reshape(m.data>m_thresh,IJ,1);
else
   mask(:,flip)=ones(IJ,1);
end
m=fmris_read_image(X_file,1,1);
X(:,flip)=reshape(m.data,IJ,1);
x=[input_thresh Inf];
mink=zeros(1,length(x));

for slice=1:numslices
   slice
   flip=3-flip;
   if slice<numslices
      if ismask
         m=fmris_read_image(m_file,slice+1,1);
         mask(:,flip)=reshape(m.data>m_thresh,IJ,1);
      else
         mask(:,flip)=ones(IJ,1);
      end
      m=fmris_read_image(X_file,slice+1,1);
      X(:,flip)=reshape(m.data,IJ,1);
   else
      mask(:,flip)=zeros(IJ,1);
      X(:,flip)=repmat(-Inf,IJ,1);
   end
   
   pc=IJ*(2-flip);
   pf=(1+pc):(IJ+pc);
   p=find(mask(pf))+pc;
   minX=X(p);
   mink00=histec(minX,x);
   
   ec=(3*IJ-2*I-2*J+1)*(2-flip);
   ef=(1+ec):(6*IJ-3*I-3*J+1+ec);
   edge=find(mask(edge1(ef)) & mask(edge2(ef)))+ec;
   minX=min([X(edge1(edge)); X(edge2(edge))]);
   mink10=histec(minX,x);
   
   tc=2*(IJ-I-J+1)*(2-flip);
   tf=(1+tc):(10*IJ-8*I-8*J+6+tc);
   tri=find(mask(tri1(tf)) & mask(tri2(tf)) & mask(tri3(tf)))+tc;
   minX=min([X(tri1(tri)); X(tri2(tri)); X(tri3(tri))]);
   mink20=histec(minX,x);
   
   tet=find(mask(tet1) & mask(tet2) & mask(tet3) & mask(tet4));
   minX=min([X(tet1(tet)); X(tet2(tet)); X(tet3(tet)); X(tet4(tet))]);
   mink30=histec(minX,x);
   
   mink=mink+mink00-mink10+mink20-mink30;
end   

minkk=fliplr(cumsum(fliplr(mink)));
ec=minkk(1:length(input_thresh));

if ismask
   delete(X_file);
   delete(m_file);
end

return

function c=histec(y,x);
if isempty(y)
   c=zeros(1,length(x));
else
   c=histc(y,x);
end
return


