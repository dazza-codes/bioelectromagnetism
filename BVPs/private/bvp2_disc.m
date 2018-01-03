function res=bvp_disc(x,f,t,dim,N,h,dummy,parvec)  % x must be first

%   This function is used by BVP2 in order to solve 
%   boundary value problems.
%
%   See BVP2

%   Copyright (c) 1999 Guenter Kneisl
%                      University of Technology, Vienna
%                      e9425595@fbma.tuwien.ac.at

if isempty(parvec)
   parstr=')';
else
   parstr=',parvec)';
end


res=zeros(dim*N,1);  % discretisation residuum

x0  =x(      1:  dim);
x1  =x(  dim+1:2*dim);
x2  =x(2*dim+1:3*dim);

xNm1=x(end-3*dim+1:end-2*dim);
xN  =x(end-2*dim+1:end-  dim);
xNp1=x(end-  dim+1:end      );

x1prime=(x2-x0)/(2*h(1));
xNprime=(xNp1-xNm1)/(2*h(end));

res(1:2*dim)=eval([f '(0,x1,x1prime,xN,xNprime,''R''' parstr] );  % residuum of the boundary condition

lastx=x0;
thisx=x1;

for i=2:N-1
   i_dim=i*dim;
   hav=h(i-1)+h(i);  % average h
   
   nextx = x(i_dim+1:i_dim+dim);
   
   thisf = eval([f '(t(i-1)+t(i+1)/2,(lastx+nextx)/2,(nextx-lastx)/hav,0,0,''f''' parstr]);
   
   res(i_dim+1:i_dim+dim)=-h(i)/hav*lastx + thisx -h(i-1)/hav*nextx + h(i-1)*h(i)/2*thisf;
   
   lastx = thisx;
   thisx = nextx;
end
