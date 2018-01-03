function J=disc_jacobian(x,f,t,dim,N,h,bvpopt,parvec)  % x must be first

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

J=sparse([],[],[],dim*N,dim*N,0);  % Jacobian of the discretisation

x0  =x(      1:  dim);             % some definitions for readability
x1  =x(  dim+1:2*dim);
x2  =x(2*dim+1:3*dim);

xNm1=x(end-3*dim+1:end-2*dim);     % x_{N-1}
xN  =x(end-2*dim+1:end-  dim);
xNp1=x(end-  dim+1:end      );

x1prime=(x2-x0)/(2*h(1));
xNprime=(xNp1-xNm1)/(2*h(end));


dRdx1prime=eval([f '(0,x1,x1prime,xN,xNprime,''dR/dx1prime''' parstr]);

dRdx0=dRdx1prime /(-2*h(1));
dRdx1=eval([f '(0,x1,x1prime,xN,xNprime,''dR/dx1''' parstr]);
dRdx2=-dRdx0;               % =dRdx1prime /(2*h(1));


dRdxNprime=eval([f '(0,x1,x1prime,xN,xNprime,''dR/dxNprime''' parstr]);

dRdxNm1= dRdxNprime /(-2*h(end));
dRdxN  = eval([f '(0,x1,x1prime,xN,xNprime,''dR/dxN''' parstr]);
dRdxNp1= -dRdxNm1;          % =dRdxNprime /(2*h(N-1));




J(1:2*dim,1:3*dim)=[dRdx0 dRdx1 dRdx2];             % Randwerte
J(1:2*dim,end-3*dim+1:end)=[dRdxNm1 dRdxN dRdxNp1];


lastx=x0;
thisx=x1;
eyedim=speye(dim);

for i=2:N-1
   i_dim=i*dim;
   hav=h(i-1)+h(i);  % average h
   
   nextx = x(i_dim+1:i_dim+dim);
   
   dfdx      = eval([f '(t(i-1)+t(i+1)/2,(lastx+nextx)/2,(nextx-lastx)/hav,0,0,''df/dx''' parstr]);
   dfdxprime = eval([f '(t(i-1)+t(i+1)/2,(lastx+nextx)/2,(nextx-lastx)/hav,0,0,''df/dxprime''' parstr]);
   
   J(i_dim+1:i_dim+dim,i_dim-2*dim+1:i_dim+dim) = ...
   [  -h(i)/hav*eyedim + h(i-1)*h(i)/2 * (dfdx/2 - dfdxprime/hav) , eyedim , ...
    -h(i-1)/hav*eyedim + h(i-1)*h(i)/2 * (dfdx/2 + dfdxprime/hav) ];
   
   lastx = thisx;
   thisx = nextx;
 
   
end

J=J';
if bvpopt.fsolve(9)~=0  %check jacobian
   J=full(J);  
end
