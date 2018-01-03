function res=bvp_disc(x,f,t,dim,N,h2,dummy,parvec)  % x must be first

%   This function is used by BVP in order to solve 
%   boundary value problems.
%
%   See BVP

%   Copyright (c) 1999 Guenter Kneisl
%                      University of Technology, Vienna
%                      e9425595@fbma.tuwien.ac.at


if isempty(parvec)
   parstr=')';
else
   parstr=',parvec)';
end


res=zeros(dim*N,1);  % discretisation residuum

x1=x(1:dim);
xN=x(end-dim+1:end);

res(1:dim)=eval([f '(0,x1,xN,''R''' parstr] );  % residuum of the boundary condition

thispt = x1;
thisf= eval([f '(t(1),thispt,0,''f''' parstr]);

for i=2:N
   i_dim=dim*i;
   
   nextpt = x(i_dim-dim+1:i_dim);
   nextf  = eval([f '(t(1),nextpt,0,''f''' parstr]);
   
   res(i_dim-dim+1:i_dim)=nextpt - thispt - h2(i-1)*(nextf + thisf);
   
   thispt=nextpt;
   thisf=nextf;
end
