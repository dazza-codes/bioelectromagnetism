function J=disc_jacobian(x,f,t,dim,N,h2,bvpopt,parvec)  % x must be first

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


J=sparse([],[],[],dim*N,dim*N,0);  % Jacobian of the discretisation

x1=x(1:dim);
xN=x(end-dim+1:end);

J(1:dim,1:dim)=eval([f '(0,x1,xN,''dR/dx1''' parstr]);             % Randwerte
J(1:dim,end-dim+1:end)=eval([f '(0,x1,xN,''dR/dxN''' parstr]);

J(dim+1:2*dim,1:dim)= - h2(1) * eval([f '(t(1),x1,0,''df/dx''' parstr]);  % die ersten dim Spalten


for i=2:N-1
   i_dim=i*dim;
   
   thisjac = eval([f '(t(i),x(i_dim-dim+1:i_dim),0,''df/dx''' parstr]);  % x(...)=thispt
   
   J(i_dim-dim+1:i_dim+dim,i_dim-dim+1:i_dim) = [-h2(i-1)*thisjac ; -h2(i)*thisjac];
end

J(end-dim+1:end,end-dim+1:end) = - h2(end) * eval([f '(t(N),xN,0,''df/dx''' parstr]);  % die letzten dim Spalten

E=zeros(dim*N,dim+1);    % Jetzt noch die Einheitsmatrizen der Diskretisierung: 
E(:,1)=-1;               % y_{j+1} - y_j - h_j/2*(f(t_j,y_j)+f(t_{j+1},y_{j+1})) = 0
E(:,end)=1;
E(1:dim,end)=0;


J=(J+spdiags(E,-dim:0,dim*N,dim*N))';  % Transponieren f. fsolve

if bvpopt.fsolve(9)~=0  % check jacobian
   J=full(J);  
end
