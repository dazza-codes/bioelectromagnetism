function out=bvpdemo4(t,x,xN,flag,parvec)
% Used by BVPDEMO
%
% homoclinic orbit


E0=parvec(1);
c=parvec(2);

switch flag
case 'f'  
  n=x(1);
  w=x(2);
  E=x(3);

  out = [w  ; (v(E)-c).*w + vprime(E).*(n-1).*n ; n-1];

case 'df/dx' 
  n=x(1);
  w=x(2);
  E=x(3);

  out=[0 1 0; vprime(E)*(2*n-1) , v(E)-c , vprime(E)*w + vdblprime(E)*(n^2-n) ;1 0 0];
  
case 'R' %boundary condition
   x1=x;
   %out=[x1(1)-1;x1(3)-E0;xN(1)-1]; %  Dirichlet: n(0)=1 E(0)=E0 n(end)=1
   out=x1-xN;                      %  periodische RDBDGen
case 'dR/dx1'   
   %out=[1 0 0;0 0 1;0 0 0];        %  Dirichlet: n(0)=1 E(0)=E0 n(end)=1
   out=eye(3);                     %  periodische RDBDGen
case 'dR/dxN'
   %out=[0 0 0;0 0 0;1 0 0];        %  Dirichlet: n(0)=1 E(0)=E0 n(end)=1
   out=-eye(3);                    %  periodische RDBDGen
otherwise
  error('unknown flag');
end




function out=v(E)
out=-E/(E+1)^2;

function out=vprime(E)
out=(E-1)/(E+1)^3;

function out=vdblprime(E)
out=2*(2-E)/(E+1)^4;


