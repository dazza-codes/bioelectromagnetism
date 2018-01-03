function out=bvpdemo6(t,x,xprime,xN,xNprime,flag)
%  Used by BVPDEMO
%
%  gravitation

switch flag
case 'f'
   r3=(x(1)^2+x(2)^2)^1.5;   % r^3
  out=[-10*x(1)/r3;-10*x(2)/r3];
case 'df/dx' 
  r=sqrt(x(1)^2+x(2)^2);
  r3=r^3;
  r5=r^5;
   
  out=10*[-1/r3 + 3*x(1)^2/r5 , 3*x(1)*x(2)/r5 ; 3*x(1)*x(2)/r5 -1/r3 + 3*x(2)^2/r5];
case 'df/dxprime' 
  out=[0 0;0 0];
case 'R' 
  x1=x;          
  x1prime=xprime; 
  %  out=[x1(1)-2 ; x1(2)-1 ; xN(1)-2 ; xNprime(2)/xNprime(1)+x1prime(2)/x1prime(1)];
  out=[x1(1)-2 ; x1(2)-1 ; xN(2)+2 ; xNprime(1)];

case 'dR/dx1'   
  out=[1 0;0 1;0 0;0 0];
case 'dR/dx1prime'   
   x1prime=xprime;
   % out=[0 0;0 0;0 0;-x1prime(2)/x1prime(1)^2 , 1/x1prime(1)];
   out=zeros(4,2);
case 'dR/dxN'   
  out=[0 0;0 0;0 1;0 0];
case 'dR/dxNprime'   
   %out=[0 0;0 0;0 0;-xNprime(2)/xNprime(1)^2 , 1/xNprime(1)];
   out=[0 0;0 0;0 0;1 0];
otherwise
  error('unknown flag');
end
