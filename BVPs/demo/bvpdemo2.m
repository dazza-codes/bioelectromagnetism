function out=bvpdemo2(t,x,xN,flag,parvec)

%  used by BVPDEMO
%  pendulum

switch flag
case 'f'  
   out=[x(2) ; -sin(x(1))];
case 'df/dx' 
   out=[0 1;-cos(x(1)) 0];
case 'R'
   x1=x;
   out=[x1(1)-xN(1);x1(2)-xN(2)];
case 'dR/dx1'   
   out=[1 0;0 1];
case 'dR/dxN'
   out=[-1 0;0 -1];
otherwise
   error('unknown flag');
end
