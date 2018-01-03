function out=bvpdemo1(t,x,xN,flag,parvec)
%  used by BVPDEMO
%
%  x'=x , x(0)+x(N)=2

switch flag
case 'f'  
   out=x;
case 'df/dx' 
   out=1;
case 'R'
   x1=x;
   out=x1+xN-1;
case 'dR/dx1'   
   out=1;
case 'dR/dxN'
   x1=x;
   out=1;
otherwise
   error('unknown flag');
end
