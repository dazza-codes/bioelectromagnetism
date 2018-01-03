function out=bvpdemo5(t,x,xprime,xN,xNprime,flag)
%  Used by BVPDEMO
%
%  cosh(t)

switch flag
case 'f'  
  out=x;
case 'df/dx' 
  out=1;
case 'df/dxprime' 
  out=0;
case 'R' 
  x1=x;             
  out=[x1-cosh(2);xN-cosh(2)];
case 'dR/dx1'   
  out=[1;0];
case 'dR/dx1prime'   
  out=[0;0];
case 'dR/dxN'   
  out=[0;1];
case 'dR/dxNprime'   
  out=[0;0];
otherwise
  error('unknown flag');
end
