function vol = anisoSapiro3(vol,numiters,deltaT,epsilon)
% vol = anisoSapiro3(vol,numiters,[deltaT],[epsilon])
%
% Anisotropic diffusion of a 3D image, using Guillermo Sapiro's method.
%
% vol - input volume
% numiters - number of iterations
% deltaT - time step (default 0.02)
% epsilon - avoids divide by zero in sapiroGrad (default 1e-6)
%
% WAP '99

if ~exist('deltaT')
  deltaT = 0.02;
end

if ~exist('epsilon')
  epsilon = 1e-6;
end

h = waitbar(0,'Performing Sapiro anisotropic smoothing');
for iter = 1:numiters
   result = sapiroGrad3(vol,epsilon);
   if iter==1
      % First time around, calculate convolved volume size
		xoff = (size(vol,1)-size(result,1))/2;
   	yoff = (size(vol,2)-size(result,2))/2;
      zoff = (size(vol,3)-size(result,3))/2;
      % Offsets should be integral, since conv throws out equal
      % number of pixels on each side of each dimension.
      xrange = (1+xoff):(size(vol,1)-xoff);
      yrange = (1+yoff):(size(vol,2)-yoff);
      zrange = (1+zoff):(size(vol,3)-zoff);
   end   
   vol(xrange,yrange,zrange) = vol(xrange,yrange,zrange) + deltaT .* result;
   waitbar(iter/numiters,h);
end
close(h);
return


function result = sapiroGrad3(v,epsilon)
% result = sapiroGrad3(v,[epsilon])
%
% Gradient function used by anisoSapiro3.m 
% fourthroot(curvature)/(gradient+epsilon)
%
% v - input volume
% epsilon - avoids divide by zero
%
% WAP '99

if ~exist('epsilon')
  epsilon = 1e-6;
end

fx = dx3(v);	fxx = dxx3(v);		fxy = dxy3(v);		fx2 = fx.^2;
fy = dy3(v);	fyy = dyy3(v);		fxz = dxz3(v);		fy2 = fy.^2;
fz = dz3(v);	fzz = dzz3(v);		fyz = dyz3(v);		fz2 = fz.^2;

grad = sqrt(fx2 + fy2 + fz2);
curvature = ...
   (fx2.*(fyy.*fzz-fyz.^2) + ...
    fy2.*(fxx.*fzz-fxz.^2) + ...
    fz2.*(fxx.*fyy-fxy.^2) + ...
    2.*fy.*fz.*(fxy.*fxz-fyz.*fxx) + ...
    2.*fx.*fz.*(fxy.*fyz-fxz.*fyy) + ...
    2.*fx.*fy.*(fxz.*fyz-fxy.*fzz)) ./ ...
   (fx2+fy2+fz2).^2;
result = fourthroot(curvature)./(grad+epsilon);
return


function b = fourthroot(a)
b = sign(a) .* (abs(a)).^(1/4);
return
