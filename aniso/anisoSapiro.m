function result = anisoSapiro(im,numiters,edges,deltaT,epsilon)
% anisoSapiro: Anisotropic diffusion of an image, using Guillermo Sapiro's
%     method.
%
% result = anisoSapiro(im,numiters,edges,deltaT,epsilon)
% 
%      result=anisoSapiro(im,numiters,deltaT,epsilon)
%      im - input image
%      numiters - number of iterations
%      edges - any valid edge handler for upConv (default is 'circular').
%      deltaT - time step (default 0.02)
%      epsilon - avoids divide by zero in sapiroGrad (default 1e-6)
%
% Example in anisoTest.m
%
% DJH '96

if ~exist('deltaT')
  deltaT = 0.02;
end

if ~exist('edges')
  edges='circular';
end

if ~exist('epsilon')
  epsilon = 1e-6;
end

result = im(:,:);
for iter = 1:numiters
  result = result + deltaT .* sapiroGrad(result,edges,epsilon);
end

return;

function result = sapiroGrad(f,edges,epsilon)
% SAPIROGRAD: Gradient function used by anisoSapiro.m 
% 1/(gradient+epsilon) (fy^2 fxx + fx^2 fyy - 2 fx fy fxy)^(1/3)
%
% result = sapiroGrad(im,edges,epsilon)
% im - input image
% edges- any vali edge handler for upConv (default is 'circular')
% epsilon - avoids divide by zero
%
% DJH '96

if ~exist('edges')
  edges='circular';
end

if ~exist('epsilon')
  epsilon = 1e-6;
end

fx = dx(f,edges);
fxx = dxx(f,edges);
fy = dy(f,edges);
fyy = dyy(f,edges);
fxy = dxy(f,edges);
fx2 = fx.^2;
fy2 = fy.^2;
grad = sqrt(fx2 + fy2);
g = 1 ./ (grad + epsilon);
tmp =fy2.*fxx + fx2.*fyy - 2.*fx.*fy.*fxy;
result = g.* sign(tmp) .* (abs(tmp)).^(1/3);

return;

% Test/debug

% Noisy step
step = 1+[ones([32 64]) ; zeros([32 64])];
noise = 0.1 * randn(size(step));
im = (step + noise);
figure(1)
displayImage(im);
truesize

sap = anisoSapiro(im,40); 
figure(2)
displayImage(sap);
truesize

figure(3)
plot([im(:,32) sap(:,32)]);

figure(4)
hist(sap(:))

