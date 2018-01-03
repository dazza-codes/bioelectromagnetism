function result = aniso3(input,psiFunction,sigma,iterations,lambda)
%
% result = aniso3(input,[psiFunction],[sigma],[iterations],[lambda])
%
% Anisotropic diffusion in 3d, following the robust statistics
% formulation laid out in Black et al, IEEE Trans on Im Proc,
% 7:421-432, 1998.
%
% Examples of how to run this function are included at the end of
% aniso3.m 
%
%   input: input volume
%   psiFunction: influence function that determines how the
%      diffusion depends on the local image gradient.  Three
%      example psi functions are:
%         linearPsi.m
%         lorentzianPsi.m.
%         tukeyPsi.m
%      Default is tukeyPsi, that makes this the same as one of
%      original two algorithms proposed by Perona et al.  But
%      tukeyPsi is often a better choice because it does a better
%      job of maintaining the sharpness of an edge.  LinearPsi
%      gives standard linear diffusion, i.e., shift-invariant
%      convolution by a Gaussian. 
%   sigma: scale parameter on the psiFunction.  Choose this
%      number to be bigger than the noise but small than the real
%      discontinuties. [Default = 1]
%   iterations: number of iterations [Default = 10]
%   lambda: rate parameter [Default = 1/4] To approximage a
%      continuous-time PDE, make lambda small and increase the
%      number of iterations.
%
% DJH, 8/2000

if ~exist('psiFunction','var')
   psiFunction = 'tukeyPsi';
end
if ~exist('sigma','var')
   sigma = 1;
end
if ~exist('iterations','var')
   iterations = 10;
end
if ~exist('lambda','var')
   lambda = 0.25;
end
lambda = lambda/6;

% Initialize result
result = input;

% Indices for the center pixel and the 4 nearest neighbors
% (north, south, east, west)
[m,n,p] = size(input);
rowC = [1:m]; rowN = [1 1:m-1]; rowS = [2:m m];
colC = [1:n]; colE = [1 1:n-1]; colW = [2:n n];
sliC = [1:p]; sliD = [1 1:p-1]; sliU = [2:p p];

for i = 1:iterations
   % Compute difference between center pixel and each of the 4
   % nearest neighbors.
   north = result(rowN,colC,sliC)-result(rowC,colC,sliC);
   south = result(rowS,colC,sliC)-result(rowC,colC,sliC);
   east = result(rowC,colE,sliC)-result(rowC,colC,sliC);
   west = result(rowC,colW,sliC)-result(rowC,colC,sliC);
   up = result(rowC,colC,sliU)-result(rowC,colC,sliC);
   down = result(rowC,colC,sliD)-result(rowC,colC,sliC);
   % Evaluate the psiFunction for each of the neighbor
   % differences and add them together.  If the local gradient is
   % small, then the psiFunction should increase roughly linearly
   % with the neighbor difference.  If the local gradient is large
   % then the psiFunction should be zero (or close to zero) so
   % that large gradients are ignored/treated as outliers/stop the
   % diffusion.
   psi = eval([psiFunction,'(north,sigma)']) + ...
      eval([psiFunction,'(south,sigma)']) + ...
      eval([psiFunction,'(east,sigma)']) + ...
      eval([psiFunction,'(west,sigma)']) + ...
      eval([psiFunction,'(up,sigma)']) + ...
      eval([psiFunction,'(down,sigma)']);
   % Update result
   result = result + lambda * psi;
end;
return


function y = linearPsi(x,sigma)
y = 2*x;
return;


function y = lorentzianPsi(x,sigma)
y = (2*x)./(2*sigma^2 + abs(x).^2);
return


function y = tukeyPsi(x,sigma)
y = zeros(size(x));
id = (x > -sigma) & (x < sigma);
xid = x(id);
y(id) = xid.*((1-(xid/sigma).^2).^2);
return


% Test, debug, examples

% Noisy step
step = zeros([32 32 32]);
step(:,:,[1:16]) = ones([32 32 16]); 
noise = 0.1 * randn(size(step));
vol = (1 + step + noise);

% Small number of iterations.  Not much difference.
resultLin10 = aniso3(vol,'linearPsi',0,10);
resultLor10 = aniso3(vol,'lorentzianPsi',0.5,10);
resultTuk10 = aniso3(vol,'tukeyPsi',0.5,10);
figure(1); clf; 
plot([squeeze(vol(16,16,:)),...
      squeeze(resultLin10(16,16,:)),...
      squeeze(resultLor10(16,16,:)),...
      squeeze(resultTuk10(16,16,:))]);

% More iterations.  Note that tukeyPsi is much more robust
% because it is a fully redescending influence function (see
% Black et al).  The upshot of this is that it maintains the
% sharpness of the edge even after a large number of iterations
% whereas lorentzianPsi does not.
resultLin100 = aniso3(vol,'linearPsi',0,100);
resultLor100 = aniso3(vol,'lorentzianPsi',0.5,100);
resultTuk100 = aniso3(vol,'tukeyPsi',0.5,100);
figure(2); clf; 
plot([squeeze(vol(16,16,:)),...
      squeeze(resultLin100(16,16,:)),...
      squeeze(resultLor100(16,16,:)),...
      squeeze(resultTuk100(16,16,:))]);
