function outmat=normalize(inmat,bval,tval)
% function outval=normalize(inmat)
% Normalises the matrix in inmat to between bval (for example 0)
% and tval (for example 1)
if (nargin==1)
   bval=0;
   tval=1;
end

% get dimensions of input matrix
[x y z]=size(inmat);

minval=min(min(min(inmat))); % only want matrices up to 3D

% do normalisation
inmat=inmat-minval;
maxval=max(max(max(inmat)));
if (maxval==0) 
   maxval=0.000001;
end

inmat=inmat*(tval-bval)/maxval;
outmat=inmat+bval;
