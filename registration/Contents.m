% Image Registration Tools.
% Last update:   2/9/99
% Authors: David Heeger, heeger@stanford.edu
%          Oscar Nestares, oscar@white.stanford.edu
% 
% See README file for brief description.
% Type "help <command-name>" to get help on individual commands.
%
% ---------------------------------------------------------------
% circularShift3 - shift a 3D array by integer number of samples.
% computeDerivatives2 - internal function that computes x-, y-,
%			   and t-derivatives of a pair of images.
% computeDerivatives3 - internal function that computes x-, y-,
%			   z-, and t-derivatives of a pair of
%			   volumes.
% convXYsep - 2D convolution on each sub-image of a volume.
% convZ - 1D convolution in the Z direction on a volume.
% 
% estMotion2 - 2D motion estimation for a pair of images.  Affine
%    or rigid (translation and rotation) motion.  Least-squares or
%    robust.
% estMotion3 - 3D motion estimation for a pair of volumes.
% estMotionIter2 - iterative 2D motion estimation
% estMotionIter3 - iterative 3D motion estimation
% estMotionMulti2 - multiscale/iterative 2D motion estimation
% estMotionMulti3 - multiscale/iterative 3D motion estimation
% 
% reduce - separable convolution and subsampling by a factor of 2. 
% robustMest - robust M-estimator using Beaton/Tukey's weighting
%    function.
% translateAffine2 - Given an 2D affine transform estimated from
%    a section of an image, compute the alignment matrix that should
%    be applied to the original entire image. 
% warpAffine2 - warp an image according to a 2D affine transform.
% warpAffine3 - warp a volume according to a 3D affine transform.
