function [F, I] = sf_gen_map(intri, intro, mode, p1, p2, p3, p4)
% SF_GEN_MAP	Generate mappings between representations
%		[F, I] = SF_GEN_MAP(intri, intro [, mode [, p1 ...]])
%		computes forward and inverse array of mapping coefficients 
%		between INTRinsic Input and Output coordinates
%		  MODE = {'ft', 'ftx', 'fty', 'ftxy'}
%		  pn = mode-dependent parameters:

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - the two matrices are conjugate transposes of each other
% - sum(abs(F)) = sum(abs(F')) = area of domain

%%% THINGS TO DO
% - add wrapping on either/both axes to prevent Gibbs
% - use sparse matrices for localized transforms when possible
% - add arguments to accept and incorporate area estimates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 2) help sf_gen_map; return; end
if (nargin < 3) mode = 'ft'; end

%%% check for compatible array sizes
if     size(intri, 2) ~= 2
    error(['intri must have 2 columns, not ', num2str(size(intri,2)), '.']);
elseif size(intro, 2) ~= 2
    error(['intro must have 2 columns, not ', num2str(size(intro,2)), '.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute coefficients for Fourier transform
%%%   I = F' ( ' is conj transpose in MATLAB)

if strcmp(mode, 'ft')
    %%% standard Fourier transform
    F = exp(-i * 2 * pi * intro * intri');
%   I = exp( i * 2 * pi * intri * intro');
    I = F';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    error(['Unknown mode (', mode, ').']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

    
