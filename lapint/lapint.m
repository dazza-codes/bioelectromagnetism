function [int, keepindex, repindex] = lapint(lap, index);

% LAPINT: Computes the zero laplacian interpolation matrix
% 
% Useage:   [int, keepindex, repindex] = lapint(lap, index)
% 
% lap is the laplacian matrix for the full mesh (see lapcal)
% int is the matrix which interpolates the points in 'index'
% to the full mesh.  'index' is a row vector of indices 
% into a subset of the vertices used to calculate
% 'lap'.  This subset is where the electric potential
% is known and corresponds to the given electrode
% vertices (eg, index = dsearchn(scalpvert,elecvert)';).
% 
% Where the index contains repeated indices, only the
% unique indices are kept.  The 'keepindex' array can
% be used to select these.  The 'repindex' array is 
% the repeated indices. If they are electrodes, 
% remove them from further interpolation calculations.
% 
% Interpolations can be done using matrix 'int':
% 
% [int, keepindex, repindex] = lapint(lap,index);
% if isempty(repindex),
%   Vint = int * Vknown;
% else
%   Vint = int * Vknown(keepindex);
% end
% 
% This implements interpolation method B (p. 336) of 
% Oostendorp T, Oosterom A & Huiskamp G (1989),
% Interpolation on a triangulated 3D surface.
% Journal of Computational Physics, 80: 331-343.
%

% Licence:  GNU GPL, no implied or express warranties
% History:  (c) 04/2002 Robert Oostenveld
%           - agreed to release into eeg_toolbox under GNU GPL
%           04/2002, Darren.Weber@flinders.edu.au
%           - introduced check for index replications and
%             adjusted calculations/output accordingly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(lap,1)~=size(lap,2), error('laplacian matrix is not square'); end

if size(index,1)~=1, index = index'; end

% Remove any replicate indices from 'index'
[KnownIndex, i, i] = union(index,index);
if ~isempty(i),
    keepindex = sort(i);
    repindex = setdiff(1:length(index),sort(i));
else
    keepindex = 1:length(index);
    repindex = [];
end
KnownIndex = index(keepindex); % unsort KnownIndex

k = length(KnownIndex);
n = length(lap);

% find 'unknown' indices of lap matrix
UnknownIndex = setdiff(1:n, KnownIndex);

% reshuffle rows & columns of lap matrix
lapi = [KnownIndex, UnknownIndex];
lap = lap(lapi, :); % rows
lap = lap(:, lapi); % columns

% Segregate known/unknown portions of lap
L11 = lap(1:k    ,1:k    );
L12 = lap(1:k    ,(k+1):n);
L21 = lap((k+1):n,1:k    );
L22 = lap((k+1):n,(k+1):n);

int = -pinv([L12; L22]) * [L11; L21];

% append the interpolating piece to the identity matrix 
% these take care of the known potentials
int = [eye(k); int];

% reshuffle the columns of the interpolating matrix
[tmp, order] = sort(KnownIndex);
int = int(:,order);

% reshuffle the rows of the interpolating matrix
[tmp, order] = sort([KnownIndex, UnknownIndex]);
int = int(order, :);

return
