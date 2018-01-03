function sf_spy(S, varargin)
% SF_SPY	Visualize sparsity structure (color points by magnitude)
%		SF_SPY(S [, label], [, mode] [, threshold])
%                 LABEL = string for figure title (appended)
%                 MODE = {'abs', 'mag', 'pha', 'real', 'imag'}
%		  THRESHOLD - ignore smaller values in S 

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES:
% - to see distribution of values:  plot(full(sort(abs(S(:)))))

%%% THINGS TO DO:
% ? min element color should be different from 0 in sparse matrices
% ? deal with multiple array arguments
%   change (S, la) to (S1, la1) and maintain count
% - tweak crossover point between full/sparse versions
%   for 0.1, full is faster than sparse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments 
if (nargin < 1) 
    help sf_spy; return; 
end;

%%% interpret optional arguments
la = '';
tr = 0;

for j = 1:nargin-1
    aj = varargin{j};
    if ischar(aj)
    	switch aj
    	    case 'abs',	    S = abs(S);
	    case 'mag',	    S = abs(S);
	    case 'pha',	    S = angle(S);
	    case 'real',    S = real(S);
	    case 'imag',    S = imag(S);
	    otherwise,	    la = aj;
	end
    elseif prod(size(aj)) == 1	tr = aj;
    else
	error(['Unable to process argument.']);
    end
    clear aj;
end

if ~isreal(S) error(['Input matrix is complex.']); end

%%% compute figure size
xsize = 400;
ysize = 400;
if prod(size(S)) ~= 1
    if (size(S, 1) > size(S, 2))	ysize = ysize * size(S, 1) / size(S, 2);
    else				xsize = xsize * size(S, 2) / size(S, 1);
    end
end

figure(	'Name',		    	'Matrix structure', ...
	'Position',	    	[0, 0, xsize, ysize], ...
	'PaperPositionMode',	'auto');
S = S .* (abs(S) > tr);

disp(['  min= ', num2str(min(min(S))), ...
       ' thr= ', num2str(tr), ...
       ' max= ', num2str(max(max(S)))]);
disp(['  ', num2str(100*nnz(S)/prod(size(S))), '% full - ', ...
	    num2str(nnz(S), 6), ' of ', num2str(prod(size(S)), 6)]);

if (min(size(S)) > 1)
    %%% display matrix structure
    [I, J, V] = find(S);
    if length(V) < 0.1 * prod(size(S))
	%%% plot each patch
	I = I'; J = J'; V = V';	
	I = [ I; I; I+1; I+1 ];	
	J = [ J; J+1; J+1; J ];
	disp('patch');
	patch(J, I, V);
    else
	%%% pad matrix with zeros and plot
	S = [ S, zeros(size(S, 1), 1); zeros(1, size(S, 2)), 0 ];
	disp('pcolor');
	pcolor(S);
    end
    axis('equal', 'ij');
    colorbar;
    title(la);

else
    %%% display vector structure (plus sorted distribution
    disp('plot')
    plot(1:length(S),      S ,'r+', ...
	 1:length(S), sort(S),'g-');
end


%%% adjust options and add tools
colormap('gray');
shading flat
sf_tool cmap misc zoom;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

