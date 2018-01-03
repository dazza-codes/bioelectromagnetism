function [nc,ns,nb,ne,nf] = sf_sub_surf(oc,os,ob,oe,of, varargin)
% SF_SUB_SURF	Subdivide surface (add points, edge, faces)
%   	    	[nc,ns,nb,ne,nf] = sf_sub_surf(oc,os,ob,oe,of [,mode[,mval]])
%   	    	add points to old surface (OC,OS,OB,OE,OF) 
%   	    	  to produce new surface (NC,NS,NB,NE,NF)
%   	MODE 
%    	  'edge' or 'face' - divide longest edges or largest faces (def=edge)
%   	  'extr' or 'intr' - divide based on extr/intr coordinates (def=extr)
%   	  'plot'    	   - plot progress
%   	MODE (with corresponding MVAL)
%   	  'eadd','fadd','padd' - stop when # of new edge/face/pt >= MVAL
%    	  'efac','ffac','pfac' - stop when total # of edge/face/pt 
%   	    	    	    	 increase by factor >= MVAL (1.2 = 20% increase)
%    	  'etot','ftot','ptot' - stop when total # of edge/face/pt >= MVAL
%   	  'mlen','marea'       - stop when max elen/farea < MVAL
%   	  'pass'    	       - stop when # passes >= MVAL

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - face methods aren't good 
%   - faces with small area but long edges can be split poorly

%%% THINGS TO DO
% - assign border based on # of adjacent edges/faces, or status of adjacent points
% ? use BORDER for something useful
% ? improve memory allocation (realloc's use 20% of time)
%   ? replace with C function (dummy values for removed items, enlarge in blocks)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if nargin < 5
    help sf_sub_surf; return;
elseif size(oc,2) ~= 2
    error(['Invalid number of columns (', num2str(size(oc, 2)), ') in OC.']);
elseif size(os,2) ~= 3
    error(['Invalid number of columns (', num2str(size(os, 2)), ') in OS.']);
elseif size(oc,1) ~= size(os,1)
    error(['OC and OS do not have equal number of rows (', ...
    	   num2str(size(oc,1)), ' ~= ', num2str(size(os,1)), ').']);
elseif size(ob,2) ~= 1
    error(['Invalid number of columns (', num2str(size(ob, 2)), ') in OB.']);
elseif size(oc,1) ~= size(ob,1)
    error(['OC and OB do not have equal number of rows (', ...
    	   num2str(size(oc,1)), ' ~= ', num2str(size(ob,1)), ').']);
elseif size(oe,2) ~= 2
    error(['Invalid number of columns (', num2str(size(oe, 2)), ') in OE.']);
elseif size(of,2) ~= 3
    error(['Invalid number of columns (', num2str(size(of, 2)), ') in OF.']);
end

%%% process variable arguments
edgeflg = 1;
extrflg = 0;
plotflg = 0;
mode = 'pfac';
mval = 1.1;

for j=1:nargin-5
    aj = varargin{j};
    switch aj
    	case 'edge', edgeflg=1;
    	case 'face', edgeflg=0;
    	case 'extr', extrflg=1;
    	case 'intr', extrflg=0;
    	case 'plot', plotflg=1;
    	case {'eadd','fadd','padd','efac','ffac','pfac', ...
    	      'etot','ftot','ptot','mlen','marea','pass'},
    	    mode = aj;
    	    mval = varargin{j+1};
    	    if ~isnumeric(mval) | prod(size(mval)) ~= 1
    	    	error(['Invalid argument for mode ', aj, '.']);
    	    end
   end % switch aj
end % for j

  	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

nc = oc; ns=os; nb=ob; ne=oe; nf=of;
pass = 0;
fh1 = [];
if plotflg, sf_plot(nc,zeros(size(nc,1),1)); end

%%% precompute length/area as necessary
if edgeflg
    if extrflg	elen = sf_sub_surf_elen(ns,ne);
    else    	elen = sf_sub_surf_elen(nc,ne); end
else
    if extrflg  flen = sf_sub_surf_flen(ns,nf);
    else    	flen = sf_sub_surf_flen(nc,nf); end
    farea = sf_sub_surf_farea(flen);
end % if edgeflg

%%% convert modes to simplify tests
switch mode
    case {'etot','ftot','ptot','pass'},
    case 'eadd',    mode = 'etot'; mval = mval + size(ne,1);
    case 'fadd',    mode = 'ftot'; mval = mval + size(nf,1);
    case 'padd',    mode = 'ptot'; mval = mval + size(nc,1);
    case 'efac',    mode = 'etot'; mval = mval * size(ne,1);
    case 'ffac',    mode = 'ftot'; mval = mval * size(nf,1);
    case 'pfac',    mode = 'ptot'; mval = mval * size(nc,1);
    case 'marea',
    	if  edgeflg, error(['Mode marea not allowed for edge-based division.']); end
    case 'mlen',    
        if ~edgeflg, error(['Mode mlen not allowed for face-based division.']); end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% subdivision loop starts here

while 1
    
    %%% find longest edge or largest face to subdivide
    if edgeflg
    	[ev,ei] = max(elen);
    	e0 = ne(ei,:);
    else
    	[fv,fi] = max(farea);
    	% find longest edge of largest face
    	[ev,ei] = max(flen(fi,:));
    	e0 = nf(fi,[ei mod(ei,3)+1]);
    	ei = find(any(ne == e0(1),2) & any(ne == e0(2),2));
    end
    
    %%% check whether subdivision is finished
    switch mode
    	case  'etot', if size(ne,1) >= mval, break; end
    	case  'ftot', if size(nf,1) >= mval, break; end
    	case  'ptot', if size(nc,1) >= mval, break; end
    	case  'mlen', if         ev <  mval, break; end
    	case 'marea', if         fv <  mval, break; end
    	case  'pass', if       pass >= mval, break; end
    end % switch MODE
    	    
    %%% remove 1 edge
    ne(ei,:) = [];
    
    %%% find and remove 2 faces along edge (1 at border)
    %   (40% of time was spent in next line)
%    fi = find(any(nf == e0(1),2) & any(nf == e0(2),2));
    fi = sf_find_vals(e0,nf);
    f0 = nf(fi,:);
    nf(fi,:) = [];

    %%% add 1 point at midpoint of edge (intr and extr coordinates)
    p1 = 1+size(nc,1);
    nc(p1,:) = sum(nc(e0,:))/2;
    ns(p1,:) = sum(ns(e0,:))/2;
    nb(p1,:) = length(fi) < 2;
    
    %%% add 4 edges (3 at border) (net gain of 3 or 2)
    %%% (2 to old edge, 1-2 to far corners of faces)
    e1 = [e0'; f0(f0 ~= e0(1) & f0 ~= e0(2))];
    e1 = [e1 p1(ones(size(e1,1),1))];
    ne = [ne; e1];

    %%% add 4 faces (2 at border) (net gain of 2 or 1)
    f1 = f0; f1(f1==e0(1)) = p1;
    f2 = f0; f2(f2==e0(2)) = p1;
    nf = [nf; f1;f2];

    %%% update length/area (remove old, add new)
    if edgeflg
    	elen(ei) = [];
    	if extrflg  elen = [elen; sf_sub_surf_elen(ns,e1)];
    	else	    elen = [elen; sf_sub_surf_elen(nc,e1)]; end
    else
    	if extrflg  flen1 = sf_sub_surf_flen(ns,[f1;f2]);
    	else	    flen1 = sf_sub_surf_flen(nc,[f1;f2]); end
    	 flen(fi,:) = [];  flen = [ flen; flen1];
    	farea(fi,:) = []; farea = [farea; sf_sub_surf_farea(flen1)];
    end % if edgeflg
    
    pass = pass + 1;

    % display current status
    if plotflg, fh1 = sf_plot(nc,nf,fh1); pause(0); end

end % while 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting functions start here

function elen  = sf_sub_surf_elen(p,e)
    elen = sqrt(sum(diff(reshape(p(e             ,:), ...
    	    	    	    	 [size(e)     size(p,2)]),1,2).^2,3));

function flen  = sf_sub_surf_flen(p,f)
    flen = sqrt(sum(diff(reshape(p(f(:,[1 2 3 1]),:), ...
    	    	    	    	 [size(f,1) 4 size(p,2)]),1,2).^2,3));

function farea = sf_sub_surf_farea(flen)
    fs = sum(flen,2)/2;
    farea = sqrt(prod([fs, fs(:,[1 1 1]) - flen],2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

