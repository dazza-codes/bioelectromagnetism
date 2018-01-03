function [intr, extr, edge, face] = sf_gen_surf(num, org, iran, ijit, eran)
% SF_GEN_SURF	Generate surface (both intrinsic and extrinsic)
%		[intr, extr, edge, face] = sf_gen_surf(num, org, iran, ijit, eran, ejit)
%		generates (roughly) NUM*NUM points in (-IRAN/2 ... +IRAN/2)^2 
%		ORG = 
%   	    	    'icos'	- tesselated icosohedron
%   	    	    'isph'  	- sphere (radial proj from icos)
%		IJIT = variance of gaussian jitter for each point
%			 ( multiplied by iran/num )
%		defaults:  IRAN=1.0  IJIT=0.01  ERAN=1.0
%		edges and faces generated if specified

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - face vertices are clockwise from exterior, not in order
% - intrinsic coordinates are completed at end from spherical coordinates
%   (replaces coordinates computed for each new vertex)
% - code for edge-based faces in icos is commented with %% (untested)
% ICOS,ISPH: pi/10 rotation aligns with SF_GEN_EXTR(sphere)
% - to look at intermediate results
%	sf_plot(extr(1:vn, :), edge(1:en, :));
%	sf_plot(extr(1:vn, :), face(1:fn, :));

%%% THINGS TO DO
% - jitter should affect intrinsic coords, and be reflected in extrinsic
%   - currently extr do not correspond to intr for nonzero jitter
%   - THIS SHOULD BE FIXED so sf_gen_surf works like sf_gen_intr,sf_gen_extr
% ? option to use old intrinsic coordinates
% - add extrinsic jitter, which will be converted to intrinsic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 2) help sf_gen_surf; return; end
if (nargin < 5) eran	= 1.00; end
if (nargin < 4) ijit 	= 0.01; end
if (nargin < 3) iran	= 1.00; end
if ~ischar(org) error(['ORG argument must be a string']); end

disp(['WARNING: jitter implemented incorrectly.']);

edge = []; 
face = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate appropriate surface

switch org
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'icos','isph'},
	%%% generate corner vertices
	cv = zeros(12, 3);
	cv(   1, :) = [0, 0, 1];
	cv(2: 6, :) = [  0.8945 * sin(linspace(-4*pi/5, 4*pi/5, 5)' + pi/10), ...
			 0.8945 * cos(linspace(-4*pi/5, 4*pi/5, 5)' + pi/10), ...
			 0.4472 * linspace(1, 1, 5)'];
	cv(7:11, :) = [  0.8945 * sin(linspace(-4*pi/5, 4*pi/5, 5)' + 3*pi/10), ...
			 0.8945 * cos(linspace(-4*pi/5, 4*pi/5, 5)' + 3*pi/10), ...
			-0.4472 * linspace(1, 1, 5)'];
	cv(  12, :) = [0, 0, -1];
	cv = cv / 2;

	%%% allocate storage for vertices, edges, faces
	r = round(num/sqrt(10));
	num = sqrt(r*r*10+2);
	intr = zeros(r*r*10+2, 2);
	extr = zeros(r*r*10+2, 3);
	edge = zeros(r*r*30, 2);
	face = zeros(r*r*20, 3);
	vn = 0; en = 0; fn = 0;
	vf = 1; ef = 1; 

	%%% add N-pole vertex
	intr(vn+1, :) = [0, 0.5];
	extr(vn+1, :) = cv(1, :);
	vn = vn+1;

	%%% for each ring below N-pole and through N-tropic
	for i = 1:r
	    vpf = vf; vpn = vf; vf = vn+1; 
    %%	epn = ef+1; ef = en+1;
	    v2 = cv(1, :) + (cv(6, :) - cv(1, :))*i/r;
	    intr(vn+1, :) = [-0.5, 0.5 - i/(3*r)];
	    extr(vn+1, :) = v2;
	    edge(en+1, :) = [vpn,vn+1];
	    vn = vn+1; en = en+1;

	    for j = 2:6
		v1 = v2; v2 = cv(1, :) + (cv(j, :) - cv(1, :))*i/r;
		for k = 1:i-1
		    intr(vn+1, :) = [-0.5 + (j-1)/5 + ((k/i) -1)/5, 0.5 - i/(3*r)];
		    extr(vn+1, :) = v1 + (v2 - v1)*k/i;
		    edge(en+1, :) = [vn,vn+1];
		    edge(en+2, :) = [vpn,vn+1];
		    % switch to clockwise vertices
		    face(fn+1, :) = [vpn, vn+1, vn];
    %%		face(fn+1, :) = [en,en+1,en+2];
		    vpn = vpn+1;
		    if (vpn < vf)
			edge(en+3, :) = [vpn, vn+1];
			face(fn+2, :) = [vpn-1, vpn, vn+1];
		    else
			edge(en+3, :) = [vpf, vn];
			% switch to clockwise vertices
			face(fn+2, :) = [vpn-1, vpf, vn+1];
		    end
    %%		face(fn+2, :) = [epn,en+2,en+3];
    %%		epn = epn+3; 
		    vn = vn+1; en = en+3; fn = fn+2;
		end

		if (j ~= 6)
		    intr(vn+1, :) = [-0.5 + (j-1)/5, 0.5 - i/(3*r)];
		    extr(vn+1, :) = v2;
		    edge(en+1, :) = [vn, vn+1];
		    edge(en+2, :) = [vpn, vn+1];
		    % switch to clockwise vertices
		    face(fn+1, :) = [vpn, vn+1, vn];
    %%		epn = epn-1; 
		    vn = vn+1; en = en+2; fn = fn+1;
		else
		    edge(en+1, :) = [vf, vn];
		    face(fn+1, :) = [vpf, vf, vn];
		    en = en+1; fn = fn+1;
		end
    %%	    face(fn+1, :) = [en, en+1, ef];
	    end
	end

	%%% for each ring below N-tropic and through S-tropic
	for i = 1:r
	    vpf = vf; vpn = vf; vf = vn+1; 
    %%	epn = ef+1; ef = en+1;
	    v2 = cv(6, :) + (cv(11, :) - cv(6, :))*i/r;
	    intr(vn+1, :) = [-0.5, (1/6) - i/(3*r)];
	    extr(vn+1, :) = v2;
	    edge(en+1, :) = [vpn, vn+1];
	    vn = vn+1; en = en+1;

	    for j = 2:6
		vpn = vpn+1;
		edge(en+1, :) = [vpn, vn];
		face(fn+1, :) = [vpn-1, vpn, vn];
    %%	    face(fn+1, :) = [epn, en, en+1];
    %%	    epn = epn+3; 
		en = en+1; fn = fn+1;
		v1 = v2; v2 = cv(j, :) + (cv(rem(j+2, 5)+7, :) - cv(j, :))*i/r;

		for k = 1:r-i
		    intr(vn+1, :) = [-0.5 + (j-1)/5 + ((k/r) -1)/5, (1/6) - i/(3*r)];
		    extr(vn+1, :) = v1 + (v2 - v1)*k/(r-i);
		    edge(en+1, :) = [vn, vn+1];
		    edge(en+2, :) = [vpn, vn+1];
		    % switch to clockwise vertices
		    face(fn+1, :) = [vpn, vn+1, vn];
    %%		face(fn+1, :) = [en,en+1,en+2];
		    vpn = vpn+1;
		    if (vpn < vf)
			edge(en+3, :) = [vpn, vn+1];
			face(fn+2, :) = [vpn-1, vpn, vn+1];
		    else
			edge(en+3, :) = [vpf, vn+1];
			% switch to clockwise vertices
			face(fn+2, :) = [vpn-1, vpf, vn+1];
		    end
    %%		face(fn+2, :) = [epn,en+2,en+3];
    %%		epn = epn+3; 
		    vn = vn+1; en = en+3; fn = fn+2;
		end

		v1 = v2; v2 = cv(j, :) + (cv(j+5, :) - cv(j, :))*i/r;

		for k = r-i+1:r-1
		    intr(vn+1, :) = [-0.5 + (j-1)/5 + ((k/r) -1)/5, (1/6) - i/(3*r)];
		    extr(vn+1, :) = v1 + (v2 - v1)*(k+i-r)/i;
		    edge(en+1, :) = [vn, vn+1];
		    edge(en+2, :) = [vpn, vn+1];
		    % switch to clockwise vertices
		    face(fn+1, :) = [vpn, vn+1, vn];
    %%		face(fn+1, :) = [en,en+1,en+2];
		    vpn = vpn+1;
		    if (vpn < vf)
			edge(en+3, :) = [vpn, vn+1];
			face(fn+2, :) = [vpn-1, vpn, vn+1];
		    else
			edge(en+3, :) = [vpf, vn+1];
			% switch to clockwise vertices
			face(fn+2, :) = [vpn-1, vpf, vn+1];
		    end
    %%		    face(fn+2, :) = [epn,en+2,en+3];
    %%		epn = epn+3; 
		    vn = vn+1; en = en+3; fn = fn+2;
		end

		if (j ~= 6)
		    intr(vn+1, :) = [-0.5 + (j-1)/5, (1/6) - i/(3*r)];
		    extr(vn+1, :) = v2;
		    edge(en+1, :) = [vn, vn+1];
		    edge(en+2, :) = [vpn, vn+1];
		    % switch to clockwise vertices
		    face(fn+1, :) = [vpn, vn+1, vn];
    %%		face(fn+1, :) = [en, en+1, en+2];
    %%		if (i == 1) epn = epn-1; end
		    vn = vn+1; en = en+2; fn=fn+1;
		else
		    edge(en+1, :) = [vf, vn];
		    face(fn+1, :) = [vpf, vf, vn];
    %%		face(fn+1, :) = [en, en+1, ef];
    %%		ef = ef+1; 
		    en=en+1; fn=fn+1;
		end
	    end
	end

	%%% for each ring below S-tropic and above S-pole
	for i = r-1:-1:1
	    vpf = vf; vpn = vf; vf = vn+1; 
    %%	epn = ef+1; ef = en+1;
	    v2 = cv(12, :) + (cv(11, :) - cv(12, :))*i/r;
	    intr(vn+1, :) = [-0.5, -(1/6) - (r-i)/(3*r)];
	    extr(vn+1, :) = v2;
	    edge(en+1, :) = [vpn, vn+1];
	    vn = vn+1; en = en+1;

	    for j = 7:11
		vpn = vpn+1;
		edge(en+1, :) = [vpn, vn];
		face(fn+1, :) = [vpn-1, vpn, vn];
    %%	    face(fn+1, :) = [epn, en, en+1];
    %%	    epn = epn+3; 
		en = en+1; fn = fn+1;
		v1 = v2; v2 = cv(12, :) + (cv(j, :) - cv(12, :))*i/r;

		for k = 1:i-1
		    intr(vn+1, :) = [-0.5 + (j-6)/5 + ((k/i)-1)/5, -(1/6) - (r-i)/(3*r)];
		    extr(vn+1, :) = v1 + (v2 - v1)*k/i;
		    edge(en+1, :) = [vn, vn+1];
		    edge(en+2, :) = [vpn, vn+1];
		    % switch to clockwise vertices
		    face(fn+1, :) = [vpn, vn+1, vn];
    %%		face(fn+1, :) = [en,en+1,en+2];
		    vpn = vpn+1;
		    edge(en+3, :) = [vpn, vn+1];
		    face(fn+2, :) = [vpn-1, vpn, vn+1];
    %%		face(fn+2, :) = [epn,en+2,en+3];
    %%		epn = epn+3; 
		    vn = vn+1; en = en+3; fn = fn+2;
		end

		if (j ~= 11)
		    intr(vn+1, :) = [-0.5 + (j-6)/5, -(1/6) - (r-i)/(3*r)];
		    extr(vn+1, :) = v2;
		    edge(en+1, :) = [vn, vn+1];
		    edge(en+2, :) = [vpn, vn+1];
		    % switch to clockwise vertices
		    face(fn+1, :) = [vpn, vn+1, vn];
    %%		face(fn+1, :) = [en,en+1, en+2];
		    vpn = vpn+1;
		    edge(en+3, :) = [vpn, vn+1];
		    face(fn+2, :) = [vpn-1, vpn, vn+1];
    %%		face(fn+2, :) = [epn,en+2,en+3];
    %%		if (i == r-1) epn = epn+3; else epn = epn+4; end
		    vn = vn+1; en = en+3; fn = fn+2;
		else
		    edge(en+1, :) = [vf, vn];
		    edge(en+2, :) = [vpn, vf];
		    face(fn+1, :) = [vpn, vf, vn];
		    % switch to clockwise vertices
		    face(fn+2, :) = [vpn, vpf, vf];
    %%		face(fn+1, :) = [en, en+1, en+2];
    %%		face(fn+2, :) = [epn, en+2, ef];
    %%		ef = ef+1; 
		    en = en+2; fn = fn+2;
		end
	    end
	end

	vpf = vf; vpn = vf; vf = vn+1; 
    %%  epn = ef+1; ef = en+1;
	intr(vn+1, :) = [0, -0.5];
	extr(vn+1, :) = cv(12, :);
	edge(en+1, :) = [vpn, vn+1];
	vpn = vpn+1; 
	vn = vn+1; en = en+1;
	for j=1:4
	    edge(en+1, :) = [vpn, vn];
	    face(fn+1, :) = [vpn-1, vpn, vn];
    %%	face(fn+1, :) = [epn, en, en+1];
	    vpn = vpn+1; 
    %%	epn = epn+3; 
	    en = en+1; fn = fn+1;
	end

    %%  face(fn+1, :) = [epn, en, ef];
	% switch to clockwise vertices
	face(fn+1, :) = [vpn-1, vpf, vn];
	fn = fn+1;

	% normalize radius to convert icosohedron to sphere
	if strcmp(org,'isph')
    	    rad = sqrt(sum(extr.^2,2));
    	    extr = [extr(:,1) ./ rad, extr(:,2) ./ rad, extr(:,3) ./ rad] / 2;
	end

	% compute and adjust intrinsic coordinates
	intr = [ atan2(extr(:,1),extr(:,2)) ./ (2*pi), ...
    		 asin(extr(:,3) ./ sqrt(sum(extr.^2,2))) ./ (pi) ];

	intr = intr*iran + randn(size(intr))*(ijit*iran/num);
	extr = extr*eran;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise,
    	error(['Unknown organization: ', org, '.']);

end % switch ORG


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print sampling values for both domains

disp(['sampling:  num=', num2str(num), ...
		' ran=', num2str(iran), ...
		' min=', num2str(iran/num), ...
		' max=', num2str(0.5*iran)]);
disp([' dual ->   num=', num2str(num), ...
		' ran=', num2str(num/iran), ...
		' min=', num2str(1/iran), ...
		' max=', num2str(0.5*num/iran)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate edges and faces if they don't exist

if (nargout > 1) & all(size(edge) == [0 0])
    [edge,face] = sf_delaunay(intr, 5*range/num);
end

if (nargout > 1) & all(size(edge) == [0 0]) & 0
    edge = sf_gen_edge(intr, 5*iran/num)
end
if (nargout > 2) & all(size(face) == [0 0]) & 0
    face = sf_gen_face(intr, edge)
end

disp([' intr=', num2str(length(intr)), ' extr=', num2str(length(extr)), ...
      ' edge=', num2str(length(edge)), ' face=', num2str(length(face))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

