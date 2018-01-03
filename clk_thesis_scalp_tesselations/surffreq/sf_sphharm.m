function sph = sf_sphharm(degree,order,theta,phi)
% SF_SPHHARM	Spherical harmonics of given DEGREE and ORDER
%		sh = SF_SPHHARM(degree [, order [, theta, phi]])
%		with one or two arguments, plots spherical harmonics
%		with four arguments, evaluates at given positions

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - Brechbuhler:  theta = latitude (0..pi) phi = longitude (0..2pi)
% - spherical coordinate conversions:
%     x = sin(theta)cos(phi);
%     y = sin(theta)sin(phi);
%     z = cos(theta)

%%% THINGS TO DO
% - transformation required for intrinsic coordinates
%   theta = 0:pi	(transform more difficult)
%   phi   = 0:2pi	(easy to make -pi:pi or -x:x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments

if nargin > 2
    if abs(order) > degree
	error(['abs(order) (', num2str(abs(order)), ...
	       ') is greater than degree (', num2str(degree), ').']);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for one argument, plot all orders for given degree
%%% for two arguments, plot given order and degree

if (nargin == 1) | (nargin == 2)
    [theta,phi] = meshgrid(linspace(0,pi,20),linspace(0,2*pi,20));
    x = sin(theta) .* cos(phi);
    y = sin(theta) .* sin(phi);
    z = cos(theta);
    if (nargin == 1) 
	offset = degree;
	order = -degree:degree;
    else
	offset = -order;
    end
    figure( 'Name', 	'spherical harmonics', ...
	    'Position',	[0, 200, 400, 200 * length(order)]);
    for j = order
	sph = sf_sphharm(degree, j, theta, phi);
	subplot(length(order), 2, 1+2*(offset+j));
%	surfl(theta, phi, real(sph));
	surf(x, y, z, real(sph));
%	pcolor(real(sph));
	subplot(length(order), 2, 2+2*(offset+j));
%	surfl(theta, phi, imag(sph));
	surf(x, y, z, imag(sph));
%	pcolor(imag(sph));
    end
    sph = [];

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for four arguments, evaluate at specified values

elseif (nargin == 4)
    if size(theta) ~= size(phi)
	error(['theta and phi are not the same size (', ...
	       num2str(size(theta)), ' ~= ', num2str(size(phi)), ').']);
    end
    
    rtheta = theta(:);
    rphi = phi(:);
    
    dmo = 1; for j = 1:degree-abs(order) dmo = dmo * j; end;
    dpo = 1; for j = 1:degree+abs(order) dpo = dpo * j; end;
    leg = legendre(degree, cos(rtheta))';
    sph = sqrt((2 * degree - 1) * dmo / (4 * pi * dpo)) ...
	  * leg(:, abs(order)+1) .* exp(i * abs(order) * rphi);
    
    if (order < 0)
	sph = (-1)^(-order) * conj(sph);
    end
    
    sph = reshape(sph, size(theta));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% otherwise print help message

else
    help sf_sphharm;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

