function [Hx] = eeg_interp_sph_spline_h(COS,r)

% [Hx] = eeg_interp_sph_spline_h(COS)
%
% COS is the cosine matrix from elec_cosines
%
% To solve h(x) in Eq. 5 Perrin et al. (1989), which has pseudo-code
% something like this:
%
% h(x) = (1/r^2) * 1/(4*pi) * ...
%        for n=1:inf, sum = sum + ...
%        ( ( (2*n+1)/(n^(m-1) * (n+1)^(m-1)) ) * Pn(x) );
%
% where x is the cosine between two points on a sphere,
% m is a constant > 1 and Pn(x) is the nth degree Legendre
% polynomial.  Perrin et al. (1989) evaluated m=1:6 and 
% recommend m=4, for which the first 20 terms of Pn(x) 
% are sufficient to obtain a precision of 10^-6 for h(x)
%
% Notes:    This function solves h(x), Eq. 5 in
%           Perrin et al (1989).  Electroenceph. & Clin. 
%             Neurophysiology, 72: 184-187. Corrigenda (1990),
%             Electroenceph. & Clin. Neurophysiology, 76: 565.
%             (see comments in the .m file for details).
%
% see elec_cosines, eeg_interp_sph_spline, LegendreP

% $Revision: 1.1 $ $Date: 2008/05/04 18:25:00 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2001 Darren.Weber_at_radiology.ucsf.edu, with
%                   mathematical advice from
%                   Dr. Murk Bottema (Flinders University of SA)
%           10/2003 Darren.Weber_at_radiology.ucsf.edu, with
%                   mathematical advice and LegendreP function from 
%                   Dr. Tom.Ferree_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegversion = '$Revision: 1.1 $';
fprintf('EEG_INTERP_SPH_SPLINE_H [v %s]\n',eegversion(11:15)); tic

m = 4;

% N=9 works fine for 32 electrodes, but for 128 electrodes it should be larger
if max(size(COS))<=32
  N = 9;
elseif max(size(COS))<=64
  N = 14;
elseif max(size(COS))<=128
  N = 20;
else
  N = 32;
end

% first calculate the series to be multiplied by Pn(x)
n = [1:N]';
Series = (2*n+1) ./ [ [ n.^(m-1) ] .* [ (n+1).^(m-1) ]  ];

%fprintf('%12.10f  ',Series) % note how Series diminishes quickly
%0.3750000000  0.0231481481  0.0040509259  0.0011250000  0.0004074074
%0.0001754670  0.0000854136  0.0000455461  0.0000260631  0.0000157776
%0.0000100001  0.0000065852  0.0000044787  0.0000031314  0.0000022425
%0.0000016399  0.0000012215  0.0000009250  0.0000007107  0.0000005534

% Perrin et al. (1989) recommend tabulation of h(x) for
% x = linspace(-1,1,2000) to be used as a lookup given actual
% values for cos(Ei,Ej).
if isempty(COS),
  error('...cosine matrix is empty.\n');
  %msg = sprintf('...cosine matrix empty, generating COS.\n');
  %warning(msg);
  %COS = linspace(-1,1,2000)';
end

Hx = zeros(size(COS));

% *********************************************************
% What about this adjustment!
%
% Perrin et al. (1989) assume r = 1, so we adjust Eq. 5:
% del_perp^2[Pn(x)] = [ -n*(n+1)/r^2 ] Pn(x)
CONST = (1/r^2) * 1/(4*pi);
% *********************************************************

for i = 1:size(COS,1),
  for j = 1:size(COS,2),
    
    % Px is an Nx1 column array, starting at P0x (we only need P1x - PNx)
    Px = LegendreP(N,COS(i,j));
    
    Hx(i,j) = CONST * dot( Series, Px(2:N+1) );
    
  end
end

t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return
