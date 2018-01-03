function [FV,C] = eeg_interp_sph_spline_scd_sig(Zi,Ei)
%
% SIGNE: modified this function to calculate the scd only at the
% 	 electrode locations
% eeg_interp_sph_spline_scd - Spherical Spline Scalp Current Density
%
% Useage:   [FV,C] = eeg_interp_sph_spline_scd(Zi,Ei)
%
% This function calls eeg_interp_sph_spline,
% and replaces the relevant sections from Perrin et al. (1990)
% to calculate the scalp current density (SCD).
%
% FV => interpolated spherical surface
%
% FV.faces    => triangulation of FV.vertices 
% FV.vertices => cartesian coordinates (Nx3)
% FV.Cdata    => spherical spline SCD at FV.vertices
% 
% C => interpolation coefficients of Ei (includes co = C(1))
% 
% Notes:    This function calculates the spherical spline SCD of 
%           Perrin et al (1989).  Electroenceph. & Clin. 
%             Neurophysiology, 72: 184-187. Corrigenda (1990),
%             Electroenceph. & Clin. Neurophysiology, 76: 565.
%             (see comments in the .m file for details).

% $Revision: 1.1.1.1 $ $Date: 2007/11/26 23:01:43 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2003 Darren.Weber_at_radiology.ucsf.edu, with
%                   adapted from eeg_interp_sph_spline
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

keyboard

COS = elec_cosines(Ei,Ei); % where Ei is Nx3 [X Y Z]
G = eeg_interp_sph_spline_g(COS);
C = eeg_interp_sph_spline_c(Zi,G);

Co = C(1);
Ci = C(2:end);

eegversion = '$Revision: 1.1.1.1 $';
fprintf('EEG_INTERP_SPH_SPLINE_SCD [v %s]\n',eegversion(11:15)); tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that C is solved, we can obtain interpolated scalp current density
% at Ej (Eq. 5)

% create spherical interpolation positions
%fprintf('...generating spherical interpolation points\n');
%FV = sphere_tri('ico',4,r);

% cosines between measurement electrodes and interpolation points
%Cos = elec_cosines(Ei,Ei);%FV.vertices);

% Calculate h(x)
Hx  = h(COS,r);  % nElectrodes x NinterpolationPoints

% Solve Eq 5.
FV.Cdata = Ci' * Hx;

t=toc; fprintf('...done (%6.2f sec)\n',t);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Hx] = h(X,r)
% Solve h(x) in Eq. 5 Perrin et al. (1989)
%
% h(x) = (1/r^2) * 1/(4*pi) *
%  for n=1:inf, sum = sum + ( ( (2*n+1)/(n^(m-1) * (n+1)^(m-1)) ) * Pn(x) );
% where x is the cosine between two points on a sphere,
% m is a constant > 1 and Pn(x) is the nth degree Legendre
% polynomial.  Perrin et al. (1989) evaluated m=1:6 and 
% recommend m=4, for which the first 20 terms of Pn(x) 
% are sufficient to obtain a precision of 10^-6 for h(x)
%
% Perrin et al. (1989) assume r = 1, so we adjust the eq. above Eq. 5:
% del_perp^2[Pn(x)] = [ -n*(n+1)/r^2 ] Pn(x)

fprintf('...calculating h(x)...'); tic

m = 4;
N = 20;

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
if isempty(X),
  error('...cosine matrix is empty.\n');
  %msg = sprintf('...cosine matrix empty, generating X.\n');
  %warning(msg);
  %X = linspace(-1,1,2000)';
end

Hx = zeros(size(X));

for i = 1:size(X,1),
  for j = 1:size(X,2),
    
    % P20x is a 20x1 column array, starting at P0x
    P20x = LegendreP(N,X(i,j));  % see function below
    
    Hx(i,j) = (1/(4*pi)) * ( Series' * P20x(2:end) );
    
  end
end

t=toc; fprintf('...done (%6.2f sec)\n',t);

return
