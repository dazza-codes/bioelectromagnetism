function [FV,C] = eeg_interp_sph_spline_scd(Zi,Ei)

% eeg_interp_sph_spline_scd - Spherical Spline Scalp Current Density
%
% Useage: [FV,C] = eeg_interp_sph_spline_scd(Zi,Ei)
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

% $Revision: 1.2 $ $Date: 2008/05/04 18:23:43 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2003 Darren.Weber_at_radiology.ucsf.edu, with
%                   adapted from eeg_interp_sph_spline
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COS = elec_cosines(Ei,Ei);
G = eeg_interp_sph_spline_g(COS);
C = eeg_interp_sph_spline_c(Zi,G);

Co = C(1);
Ci = C(2:end);

eegversion = '$Revision: 1.2 $';
fprintf('EEG_INTERP_SPH_SPLINE_SCD [v %s]\n',eegversion(11:15)); tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that C is solved, we can obtain interpolated scalp current density
% at Ej (Eq. 5)

r = 10.0;

% create spherical interpolation positions
fprintf('...generating spherical interpolation points\n');
FV = sphere_tri('ico',4,r);

% cosines between measurement electrodes and interpolation points
Cos = elec_cosines(Ei,FV.vertices);

% Calculate h(x)
%Hx = eeg_interp_sph_spline_h(Cos,r);  % nElectrodes x NinterpolationPoints
Hx = eeg_interp_sph_spline_h(Cos,1);  % nElectrodes x NinterpolationPoints

% Solve Eq 5.
FV.Cdata = (Ci' * Hx)';

t=toc; fprintf('...done (%6.2f sec)\n',t);

return
