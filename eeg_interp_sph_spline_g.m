function [Gx] = eeg_interp_sph_spline_g(COS)

% [Gx] = eeg_interp_sph_spline_g(COS)
%
% COS is the cosine matrix from elec_cosines
%
% Gx is the solution to Eq. 3 of Perrin et al. (1989)
%
% g(x) = 1/(4*pi) * (for n=1:inf, sum = sum + ( ( (2*n+1)/(n^m * (n+1)^m) ) * Pn(x) ) );
%
% where x is the cosine between two points on a sphere,
% m is a constant > 1 and Pn(x) is the nth degree Legendre
% polynomial.  Perrin et al. (1989) evaluated m=1:6 and 
% recommend m=4, for which the first 7 terms of Pn(x) 
% are sufficient to obtain a precision of 10^-6 for g(x)
%
% Notes:    This function solves g(x), Eq. 3 in
%           Perrin et al (1989).  Electroenceph. & Clin. 
%             Neurophysiology, 72: 184-187. Corrigenda (1990),
%             Electroenceph. & Clin. Neurophysiology, 76: 565.
%             (see comments in the .m file for details).
%
% see elec_cosines, eeg_interp_sph_spline, LegendreP

% $Revision: 1.3 $ $Date: 2007/02/16 01:24:26 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2001 Darren.Weber_at_radiology.ucsf.edu, with
%                   mathematical advice from
%                   Dr. Murk Bottema (Flinders University of SA)
%           10/2003 Darren.Weber_at_radiology.ucsf.edu, with
%                   mathematical advice and LegendreP function from 
%                   Dr. Tom.Ferree_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegversion = '$Revision: 1.3 $';
fprintf('EEG_INTERP_SPH_SPLINE_G [v %s]\n',eegversion(11:15)); tic
fprintf('...calculating Legendre function of the cosine matrix\n');

% see comments at end by Robert Oostenveld
m = 4;
N = 7;

n = [1:N]';
Series = (2*n + 1) ./ (n.^m .* (n+1).^m);
%fprintf('%12.10f  ',Series) % note how Series diminishes quickly
%0.1875000000  0.0038580247  0.0003375772  0.0000562500  
%0.0000135802  0.0000041778  0.0000015252


% Perrin et al. (1989) recommend tabulation of g(x) for
% x = linspace(-1,1,2000) to be used as a lookup given actual
% values for cos(Ei,Ej).
if isempty(COS),
  error('...cosine matrix is empty.\n');
  %msg = sprintf('...cosine matrix empty, generating COS.\n');
  %warning(msg);
  %COS = linspace(-1,1,2000)';
end

Gx = zeros(size(COS));

CONST = 1/(4*pi);

for i = 1:size(COS,1),
  for j = 1:size(COS,2),
    
    % P7x is an 8x1 column array, starting at P0x (we only need P1x - P7x)
    P7x = LegendreP(N,COS(i,j));
    
    Gx(i,j) = CONST * dot( Series, P7x(2:N+1) );
    
  end
end

t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return



% --------------------------------------------------
% Comments from Dr. Robert Oostenveld

%I have an independent implementation of the Perrin 1989 algorithm for
%spherical splines, and have encountered the same problem as you are
%describing (singularity in the inverse of Gx). The algorithm has some built
%in constants (m = 4; N = 7;) that are "optimal" according to Perrin. But
%these values were selected for a small number of electrodes (probably not
%more than 32 at t hat time) and do not apply here. I have tried to solve
%this in my own implementation with an adaptive N:

% N=9 works fine for 32 electrodes, but for 128 electrodes it should be larger
if max(size(x))<=32
  N = 9;
elseif max(size(x))<=64
  N = 14;
elseif max(size(x))<=128
  N = 20;
else
  N = 32;
end

%where size(x) corresponds to the number of electrodes. This considerably
%improved the results if the data was not too noisy. But I still did not get
%satisfactory results for noisy 128 EEG data (and it also slows down the
%computation). You could give it a try nevertheless. Also changing the
%backslash (solve by Gaussian elimination) into a regularised pinv (solve
%with SVD) could improve your results. You should have a look at the (test?)
%code that Darren has included the function and that is commented out.
%
%best regards,
%Robert
%
%----------------------------------------------------------------------
%Robert Oostenveld, PhD
%Center for Sensory-Motor Interaction (SMI)
%Aalborg University, Denmark
%
%and
%
%F.C. Donders Centre for Cognitive Neuroimaging
%University Nijmegen
%P.O. Box 9101
%NL-6500 AH Nijmegen
%The Netherlands
%
%Tel: +31 (0)24 3619695
%Fax: +31 (0)24 3610989
%
