function [V] = eeg_forward(rd,qd,radius,cond,el,Nspheres,Nchans,Ndipoles)
% EEG_FORWARD: forward calculation of electric potential
%
% Calculates the potential, given dipole positions and moments,
% electrode positions and the radii and conductivity of BEM spheres
%
% Useage: [V] = eeg_forward(rd,qd,radius,cond,el,Nspheres,Nchans,Ndipoles)
%
% rd        - vector[Ndipoles*3] with position dipoles
% qd        - vector[Ndipoles*3] with dipole moment
% radius    - vector[Nspheres] with the radii of the spheres
% cond      - vector[Nspheres] with the conductivity of the spheres
% el        - matrix[Nchans,3] with the Cartesian electrode coordinates
% V         - vector[Nchans] returning the calculated potential
% Nspheres  - number of spheres
% Nchans    - number of channels
% Ndipoles  - number of dipoles
%

% Licence:  GNU GPL, no implied or express warranties
% History:  (c) 1998 Robert Oostenveld, available under GNU licence
%               - created as a c program, part of the MSA library
%           04/2002, Darren.Weber@flinders.edu.au
%                    - attempting to adapt c to matlab code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (Ndipoles<1), error('Ndipoles < 1\n'); end

% Check that all dipoles are inside the innermost sphere
radius_rd = sqrt( sum (rd.^2,2) );
if ( max(radius_rd) > max(radius(1) ),
    error('A dipole is outside the inner sphere.\n\n');
end



% calculate the leadfield for a dipole in an homogenous sphere

% If there is only one dipole, just do it.
% If there are more dipoles, calculate the leadfield for each
% single dipole and copy-and-paste them together in one large matrix

lf = matrix(1,Nchans,1,Ndipoles*3);

if isequal(Ndipoles,1),
    fprintf('electric_potential: calculating leadfield for one dipole\n');
    lf = eeg_leadfield_sphere(rd, radius, cond, el, Nspheres, Nchans);
else
    % use a temporary matrix m to store the leadfield for each single dipole
    
    m = zeros(Nchans,3);
    
    for i=1:Ndipoles,
        
        fprintf('electric_potential: calculating leadfield for dipole %i\n', i);
        
        m = eeg_leadfield_sphere(rd3*(i-1), radius, cond, el, Nspheres, Nchans);
        
        lf(i) = m;
    end
end


% multiply leadfield with the dipole moment, to obtain the potential
for i=1:Nchans,
    
    % !!!!!!!!!!!!!!
    V(i) = lf(i) * qd(i);
    
    fprintf('electric_potential: %g\n', V(i));
end

% take care that the calculated potential is average referenced
avg = sum(V) / Nchans;

fprintf('electric_potential: subtracting average potential %g\n', avg);

V = V - avg;

return
