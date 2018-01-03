function [lf] = eeg_leadfield_sphere(rd, radius, conductivity, el, Nspheres, Nchans);
% EEG_LEADFIELD_SPHERE: electric leadfield matrix for a single dipole
%
% Calculates the electric leadfield from a single current dipole
% in a homogenous sphere or in multiple concentric spheres
%
% Useage: [lf] = eeg_leadfield_sphere(rd,radius,cond,el,Nspheres,Nchans)
%
% rd        - vector[1..3] position dipole
% radius    - vector[Nspheres] with the radii of the spheres
% cond      - vector[Nspheres] with the conductivity of the spheres
% el        - matrix[Nchans,3] with the Cartesian electrode coordinates
% V         - vector[Nchans] returning the calculated potential
% Nspheres  - number of spheres
% Nchans    - number of channels
% Ndipoles  - number of dipoles
%
% Notes: radius should increase from the inner to the
%        outer sphere.
%        cond should contain the conductivity values for
%        brain, CSF, skull & scalp (in the 4 shell model)
%        or brain, skull & scalp (in the 3 shell model).
%        Common values are:
%
%        SCALP 12 mm @ 0.3300 S/m
%        SKULL 10 mm @ 0.0042 S/m
%        BRAIN/CSF 08 mm @ 0.3300 S/m
%

% Licence:  GNU GPL, no implied or express warranties
% History:  (c) 1998 Robert Oostenveld, available under GNU licence
%               - created as a c program, part of the MSA library
%           04/2002, Darren.Weber@flinders.edu.au
%                    - attempting to adapt c to matlab code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default radius/conductivity
r = [0.0800 0.0850 0.0900 0.1000]; % radius in meters
                                   % (scalp 10cm, skull 9cm, CSF 8.5cm, brain 8cm)
c = [0.3300 0.3300 0.0042 0.3300]; % conducitivty in Siemens/meter

if (Nspheres==1),
    
    % here one could check if the dipole lies on the z-axis and
    % then take the more flexible but inefficient 4 sphere calculation
    lf = leadfield_sphere1(rd, radius(1), conductivity(1), electrode, Nchans);
    
elseif (Nspheres==2),
    
    % the leadfield is calculated using the algorithm for four
    % concentric spheres, whereby the spheres for the CSF
    % brain and skull are combined to a single sphere
    
    r(1) = radius(1);	% brain
    r(2) = radius(1);	% CSF is combined with brain
    r(3) = radius(1);	% skull
    r(4) = radius(2);	% scalp
    
    c(1) = conductivity(1);	% brain
    c(2) = conductivity(1);	% CSF is combined with brain
    c(3) = conductivity(1);	% skull
    c(4) = conductivity(2);	% scalp
    
    lf = leadfield_sphere4(rd, r, c, electrode, Nchans);
    
elseif (Nspheres==3),
    
    % the leadfield is calculated using the algorithm for four
    % concentric spheres, whereby the spheres for the CSF
    % and the brain are combined to a single sphere
    
    r(1) = radius(1);	% brain
    r(2) = radius(1);	% CSF is combined with brain
    r(3) = radius(2);	% skull
    r(4) = radius(3);	% scalp
    
    c(1) = conductivity(1); % brain
    c(2) = conductivity(1); % CSF is combined with brain
    c(3) = conductivity(2); % skull
    c(4) = conductivity(3); % scalp
    
    lf = leadfield_sphere4(rd, r, c, electrode, Nchans);
    
elseif (Nspheres==4),
  
  lf = leadfield_sphere4(rd, radius, conductivity, electrode, Nchans);
  
else
  error('No more than 4 spheres can be modelled\n\n');
end

return;






function [lf] = leadfield_sphere1(rd, radius, conductivity, electrode, Nchans)
% LEADFIELD_SPHERE1: electric leadfield matrix for a single dipole in a homogenous sphere
%
% Calculates the electric leadfield from a single current dipole 
% in a homogenous sphere for a set of electrodes on this sphere
%
% Useage: (lf) = leadfield_sphere1(rd,radius,cond,el,Nchans)
%
% rd        - vector[1:3] position of dipole in xyz
% radius    - radius of the sphere
% cond      - conductivity of the sphere
% el        - matrix[Nchans,3] with the Cartesian electrode coordinates
% Nchans    - number of channels
%
% lf        - leadfield matrix returned
%
% warning: 
%   the algorithm does not work for a dipole in the origin
%   and it does not work for an electrode on the z-axis 
%
% reference: Luetkenhoener, Habilschrift '92
%
%

  if (Nchans<=1), error('Nchans <= 1\n'); end
  if (radius<=0), error('radius <= 0\n'); end
  if (conductivity<=0), error('conductivity <=0\n'); end
  
  sl = norm(rd);

  if (s1 < 0.0001),
      
      % to avoid this annoying exception, shift the dipole to a position nearby
      % this should be taken care of more properly!!
      fprintf('warning: dipole is in origin \n');
      rd = [0 0 0.0001];
  end
  
  for i=1:Nchans,
      
      sl = norm(electrode(i));
      
      if (abs(s1-radius)/radius > 0.01),
          fprintf('warning: electrode %i not on sphere\n', i);
      end
      
      R    = rd; 		% position dipole
      dr   = ones(1,3);
      rxR  = ones(1,3);
      vec1 = ones(1,3);
      vec2 = ones(1,3);
      
      for i=1:Nchans,
          
          fprintf('leadfield_sphere1: calculating electrode %i\n', i);
          
          r = electrode(i);	% position electrode (carthesian coordinates)
          
          abs_r = norm(r);
          
          SUB_VEC(&(r(1)), &(R(1)), &(dr(1)), 3, float, float, float);
          CROSSP_VEC(&(r(1)), &(R(1)), &(rxR(1)), float, float, float);
          
          fprintf('R  = %f %f %f\n', R(1), R(2), R(3));
          fprintf('dr = %f %f %f\n', dr(1), dr(2), dr(3));
          
          NORM_VEC(&(R(1)), abs_R, 3, float);
          NORM_VEC(&(dr(1)), abs_dr, 3, float);
          NORM_VEC(&(rxR(1)), abs_rxR, 3, float);
          
          s1 = abs_dr*abs_dr*abs_dr;
          s2 = abs_R*abs_R;
          
          fprintf('s1 = %f, s2 = %f\n', s1, s2);
          fprintf('abs_rxR = %f\n', abs_rxR);
          
          SCALE_VEC(&(r(1)), &(vec1(1)), 1.0/abs_r, 3, float, float);
          SCALE_VEC(&(dr(1)), &(vec2(1)), 1.0/abs_dr, 3, float, float);
          SUB1_VEC(&(vec1(1)), &(vec2(1)), 3, float, float);
          DOTP_VEC(&(vec1(1)), &(R(1)), A, 3, float, float);
          
          fprintf('A = %f\n', A);
          A /= abs_rxR*abs_rxR;
          fprintf('A = %f\n', A);
          A += 2/s1;
          fprintf('A = %f\n', A);
          
          B  = (abs_r*abs_r - s2)/s1 - 1/abs_r;
          fprintf('A = %f, B = %f\n', A, B);
          
          CROSSP_VEC(&(R(1)), &(rxR(1)), &(vec1(1)), float, float, float);
          SCALE_VEC(&(vec1(1)), &(vec1(1)), A, 3, float, float);
          SCALE_VEC(&(R(1)), &(vec2(1)), B, 3, float, float);
          
          ADD1_VEC(&(vec1(1)), &(vec2(1)), 3, float, float);
          SCALE_VEC(&(vec1(1)), &(lf(i,1)), 1/(4*M_PI*conductivity*s2), 3, float, float);
      end
      
      fprintf('leadfield_sphere1: ready\n');
  end
return



function [lf] = leadfield_sphere4(rd, radius, conductivity, electrode, Nchans)
% LEADFIELD_SPHERE4: electric leadfield matrix for a single dipole in 4 concentric spheres
% 
% Calculates the electric leadfield from a single current dipole 
% in a homogenous sphere for a set of electrodes on this sphere
% 
% Useage: [lf] = leadfield_sphere4(rd,radius,cond,el,Nchans)
%
% rd        - vector[1:3] position of dipole in xyz
% radius    - vector[1:4] radius of the spheres
% cond      - vector[1:4] conductivity of the spheres
% el        - matrix[Nchans,3] with the Cartesian electrode coordinates
% Nchans    - number of channels
%
% lf        - leadfield matrix returned
%
% warning: 
%   the algorithm does not work for a dipole in the origin
%   and it does not work for an electrode on the z-axis 
%
% reference: Luetkenhoener, Habilschrift '92
%

if (Nchans<=1), error('Nchans <= 1\n'); end

% the dipole and electrode positions are rotated, such that the dipole
% position is on the z-axis (0,0,R). After calculation of the leadfield
% matrix, this matrix undergoes the backward coordinate transformation

% this is not the most efficient: it would be easier (and faster) to leave
% leadfield matrix as it is and rotate the dipole moment, which is only a
% vector instead of a matrix, before the potential is calculated

fprintf('leadfield_sphere4: calculating rotation matrix (%f,%f,%f)\n', rd(1), rd(2), rd(3));

R = norm(rd);

if (R <= 0.0001),
    
    % to avoid this annoying exception, shift the dipole to a position nearby
    % this should be taken care of more properly!!
    fprintf('warning: dipole is in origin \n');
    rd = [0 0 0.0001];
end

% for a dipole in the origin I currently (981112) believe that the
% leadfield is the same to a single sphere leadfield, so that function
% could be called here instead (but it regretfully also does not work
% for a dipole in the origin currently)

A   = zeros(3);	% rotation matrix
Atr = A';     	% inverse rotation matrix, which equals transposed matrix

val = sqrt( sum( rd(1:2).^2, 2) );

if and(val < 0.0001, rd(3)>0),
    fprintf('leadfield_sphere4: no rotation of dipole towards z-axis necessary\n');
    A = eye(3);
    Atr = A';
    
elseif and(val < 0.0001, rd(3)<0),
    fprintf('leadfield_sphere4: 180 deg. rotation of dipole towards z-axis necessary\n');
    % rotate around x-axis
    A = [1 0 0; 0 -1 0; 0 0 -1];
    Atr = A';
    
else
    fprintf('leadfield_sphere4: rotating dipole towards z-axis\n');
    
    A(1,1) = rd(1) * rd(3) / (R * val);
    A(1,2) = rd(2) * rd(3) / (R * val);
    A(1,3) = -1.0  * val   /  R;
    
    A(2,1) = -1.0 * rd(2) / val;
    A(2,2) =        rd(1) / val;
    A(2,3) =                  0;
    
    A(3,1:3) = rd ./ R;
    
    Atr = A';
end

fprintf('leadfield_sphere4: A =\n %f %f %f\n%f %f %f\n%f %f %f\n', A);

excentricity = R / radius(4);

fprintf('leadfield_sphere4: excentricity = %f\n', excentricity);

for i=1:Nchans,
    
    fprintf('leadfield_sphere4: calculating leadfield for electrode %i\n', i);
    
    r = norm(electrode(i,:));
    if (abs(r-radius(4))/radius(4) > 0.5),
        fprintf('warning: electrode %d not on outer sphere (%f, %f)\n', i, r, radius(4));
    end
    
    lf(i) = zeros(1,3);
    
    % rotate electrode position using matrix multiplication (x,y,z)' = A * (x,y,z)
    sph = [A * electrode(i)']';
    
    % calculate the position of the electrode in conventional spherical coordinates:
    % (r, phi, theta) =  (radius, angle from x-axis, angle from z-axis)
    %
    % an alternative spherical coordinate system which is also used in the msa library 
    % for example in sph2cart and cart2sph is:
    % (az, el, r) = (azimult/angle from x-axis, elevation/angle from xy-plane, radius)
    
    
    
    ****************************************************
    ******* Not sure that the original cart2sph creates 
    ******* the same values as the matlab
    ****************************************************
    [theta, phi, r] = cart2sph(sph);
    %r   = sph(3);
    %phi   = sph(1);
    %theta = M_PI_2 - sph(2);
    
    
    
    
    
    
    
    cos_theta  = cos(theta);
    
    fprintf('leadfield_sphere4: r = %f, phi = %f, theta = %f\n', r, phi, theta);
    fprintf('leadfield_sphere4: theta = %g, cos_theta = %g\n', theta, cos_theta);
    
    % calculate summation term of x and y-component of leadfield
    j = 0;
    n = 0;
    
    while and((n < plgndr_ntot),(j < plgndr_nval)),
        
        n = n + 1;
        
        % of course this could be much more efficient...
        % but it is nice enough that it works!
        
        exp = 2*n+1;
        o12 = conductivity(1)/conductivity(2);
        o23 = conductivity(2)/conductivity(3);
        o34 = conductivity(3)/conductivity(4);
        
        C = ((n*o12+n+1)*(n*o23+n+1)+n*(n+1)*(o12-1)*(o23-1)*pow(radius(1)/radius(2),exp))* ...
			((n*o34+n+1)+(n+1)*(o34-1)*(pow(radius(3)/radius(4),exp))) ...
			+((o12-1)*((n+1)*o23+n)*pow(radius(1)/radius(3),2*n+1)+(n*o12+n+1)*(o23-1)*pow(radius(2)/radius(3),exp))* ...
			(n+1)*(n*(o34-1)+((n+1)*o34+n)*pow(radius(3)/radius(4),exp));
        
        val = (pow(2*n+1, 4) / (n*C)) * pow(excentricity, n-1) * plgndr(n, 1, cos_theta);
        
        lf(i,1) = lf(i,1) + val;
        
        if or(abs(val)<eps, abs(val/lf(i,1)) < plgndr_accuracy),
            j = j + 1;
        else
            j = 0;
        end
        fprintf('leadfield_sphere4: n = %i, j = %i, C = %g, val = %g, lf = %g\n', n, j, C, val, lf(i,1));
    end
    
    fprintf('leadfield_sphere4: %i terms needed for convergence of x/y-val\n', n);
    
    % calculate summation term of z-component of leadfield
    j = 0;
    n = 0;
    
    while and((n < plgndr_ntot),(j < plgndr_nval)),
        
        n = n + 1;
        
        % of course this could be much more efficient...
        % but it is nice enough that it works!
        
        exp = 2*n+1;
        o12 = conductivity(1)/conductivity(2);
        o23 = conductivity(2)/conductivity(3);
        o34 = conductivity(3)/conductivity(4);
        
        C = ((n*o12+n+1)*(n*o23+n+1)+n*(n+1)*(o12-1)*(o23-1)*pow(radius(1)/radius(2),exp))* ...
			((n*o34+n+1)+(n+1)*(o34-1)*(pow(radius(3)/radius(4),exp))) ...
			+((o12-1)*((n+1)*o23+n)*pow(radius(1)/radius(3),2*n+1)+(n*o12+n+1)*(o23-1)*pow(radius(2)/radius(3),exp))* 
			(n+1)*(n*(o34-1)+((n+1)*o34+n)*pow(radius(3)/radius(4),exp));
        
        val = (pow(2*n+1, 4) / C) * pow(excentricity, n-1) * plgndr(n, 0, cos_theta);
        
        lf(i,3) = lf(i,3) + val;
        
        if or(abs(val)<eps, abs(val/lf(i,1))<plgndr_accuracy),
            j = j + 1;
        else
            j = 0;
        end
        fprintf('leadfield_sphere4: n = %i, j = %i, C = %g, val = %g, lf = %g\n', n, j, C, val, lf(i,3));
    end
    
    fprintf('leadfield_sphere4: %i terms needed for convergence of z-val\n', n);
    
    lf(i,3) = lf(i,3) / 4*M_PI*conductivity(4)*radius(4)*radius(4);
    lf(i,2) = sin(phi)/(4*M_PI*conductivity(4)*radius(4)*radius(4)) * lf(i,1);
    lf(i,1) = lf(i,1) * cos(phi)/(4*M_PI*conductivity(4)*radius(4)*radius(4));
    
    % this is strange -> signs need correction
    lf(i,:) = lf(i,:) * -1;
    
    % rotate leadfield back using matrix multiplication (x,y,z) = Atr * (x,y,z)'
    fprintf('leadfield_sphere4: lf before rotation = %f %f %f\n', lf(i,:));
    lf(i,:) = [Atr * lf(i,:)']';
    fprintf('leadfield_sphere4: lf after rotation  = %f %f %f\n', lf(i,:));
end

% leadfield has been calculated for this electrode
fprintf('leadfield_sphere4: ready\n');
return;
