function FV = source_dipole_sphere(dipole_moment,dipole_origin,dipole_orient,sphere_radius,LMAX)

% source_dipole_sphere - calculate surface potential on a sphere containing a dipole
%
% Calculates surface voltage, for a single sphere containing a radial
% dipole.  This function is limited to radial dipoles (a dipole at the
% origin can have any orientation).
%
% The sphere is assigned a homogenous, isotropic conductivity, 
% sigma = 0.33 Siemens/meter, which is commonly used to estimate brain
% conductivity.
%
% Also, the distance of dipole source to sink is 1e-12 (1 pico meter).
%
% FV = source_dipole_sphere([dipole_moment],[dipole_origin],...
%                           [dipole_orient],[sphere_radius],[Lmax]);
% -------------------------------------------------------------------------
% Inputs
% 
% dipole_moment - Ampere meters, 10e-9 by default, which approximates a
% small patch of cortex (see Hamalainen et al [1993], Reviews of Modern
% Physics, 65(2):413-497).
%
% dipole_origin - The XYZ Cartesian coordinates of the midpoint between the
% source and sink of the dipole.  The orientation of the dipole will be
% radial, with source and sink arranged either side of midpoint.  See also
% the sphere radius specification below.
%
% dipole_orient - Only if the dipole_origin = [0,0,0] (the default) is it
% possible to specify the dipole orientation; eg, [1 0 0] is along the
% x-axis, [0 1 0] is along the y-axis, etc.  If the dipole is located away
% from the origin, the orientation is along the radial line from the origin
% to the dipole.
%
% sphere_radius - in meters; the default 0.1 approximates a human head.
% Note that any specification of dipole origin should be within the sphere
% radius.
% 
% Lmax - the order of magnitude for the sum of Legendres, which are used
% for the analytic solution of potential at a spherical surface.  The
% default value is 20, determined after examination of values from 1:250
% (20 is conservative, as convergence was observed at Lmax < 5).
%
% -------------------------------------------------------------------------
% Returns the FV struct, a spherical tesselation
%
% FV = 
%    vertices: [2562x3 double]
%       faces: [5120x3 double]
%       Cdata: [2562x1 double] - surface potential (Volts)
%
% To visualize the result, use:
%
% patch('vertices',FV.vertices,'faces',FV.faces,...
%       'FaceVertexCData',FV.Cdata,...
%       'Facecolor','interp');
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:36 $

% Licence:  GNU GPL, no implied or express warranties
% History:  12/2003, Darren.Weber_at_radiology.ucsf.edu
%                    developed from notes and advice from 
%                    Dr. Tom.Ferree_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% LMAX is a parameter used in the Legendre approximation
if ~exist('LMAX','var'), LMAX = 20; end
if isempty(LMAX),        LMAX = 20; end

% create a spherical surface with a radius that approximates the human
% head, 0.1 meters (10 cm)
if ~exist('sphere_radius','var'), sphere_radius = 0.1; end
if isempty(sphere_radius),        sphere_radius = 0.1; end

% sphere_tri is an external function
FV = sphere_tri('ico',4,sphere_radius);

% assign a homogenous, isotropic conductivity to this spherical region
sigma = 0.33; % Siemens/meter, common brain value
fprintf('...assigned brain conductivity (0.33 S/m) to volume\n');

% distance of dipole source to sink, could approximate the cortical
% thickness
dipole_distance = 1e-12; % 1 pico meter
fprintf('...dipole has source and sink %g meters apart\n',dipole_distance);

% dipole moment, Ampere meters (10 nA m)
if ~exist('dipole_moment','var'), dipole_moment = 10e-9; end
if isempty(dipole_moment),        dipole_moment = 10e-9; end
fprintf('...dipole moment = %g Ampere meters\n',dipole_moment);

% dipole current, Amperes
dipole_current = dipole_moment / dipole_distance;
fprintf('...dipole current = %g Amperes\n',dipole_current);

% Define the vector that lies at the midpoint between the source and sink.
% The origin is the default.
if ~exist('dipole_origin','var'), dipole_origin = [0,0,0]; end
if isempty(dipole_origin),        dipole_origin = [0,0,0]; end
fprintf('...dipole origin = [ %g, %g, %g ] meters\n',dipole_origin);

% find the magnitude (ie, radius) of the dipole (do not use the matlab norm
% function, it doesn't calculate the same value)
dipole_radius = sqrt( sum( dipole_origin .^2 ) );

% Calculate the radial dipole unit vector, which is used to find the cosine
% of the angle between the dipole and surface locations below
if dipole_radius,
  
  % the dipole is not located at the origin
  dipole_unit_vector = dipole_origin / dipole_radius;
else
  
  % the dipole is located at the origin, so the radial orientation is not
  % clearly defined and we need another input parameter to define the
  % dipole unit vector
  
  % dipole geometry - oriented along z axis
  if ~exist('dipole_orient','var'), dipole_orient = [0 0 1]; end
  if isempty(dipole_orient),        dipole_orient = [0 0 1]; end
  
  % now calculate the dipole unit vector
  dipole_unit_vector = dipole_orient / sqrt( sum( dipole_orient .^ 2 ));
  
end
fprintf('...dipole unit vector = [ %g, %g, %g ]\n',dipole_unit_vector);

% calculate location of source and sink, just for fun!
source = dipole_origin + (dipole_unit_vector * dipole_distance);
sink   = dipole_origin - (dipole_unit_vector * dipole_distance);
fprintf('...source vector = [ %+g, %+g, %+g ] meters\n', source);
fprintf('...sink vector   = [ %+g, %+g, %+g ] meters\n', sink);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve voltage on surface
  
fprintf('...calculating surface voltage...'); tic

% initialise the surface potentials
FV.Cdata = zeros(length(FV.vertices),1);

if dipole_radius,
  
  DIPOLE_CONSTANT = ( -dipole_moment / (4 * pi * sigma) ) * (1 / (dipole_radius^2) );
  
  for vertex = 1:length(FV.vertices),
    
    % extract spherical surface coordinates
    surface_vector = FV.vertices(vertex,:);
    
    surface_radius = sqrt( sum ( surface_vector .^2 ) );
    
    if surface_radius < dipole_radius,
      error('surface radius is smaller than dipole radius.');
    end
    
    if surface_radius,
      surface_unit_vector = surface_vector / surface_radius;
    else
      surface_unit_vector = [0 0 0];
      warning('surface radius is zero.');
    end
    
    cos_theta = dot(dipole_unit_vector, surface_unit_vector);
    
    summation = 0;
    
    for L = 1:LMAX,
      
      TERM1 = 2 * L + 1;
      
      TERM2 = ( dipole_radius / surface_radius )^(L + 1);
      
      PnX = LegendreP(L,cos_theta);
      if length(PnX) > 1,
        PnX_sum = sum(PnX(2:end));
      else
        PnX_sum = 0;
      end
      summation = summation + ( TERM1 * TERM2 * PnX_sum );
      
    end
    
    FV.Cdata(vertex) = DIPOLE_CONSTANT * summation;
    
  end % for vertex
  
else
  
  % special case where dipole is at the origin.
  
  DIPOLE_CONSTANT = ( -dipole_moment / (4 * pi * sigma) );
  
  for vertex = 1:length(FV.vertices),
    
    surface_vector = FV.vertices(vertex,:);
    
    surface_radius = sqrt( sum ( surface_vector .^2 ) );
    
    if surface_radius < dipole_radius,
      error('surface radius is smaller than dipole radius.');
    end
    
    if surface_radius,
      surface_unit_vector = surface_vector / surface_radius;
    else
      surface_unit_vector = [0 0 0];
      warning('surface radius is zero.');
    end
    
    cos_theta = dot(dipole_unit_vector, surface_unit_vector);
    
    FV.Cdata(vertex) = DIPOLE_CONSTANT * (3 / (surface_radius^2) ) * cos_theta;
    
  end % for vertex
  
end

t=toc; fprintf('done (%5.2f sec)\n',t);

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS...

function p=LegendreP(nmax,x)

% function p=LegendreP(nmax,x)
% Computes all Legendre polynomials from n=0 to n=nmax, 
% them in an array from n=1 to nmax+1
% Uses the recursion relation in Numerical Recipes.
% T Ferree

p=zeros(nmax+1,1);

if nmax>=0
  p0=1.0;
  p(1)=p0;
end

if nmax>=1 
  p1=x;
  p(2)=p1;
end

if nmax>=2
  for j=2:nmax
    c1=(2.0*j-1.0)/(1.0*j);
    c2=(1.0*j-1.0)/(1.0*j);
    p(j+1)=c1*x*p(j-1+1)-c2*p(j-2+1);
  end
end

return;
