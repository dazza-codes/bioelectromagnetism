function [out1,out2,out3] = bvpdemo3(t,y,flag,E0,delta)

%   definition of an ODE-system used by BVPDEMO
%
%   See BVP BVPDEMO


if nargin < 3 | isempty(flag)    % Return dy/dt = F(t,y).
    
  n=y(1);
  E=y(2);
  
  v=-E./(E+1).^2;
  c=-E0./(E0+1).^2;

  
  out1 = [(v-c).*n ; n-1];

    
else
   switch(flag)
   case 'init'                       % Return default [tspan,y0,options].
      
      c=-E0./(E0+1).^2;
      vprime=(-1+E0)./(1+E0).^3;

      
      out1 = [0 200];                            % < Insert tspan here. >;
      out2 = [1;E0] + delta*[-sqrt(vprime);-1];  % Start auf instabiler MF
      %out3 = odeset('stats','on','events','on','OutputFcn','Odephas2',...
      %              'AbsTol',1e-8,'RelTol',1e-5);  
      out3 = odeset('events','on','AbsTol',1e-8,'RelTol',1e-5);  

      
      
   case 'events'                     % Return event vector and information.
      distvec=y-[1;E0];
      
      out1 = sqrt(distvec'*distvec) - delta;   %< Insert event function vector here. >
      out2 = 1;                %< Insert logical ISTERMINAL vector here.>;
      out3 = -1;               %< Insert DIRECTION vector here.>;
          
   otherwise
      error(['Unknown flag ''' flag '''.']);
   end
end
