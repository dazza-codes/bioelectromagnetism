function bvpdemo
%
%  BVPDEMO demonstrates the usage of BVP and BVP2
%
%  See also BVP  BVP2  BVPSET  BVPPLOT  BVPPHAS2  BVPPHAS3

%  Copyright (c) 1999 Guenter Kneisl
%                     University of Technology, Vienna
%                     e9425595@fbma.tuwien.ac.at

clc;home;
fprintf('\n BVP and BVP2 are tools for solving nonlinear boundary \n');
fprintf(' value problems (of first, and second order, respectively)\n');
fprintf(' of arbitrary dimension on nonequidistant grids.\n\n');
fprintf(' Now and in the sequel press any key to proceed\n\n');
pause

%simplebvp
%pendulum
homoclinic

simplebvp2
gravitation


function simplebvp
clc;home;
figure
hold on;
fprintf('\n We start with a very simple first order problem : \n\n');
fprintf(' x''=x   ,   x(0)+x(2)=1 on a nonequidistant grid\n\n');
fprintf(' The solution can be easily found by hand.\n\n');
text(0.2,0.8,'1/(1+e^2) * e^t');

t0=linspace(0,2,50);
plot(t0,1/(1+exp(2))*exp(t0));
pause



fprintf(' Now we use BVP to solve the problem numerically \n on a nonequidistant grid\n\n');
fprintf(' >>t=2*([0:0.05:1].^1.5)'';     -> choose time grid\n >>x0=t;     -> choose initial approximation');
fprintf('\n >>z=bvp(''bvpdemo1'',t,x0)      -> solve the problem defined in ''bvpdemo1.m'';\n >>plot(t,z,''ro'');\n');


t=2*([0:0.05:1].^1.5)';      % time grid
x0=t;                   % start approximation

z=bvp('bvpdemo1',t,x0);
plot(t,z,'ro');

pause
fprintf('\n\n');   


%*****************************************************************************************

function pendulum
clc;
home;
clf;
fprintf('\n Now for a system : the pendulum.\n x''=y ,  y''=-sin(x) \n\n');
fprintf(' It is well known that the frequency of the pendumlum \n depends on the amplitude.');
fprintf(' We want to face the rather \n formidable task of determining the amplitude for a given\n');
fprintf(' frequency, i.e. we want to find a solution that is periodic\n with period T.\n\n')

fprintf(' So we have to choose x(0)=x(T) , y(0)=y(T) as boundary conditions.\n');
fprintf(' Note that this problem is not well-posed, because with phi(t),\n phi(t+a) is a solution');
fprintf(' too (for arbitrary a).\n\n');

[x,y]=meshgrid(-1:0.1:1);
contour(x,y,y.^2/2-cos(x),[-1:0.1:1]);
colormap summer;
set(gca,'xgrid','off','ygrid','off');

pause

fprintf(' We plot an initial approximation... \n\n');
t=linspace(0,2*pi+0.1,50)';      % time grid
x0=0.7*([cos(t) sin(t)]+(rand(length(t),2)-0.5));


fprintf('>>t=linspace(0,2*pi+0.1,50)'';\n\n>>x0=0.7*([cos(t) sin(t)]+(rand(length(t),2)-0.5));\n\n');
fprintf('>>bvpphas2(t,x0,''init'')      -> 2-dim. phaseplot of x0\n\n');

bvpphas2(t,x0,'init');

fprintf(' ... and start calculating\n\n');
fprintf('>> z=bvp(''bvpdemo2'',t,x0);\n\n');


z=bvp('bvpdemo2',t,x0);
bvpphas2(t,z,'');

pause
fprintf(' Now we try the same for a different period T\n');
fprintf(' and demonstrate the usage of an output function\n\n');
fprintf(' >>t=linspace(0,2*pi+0.2,50)''\n\n >>bvpopt=bvpset(''outputfcn'',''bvpphas2'');\n\n');
fprintf(' >>bvp(''bvpdemo2'',t,x0,bvpopt);\n\n');

hold on
t=linspace(0,2*pi+0.2,50)';      % time grid

bvpopt=bvpset('outputfcn','bvpphas2');

bvp('bvpdemo2',t,x0,bvpopt);

pause
fprintf('\n\n');   
%*****************************************************************************************


function homoclinic(E0,delta,N,sglstep)

clc;home;
clf
fprintf('\n The last first-order-system we deal with comes from transport theory.\n');
fprintf(' We are searching for a homoclinic orbit corresponding to a travelling wave.\n\n');

E0=2;
delta=0.01;
N=100;

[t,x,te]=ode45('bvpdemo3',[],[],'',E0,delta);   % calculate hom. Orbit 

t15=interp1(t,linspace(1,length(t),15)');   % sample 15 points from the time-grid disribution
tN=interp1(t15,linspace(1,15,N)','spline'); % keep the distribution, but use N points
y=interp1(t,x,tN);                          % interpolate the hom. orbit onto new timegrid
h=diff(tN);                                 % stepsize vector

n=y(:,1);                                   % just for interpretation reasons
w=[0 ; diff(n)./h];
E=y(:,2);
c=-E0/(E0+1)^2;   % =v(E0)



% plot first integral


[nsurf,Esurf]=meshgrid(linspace(min(n)-0.1,max(n)+0.1,20),...
                       linspace(min(E)-0.1,max(E)+0.1,20)); 
m=mesh(Esurf,nsurf,(-Esurf./(Esurf+1).^2 - c).*nsurf);
set(m,'facecolor','none');
set(gca,'xgrid','off','ygrid','off','zgrid','off');
view(24,42);
shading interp;
colormap jet;
rotate3d on;
axis tight;
axis manual;
hold on


p=plot3(x(:,2),x(:,1),[0 ; diff(x(:,1))./diff(t)],'r');          % plot hom. orbit
set(p,'linewidth',2);

fopt=foptions;   % options for fsolve
fopt(1)=1;       % statistics
fopt(2:3)=1e-4;  % tolerance
%fopt(9)=1;      % check jacobian

x0=[n w E] + 0.3*(rand(N,3)-0.5);             % start approximation

bvpopt=bvpset('fsolveopt',fopt,'outputfcn','bvpgunplot');
z=bvp('bvpdemo4',tN,x0,bvpopt,[E0 c]);  
pause
fprintf('\n\n');   
%*****************************************************************************************

function simplebvp2
clc;home;
clf
hold on
t=(2*linspace(0.05,1,20).^1.5)';
t=[-t([20:-1:1]) ; t];
x0=t;

fprintf('\n As an example for a second order equation we take\n\n');
fprintf(' x"=x , x(-2)=x(2)=cosh(2)\n\n');
fprintf(' Of course we know that cosh(t) is the solution, so we plot it\n\n');

t0=linspace(-2,2,50)';
plot(t0,cosh(t0));
axis tight
axis manual

pause

fprintf(' Now we solve the problem on a nonequidistant grid\n\n');
fprintf(' >>t=...;\n >>x0=...;\n >> [z,zprime]=bvp2(''bvpdemo5'',t,x0);\n');
fprintf(' >>plot(t,z,''ro'');\n');

[z,zprime]=bvp2('bvpdemo5',t,x0);
plot(t,z,'ro');
pause
fprintf('\n We also obtain a good approximation of the derivative sinh(t)\n\n\n');
clf
hold on
plot(t0,sinh(t0));
plot(t,zprime,'ro');
pause
fprintf('\n\n');   
%*****************************************************************************************

function gravitation
clc;home;
clf
hold on    
axis equal
colormap cool

phi=linspace(0,2*pi,10);        % plot mass
fill(0.1*cos(phi),0.1*sin(phi),'r');

[x,y]=meshgrid(-2.2:0.06:2.2); % plot potential
v=linspace(0,-3,20)';
contour(x,y,-1./sqrt(x.^2+y.^2),-v.^2);

axis tight
axis manual

fprintf(' A gravitational problem in 2 space coordinates:\n\n x"=-x/r^3\n');
fprintf(' y"=-y/r^3\n\n with r=sqrt(x^2+y^2)\n\n');
fprintf(' We choose\n\n x(0)=2 , y(0)=1 , y(T)=-2 , x''(T)=0\n\n as boundary conditions.\n\n')
pause

fprintf(' Again we use the output function BVPPHAS2 to display intermediate \n approximations.');
fprintf(' Setting the ''trace'' property to ''on'' enables you to single-step the iteration.\n\n');
fprintf(' >>t=...;  x0=... ;\n');
fprintf(' >>bvpopt=bvpset(''outputfcn'',''bvpphas2'',''trace'',''on''); -> gather bvpoptions \n\n');
fprintf(' >>z=bvp2(''bvpdemo6'',t,x0,bvpopt);\n\n');
fprintf(' The output function waits for you to start the iteration\n');

t=linspace(0,2,50)';            % start approximation
phi=linspace(pi/2-1,pi+1,length(t))';
x0=2*[cos(phi) sin(phi)];
x0=x0+(rand(length(t),2)-0.5);

fopt=foptions;
fopt(1)=1;
bvpopt=bvpset('outputfcn','bvpphas2','trace','on','fsolveopt',fopt);


z=bvp2('bvpdemo6',t,x0,bvpopt); % off we go

x=z(:,1);                       % plot arrows
y=z(:,2);
r=sqrt(x.^2+y.^2);
quiver(x,y,-x./r.^3,-y./r.^3);
pause(1);

x=interp1(x,linspace(1,length(x),400));
y=interp1(y,linspace(1,length(y),400));

for i=1:3

   head = line('color','r','marker','o','markerfacecolor','r','markersize',10,...
               'erase','xor','xdata',x(1),'ydata',y(1));
   body = line('color','r','linestyle','--','linewidth',3,'erase','none', ...
               'xdata',[],'ydata',[]);
   tail = line('color','y','linestyle','-','erase','none', ...
               'xdata',[],'ydata',[]);

   m = length(x);
   k = round(0.2*m);

   % Grow the body
   for i = 2:k+1
      j = i-1:i;
      set(head,'xdata',x(i),'ydata',y(i))
      set(body,'xdata',x(j),'ydata',y(j))
      drawnow
   end

   % Primary loop
   for i = k+2:m
      j = i-1:i;
      set(head,'xdata',x(i),'ydata',y(i))
      set(body,'xdata',x(j),'ydata',y(j))
      set(tail,'xdata',x(j-k),'ydata',y(j-k))
      drawnow
   end

   % Clean up the tail
   for i = m+1:m+k
      j = i-1:i;
      set(tail,'xdata',x(j-k),'ydata',y(j-k))
      drawnow
   end
   delete(head,body,tail);

end
fprintf('\n\n');   


%*****************************************************************************************









