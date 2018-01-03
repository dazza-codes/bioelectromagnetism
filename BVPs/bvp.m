function out=bvp(f,t,x0,bvpopt,parvec)
%   BVP Solve first order boundary value problem
%
%   X=BVP(F,T,X0[,OPTIONS,PARVEC]) solves the boundary
%   value problem
%
%      x'=f(t,x) , t \in [a,b] , R(x(a),x(b))=0
%
%   using (nonequidistant) trapezoidal rule discretisation :
%
%          x1-x0=h1/2*(f(t0,x0)+f(t1,x1))
%
%   F is the name of an m-file that defines the function f,
%   the boundary conditions and some Jacobians (for details see
%   the BVP_File_Pattern at the end of this text.
%   
%   T is a strictly monotonic sequence that defines the time grid.
%   T should be a column vector.
%
%   X0 is a starting guess and should be a (length(T) x dim)-matrix,
%   where dim is the dimension of the system to be solved.
%
%   [OPTIONS is a struct containing some approximation options like
%   tolerances and output functions. See BVPSET for details.]
%
%   [PARVEC is a parameter vector (or anything else you would like
%   to be) passed on to F. See BVP_File_Pattern below.]
%
%   ***************************************************************
%
%   The output X should be a solution and has the same form as X0.
%
%   ***************************************************************
%
%   function out=BVP_File_Pattern(t,x,XN,flag[,parvec])
%
%   sigma=parvec(1); mu=parvec(2); ...
%
%   switch flag
%   case 'f'  
%     out=<f(t,x)>;
%   case 'df/dx' 
%     out=<df/dx(t,x)>;
%   case 'R'  % boundary condition
%     x1=x;   % x1 is passed as second argument;
%     out=<R(x1,xN)>;
%   case 'dR/dx1'   
%     x1=x;
%     out=<dR/dx1(x1,xN)>;
%   case 'dR/dxN'
%     x1=x;
%     out=<dR/dxN(x1,xN)>;
%   otherwise
%     error('unknown flag');
%   end
%
%   See also BVP2  BVPSET  BVPPLOT  BVPPHAS2  BVPPHAS3

%   Copyright (c) 1999 Guenter Kneisl
%                      University of Technology, Vienna
%                      e9425595@fbma.tuwien.ac.at

error(nargchk(3,5,nargin));

if size(t,2) ~= 1
   t=t';
end
N=length(t);
dim=size(x0,2);
h2=diff(t)/2;

x=x0';
x=x(:);                    % lineare Indizierung


switch nargin               % option handling
case 3
   bvpopt=bvpset;
   parvec=[];
case 4
   parvec=[];
end

if nargout==0 & isempty(bvpopt.outputfcn)
   bvpopt.outputfcn='bvpplot';   % default output function
end
if isempty(bvpopt.outputsel)
   bvpopt.outputsel=[1:dim];
end


   



if isempty(bvpopt.outputfcn)
   x=sparsefsolve('bvp_disc',x,bvpopt.fsolve,'bvp_disc_jac',f,t,dim,N,h2,bvpopt,parvec);
else
   outp=x0(:,bvpopt.outputsel);
   stop=eval([bvpopt.outputfcn '(t,outp,''init'')']);
   
   x=sparsefsolve('bvp_disc',x,bvpopt.fsolve,'bvp_disc_jac',f,t,dim,N,h2,bvpopt,parvec);
   
   z=zeros(dim,N);                  
   z(:)=x;
   outp=z(bvpopt.outputsel,:)';
   stop=eval([bvpopt.outputfcn '(t,outp,''done'')']);
end

if nargout==1
   z=zeros(dim,N);                  
   z(:)=x;
   out=z';
end




function [x,OPTIONS] = sparsefsolve(FUN,x,OPTIONS,GRADFUN,varargin)
% SPARSEFSOLVE Solves nonlinear equations by a least squares method.
% Modification of FSOLVE suitable for sparse matrices.


% Handle undefined arguments
if nargin < 2, error('fsolve requires two input arguments');end
if nargin<4
    GRADFUN=[];
    if nargin<3
        OPTIONS=[];
    end
end

% Check for old syntax
if length(GRADFUN)  & ~isstr(GRADFUN) 
   disp('The user-supplied gradient (Jacobian) must be a string (function name).');
    error('The syntax to fsolve has been changed  - refer to the Optimization Toolbox guide');
end

if length(OPTIONS)<5; 
    OPTIONS(5)=0; 
end
% Switch methods making Gauss Newton the default method.
if OPTIONS(5)==0; OPTIONS(5)=1; else OPTIONS(5)=0; end

% Convert to inline function as needed.
if ~isempty(FUN)
  [funfcn, msg] = fcnchk(FUN,length(varargin));
  if ~isempty(msg)
    error(msg);
  end
else
  error('FUN must be a function name or valid expression.')
end

if ~isempty(GRADFUN)
  [gradfcn, msg] = fcnchk(GRADFUN,length(varargin));
  if ~isempty(msg)
    error(msg);
  end
else
  gradfcn = [];
end

[x,OPTIONS] = sparsenlsq(funfcn,x,OPTIONS,gradfcn,varargin{:});

if OPTIONS(8)>10*OPTIONS(3) & OPTIONS(1)>0
    disp('Optimizer is stuck at a minimum that is not a root')
    disp('Try again with a new starting guess')
end

% end fsolve

function [x,OPTIONS,CostFunction,JACOB] = sparsenlsq(FUN,x,OPTIONS,GRADFUN,varargin)
% SPARSENLSQ Solves non-linear least squares problems.
% Modification of NLSQ suitable for sparse matrices.


XOUT = x(:);
[nvars] = length(XOUT);
how = [];

% Global parameters for outside control of leastsq
% OPT_STOP is used for prematurely stopping the optimization
% OPT_STEP is set to 1 during major (non-gradient finding) iterations
%          set to 0 during gradient finding and 2 during line search
%          this can be useful for plotting etc.
global OPT_STOP OPT_STEP;
OPT_STEP = 1;
OPT_STOP = 0;


CostFunction = feval(FUN,x,varargin{:});
CostFunction = CostFunction(:);
OPT_STEP = 0;  % No longer a major step

nfun=length(CostFunction);
GRAD=zeros(length(XOUT),nfun);
OLDX=XOUT;
MATX=zeros(3,1);
MATL=[CostFunction'*CostFunction;0;0];
OLDF=CostFunction'*CostFunction;
FIRSTF=CostFunction'*CostFunction;
[OLDX,OLDF,OPTIONS]=lsint(XOUT,CostFunction,OPTIONS);
PCNT = 0;
EstSum=0.5;
GradFactor=1; 
CHG = 1e-7*abs(XOUT)+1e-7*ones(nvars,1);

OPTIONS(10)=1;
status=-1;


while status~=1  & OPT_STOP == 0

% Work Out Gradients
    if isempty(GRADFUN) | OPTIONS(9)
        OLDF=CostFunction;
        CHG = sign(CHG+eps).*min(max(abs(CHG),OPTIONS(16)),OPTIONS(17));
        for gcnt=1:nvars
            temp = XOUT(gcnt);
            XOUT(gcnt) = temp +CHG(gcnt);
            x(:) = XOUT;
            CostFunction(:) = feval(FUN,x,varargin{:});
            OPT_STEP = 0; % We're in gradient finding mode
            GRAD(gcnt,:)=(CostFunction-OLDF)'/(CHG(gcnt));
            XOUT(gcnt) = temp;
        end
        CostFunction = OLDF;
        OPTIONS(10)=OPTIONS(10)+nvars;
% Gradient check
        if OPTIONS(9) == 1 & ~isempty(GRADFUN)
            GRADFD = GRAD;
            x(:)=XOUT; 
            GRAD = feval(GRADFUN,x,varargin{:});
            if isa(GRADFUN,'inline')
              graderr(GRADFD, GRAD, formula(GRADFUN));
            else
              graderr(GRADFD, GRAD,  GRADFUN);
            end
 
            OPTIONS(9) = 0;
        end
    else
        x(:) = XOUT;
        OPTIONS(11)=OPTIONS(11)+1;
        GRAD = feval(GRADFUN,x,varargin{:});
    end
    % Try to set difference to 1e-8 for next iteration
    if nfun==1
        if isequal(GRAD,0)
          CHG = Inf;
        else
          CHG = nfun*1e-8./GRAD;
        end
    else
        sumabsGRAD = sum(abs(GRAD)')';
        ii = (sumabsGRAD == 0);
        CHG(ii) = Inf;
        CHG(~ii) = nfun*1e-8./sumabsGRAD(~ii);
    end

    OPT_STEP = 2; % Entering line search    

    GradF = 2*GRAD*CostFunction;
    NewF = CostFunction'*CostFunction;
%---------------Initialization of Search Direction------------------
    if status==-1
        if condest(GRAD)>1e8
            SD=-(GRAD*GRAD'+(normest(GRAD)+1)*(eye(nvars,nvars)))\(GRAD*CostFunction);  % !!!!!!!!!!!!!!
            if OPTIONS(5)==0, GradFactor=normest(GRAD)+1; end  %!!!!!!!!!!!!!!!!!!!!!!!!!!!! NORM -> NORMEST
            how='condest';
        else
%       SD=GRAD'\(GRAD'*X-F)-X;
            SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*CostFunction);
        end
        FIRSTF=NewF;
        OLDG=GRAD;
        GDOLD=GradF'*SD;
        % OPTIONS(18) controls the initial starting step-size.
        % If OPTIONS(18) has been set externally then it will
        % be non-zero, otherwise set to 1.
        if OPTIONS(18) == 0, OPTIONS(18)=1; end
        if OPTIONS(1)>0
            disp([sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),NewF,OPTIONS(18)),sprintf('%12.3g  ',GDOLD)]);
        end
        XOUT=XOUT+OPTIONS(18)*SD;
        if OPTIONS(5)==0
            newf=GRAD'*SD+CostFunction;
            GradFactor=newf'*newf;
            SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*CostFunction); 
        end
        newf=GRAD'*SD+CostFunction;
        XOUT=XOUT+OPTIONS(18)*SD;
        EstSum=newf'*newf;
        status=0;
        if OPTIONS(7)==0; PCNT=1; end
        
    else
%-------------Direction Update------------------
        gdnew=GradF'*SD;
        if OPTIONS(1)>0, 
            num=[sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),NewF,OPTIONS(18)),sprintf('%12.3g  ',gdnew)];
        end
        if gdnew>0 & NewF>FIRSTF

% Case 1: New function is bigger than last and gradient w.r.t. SD -ve
% ... interpolate. 
            how='inter';
            [stepsize]=cubici1(NewF,FIRSTF,gdnew,GDOLD,OPTIONS(18));
            OPTIONS(18)=0.9*stepsize;
        elseif NewF<FIRSTF

%  New function less than old fun. and OK for updating 
%         .... update and calculate new direction. 
            [newstep,fbest] =cubici3(NewF,FIRSTF,gdnew,GDOLD,OPTIONS(18));
            if fbest>NewF,fbest=0.9*NewF; end 
            if gdnew<0
                how='incstep';
                if newstep<OPTIONS(18),  newstep=(2*OPTIONS(18)+1e-4); how=[how,'IF']; end
                OPTIONS(18)=abs(newstep);
            else 
                if OPTIONS(18)>0.9
                    how='int_step';
                    OPTIONS(18)=min([1,abs(newstep)]);
                end
            end
% SET DIRECTION.
% Gauss-Newton Method    
            temp=1;
            if OPTIONS(5)==1 
                if OPTIONS(18)>1e-8 & condest(GRAD)<1e8
                    SD=GRAD'\(GRAD'*XOUT-CostFunction)-XOUT;
                    if SD'*GradF>eps,how='ERROR- GN not descent direction',  end
                    temp=0;
                else
                    if OPTIONS(1) > 0
                        disp('Conditioning of Gradient Poor - Switching To LM method')
                    end
                    how='CHG2LM';
                    OPTIONS(5)=0;
                    OPTIONS(18)=abs(OPTIONS(18));               
                end
            end
            
            if (temp)      
% Levenberg_marquardt Method N.B. EstSum is the estimated sum of squares.
%                                 GradFactor is the value of lambda.
% Estimated Residual:
                if EstSum>fbest
                    GradFactor=GradFactor/(1+OPTIONS(18));
                else
                    GradFactor=GradFactor+(fbest-EstSum)/(OPTIONS(18)+eps);
                end
                SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*CostFunction); 
                OPTIONS(18)=1; 
                estf=GRAD'*SD+CostFunction;
                EstSum=estf'*estf;
                if OPTIONS(1)>0, num=[num,sprintf('%12.6g ',GradFactor)]; end
            end
            gdnew=GradF'*SD;

            OLDX=XOUT;
% Save Variables
            FIRSTF=NewF;
            OLDG=GradF;
            GDOLD=gdnew;    

            % If quadratic interpolation set PCNT
            if OPTIONS(7)==0, PCNT=1; MATX=zeros(3,1);  MATL(1)=NewF; end
        else 
% Halve Step-length
            how='Red_Step';
            if NewF==FIRSTF,
                if OPTIONS(1)>0,
                    disp('No improvement in search direction: Terminating')
                end
                status=1;
            else
                OPTIONS(18)=OPTIONS(18)/8;
                if OPTIONS(18)<1e-8
                    OPTIONS(18)=-OPTIONS(18);
                end
            end
        end
        XOUT=OLDX+OPTIONS(18)*SD;
         if isinf(OPTIONS(1))
            disp([num,how])
         elseif OPTIONS(1)>0
            disp(num)
         end

    end %----------End of Direction Update-------------------
    if OPTIONS(7)==0, PCNT=1; MATX=zeros(3,1);  MATL(1)=NewF; end
% Check Termination 
    if max(abs(SD))< OPTIONS(2) & (GradF'*SD) < OPTIONS(3) & max(abs(GradF)) < 10*(OPTIONS(3)+OPTIONS(2))
        if OPTIONS(1) > 0
            disp('Optimization Terminated Successfully')  
        end
        status=1; 
    elseif OPTIONS(10)>OPTIONS(14)
        disp('maximum number of iterations has been exceeded');
        if OPTIONS(1)>0
            disp('Increase OPTIONS(14)')
        end
        status=1;
    else

% Line search using mixed polynomial interpolation and extrapolation.
        if PCNT~=0
            while PCNT > 0 & OPT_STOP == 0
                x(:) = XOUT; 
                CostFunction(:) = feval(FUN,x,varargin{:});
                OPTIONS(10)=OPTIONS(10)+1;
                NewF = CostFunction'*CostFunction;
            % <= used in case when no improvement found.
                if NewF <= OLDF'*OLDF, OX = XOUT; OLDF=CostFunction; end
                [PCNT,MATL,MATX,steplen,NewF,how]=searchq(PCNT,NewF,OLDX,MATL,MATX,SD,GDOLD,OPTIONS(18),how);
                OPTIONS(18)=steplen;
                XOUT=OLDX+steplen*SD;
                if NewF==FIRSTF,  PCNT=0; end
            end
            XOUT = OX;
            CostFunction=OLDF;
        else
            x(:)=XOUT; 
            CostFunction(:) = feval(FUN,x,varargin{:});
            OPTIONS(10)=OPTIONS(10)+1;
        end
    end
    OPT_STEP = 1; % Call the next iteration a major step for plotting etc.
    
    if ~isempty(varargin{6}.outputfcn)  %varargin{6}=bvpopt
      z=zeros(varargin{3},varargin{4});                   % =zeros(dim,N)
      z(:)=x;
      outp=z(varargin{6}.outputsel,:)';
      OPT_STOP=eval([varargin{6}.outputfcn '(varargin{2},outp,'''')']);     % Transponieren nicht vergessen
      if varargin{6}.trace(1:2)=='on'
         pause
      end
    end
end
OPTIONS(8) = NewF;
XOUT=OLDX;
x(:)=XOUT;
JACOB = GRAD.';
%--end of leastsq--
























