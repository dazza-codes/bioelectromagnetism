function opt=bvpset(varargin)

%   BVPSET set options for BVP and BVP2
%
%   MYOPTIONS=BVPSET([property,value,property,value,...]);
%
%   Properties are:
%   'OutputFcn' - Name of installable output function  [ string ]
%      At the start of approximation, the solver calls OUTPUTFCN(T,Y0,'init')
%      to initialize the output function.  After each approximation step with
%      approximation vector Y the solver calls STOP = OUTPUTFCN(T,Y,'').
%      If STOP==1 iteration is stopped and the solver calls OUTPUTFCN(T,Y,'done').
%
%   'OutputSel' - Output selection indices  [ vector of integers ]
%      This vector of indices specifies which components of the solution vector
%      are passed to the OutputFcn. OutputSel defaults to all components.
%
%   'FsolveOpt' -  Parameter struct used by the optimization routines [struct]
%      Type 'help foptions' for details.
%   
%   'Trace' - Flag that allows stepwise proceeding  ['on'/'off'] 
%      Set 'Trace' to 'on' if you want to pause after each iteration step.
%      If no output function is defined, 'Trace' is ignored.
%
%   Example:
%    my_fsolve_options=foptions;
%    my_fsolve_options(1)=1;   %   display some statistics
%    my_bvp_options=bvpset('outputfcn','bvpphas2','outputsel',[1 3],...
%                          'fsolveopt',my_fsolve_options);
%    z=bvp('mybvpfile',t,x0,my_bvp_options,parametervector);
%
%   See also BVP  BVP2  BVPPLOT  BVPPHAS2  BVPPHAS3  FOPTIONS

%   Copyright (c) 1999 Guenter Kneisl
%                      University of Technology, Vienna
%                      e9425595@fbma.tuwien.ac.at


opt.outputfcn='';
opt.outputsel=[];
opt.fsolve=foptions;
opt.trace='off';

for i=1:2:length(varargin)
   switch lower(varargin{i})    
   case 'outputfcn'  % output function
       opt.outputfcn=varargin{i+1};
   case 'outputsel'  % output
       opt.outputsel=varargin{i+1};
   case 'fsolveopt' % stats
      opt.fsolve=varargin{i+1};
   case 'trace'  % trace iterations
      opt.trace=varargin{i+1};
   otherwise
      error(['Invalid property : ''' varargin{i} '''']);
   end
end
