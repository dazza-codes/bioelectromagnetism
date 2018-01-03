
% STARTUP	Define a few surfaces for testing, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['  startup.m ']);

clear
more on

%%% check screen color resolution
disp(['  running in ',num2str(get(0,'screendepth')), '-bit mode']);

%%% modify search path
disp(['  appending local search path...']);

path(path,[pwd,'/misc']);
path(path,[pwd,'/surffreq']);
cd surffreq

% suppress "Flop counts are no longer available" warning
warning off MATLAB:flops:UnavailableFunction

% list current variables, long form
whos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
