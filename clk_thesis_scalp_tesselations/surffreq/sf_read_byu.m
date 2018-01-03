function [extr] = sf_read_byu(filename)
% SF_READ_BYU	Read data from MOVIE.BYU file
%   	    	extr = SF_READ_BYU(filename)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES

%%% THINGS TO DO
% - read connectivity if it would be useful

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments 
if nargin < 1
    help sf_read_byu; return
end

fid = fopen(filename, 'r');
counts = fscanf(fid,'%d',4);

fscanf(fid,'%d',2*counts(1));
extr = fscanf(fid, '%g', [3 counts(2)])';

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

