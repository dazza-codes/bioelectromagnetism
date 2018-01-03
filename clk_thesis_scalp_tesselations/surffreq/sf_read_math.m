function [extr] = sf_read_math(filename1, filename2)
% SF_READ_MATH	Read data file (from Mathematica)
%   	    	extr = SF_READ_MATH(file1 [, file2])

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES

%%% THINGS TO DO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments 
if nargin < 1
    help sf_read_math; return
end

if nargin > 0
    fid = fopen(filename1);
    extr =            fscanf(fid,'%g %g %g',[3 inf])';
    fclose(fid);
end
if nargin > 1
    fid = fopen(filename2);
    extr = extr + i * fscanf(fid,'%g %g %g',[3 inf])';
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

