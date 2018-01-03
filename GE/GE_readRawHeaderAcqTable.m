function acq_tab = GE_readRawHeaderAcqTable(fid,sliceno)
%
% acq_tab = GE_readRawHeaderAcqTable(fid,sliceno)
%
% Loads the acq_table located in a file with file id fid
% and returns it as a structure. 
%
% Souheil J. Inati
% Dartmouth College
% August 2000
% souheil.inati@dartmouth.edu


% define the array of structures and read in the data
for i=1:sliceno
  acq_tab(i) = struct( ...
  'pass_number',   fread(fid,1,'int16'),   ... %which pass this slice is in
  'slice_in_pass', fread(fid,1,'int16'),   ... % which slice in this pass 
  'gw_point1',     fread(fid,3,'float32'), ... % corner points of image
  'gw_point2',     fread(fid,3,'float32'), ...
  'gw_point3',     fread(fid,3,'float32')  ...
  );
end

% The C stuff from GE
%#define RDB_MAX_SLICES                 512
%typedef struct
%{
%        short   pass_number;           /* which pass this slice is in */
%        short   slice_in_pass;         /* which slice in this pass */
%        float   gw_point1[3];          /* corner points of image */
%        float   gw_point2[3];
%        float   gw_point3[3];
%}       RDB_SLICE_INFO_ENTRY;
%
%typedef RDB_SLICE_INFO_ENTRY RDB_DATA_ACQ_TAB[RDB_MAX_SLICES];
