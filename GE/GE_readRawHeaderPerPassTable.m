function pp_tab = GE_readRawHeaderPerPassTable(fid,sliceno)
%
%  pp_tab = GE_readRawHeaderPerPassTable(fid.sliceno)
%
% Loads the per pass table  located in a file with filed id fid
% and returns it as a structure. 
%
% Souheil J. Inati
% Dartmouth College
% August 2000
% souheil.inati@dartmouth.edu

% define the array of structures and read in the data
for i=1:sliceno
  pp_tab(i) = struct( ...
  'bam_modifier',   fread(fid,1,'int32'),   ...
  'bam_address',    fread(fid,1,'int32'));
end

% The C stuff from GE
%#define RDB_MAX_SLICES                 512
%typedef struct
%{
%        long    bam_modifier;
%        long    bam_address;
%}       VME_ADDRESS;
%
%typedef struct
%{
%        VME_ADDRESS dab_bam[RDB_MAX_DABS];
%}       RDB_PASS_INFO_ENTRY;
%
%typedef RDB_PASS_INFO_ENTRY RDB_PER_PASS_TAB[RDB_MAX_PASSES];
