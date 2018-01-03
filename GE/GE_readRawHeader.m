function [rdb_hdr,acq_tab,ex_hdr,se_hdr,im_hdr,sizeofpoolheader] = ...
    GE_readRawHeader(PFileName, byte_align)
%
% [rdb_hdr,acq_tab,ex_hdr,se_hdr,im_hdr,sizeofpoolheader] =
% GE_readRawHeader(PFileName, byte_align)
%
% reads the headers from a P file. 
% byte_align = 0 for 5X
% byte_align = 1 for LX
%
% Souheil J. Inati
% Dartmouth College
% April 2001
% souheil.inati@dartmouth.edu

% Define the header offsets and sizes
hdr_sz = zeros(10,1);   %header section sizes
hdr_sz(1) = 2048;       %rdb_hdr_rec
hdr_sz(2) = 4096;       %rdb_hdr_per_pass
hdr_sz(3) = 4096;       %rdb_hdr_unlock_raw
hdr_sz(4) = 20480;      %rdb_hdr_data_acq_tab
hdr_sz(5) = 2052;       %rdb_hdr_nex_tab
hdr_sz(6) = 2052;       %rdb_hdr_nex_abort_tab
hdr_sz(7) = 2048;	%rdb_hdr_tool
if (byte_align == 0)
  hdr_sz(8) = 1024;       %rdb_hdr_exam
  hdr_sz(9) = 1020;       %rdb_hdr_series
  hdr_sz(10) = 1022;      %rdb_hdr_image
elseif (byte_align == 1)
  hdr_sz(8) = 1040;       %rdb_hdr_exam
  hdr_sz(9) = 1028;       %rdb_hdr_series
  hdr_sz(10) = 1044;      %rdb_hdr_image
end
rdb_hdr_offset = 0;
acq_tab_offset = sum(hdr_sz(1:3));
ex_hdr_offset = sum(hdr_sz(1:7));
se_hdr_offset = sum(hdr_sz(1:8));
im_hdr_offset = sum(hdr_sz(1:9));
sizeofpoolheader = sum(hdr_sz);

%%%% Open the PFile %%%%%
fid=fopen(PFileName,'r','b');         %note: 'b' = Big-Endian format
if fid == -1
  msg = sprintf('GE_readRawHeader: Could not open %s file',PFileName);
  error(msg);
end

%%%% Read the file Header info %%%%
% Load the rdb_header_rec
fseek(fid,rdb_hdr_offset,-1);
rdb_hdr = GE_readRawHeaderRdbRec(fid);

% Load the acquisition table
nslice = rdb_hdr.rdb_hdr_nslices;
fseek(fid,acq_tab_offset,-1);
acq_tab = GE_readRawHeaderAcqTable(fid,nslice);

% Load the exam header
fseek(fid,ex_hdr_offset,-1);
ex_hdr = GE_readHeaderExam(fid, byte_align);

% Load the series header
fseek(fid,se_hdr_offset,-1);
se_hdr = GE_readHeaderSeries(fid, byte_align);

% Load the image header
fseek(fid,im_hdr_offset,-1);
im_hdr = GE_readHeaderImage(fid, byte_align);

fclose(fid);
return;







