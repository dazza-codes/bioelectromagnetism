function [ ge ] = ge_hdr_read(imageFileName)

% ge_hdr_read - Reads the header info from a GE LX2 or 5.X file
% 
% [ge,im_offset] = ge_hdr_read('imageFileName')
%
% ge is a struct with fields:
% 
% ge.hdr.suite  - suite header
% ge.hdr.exam   - exam header
% ge.hdr.series - series header
% ge.hdr.image  - image header
% ge.hdr.pixel  - pixel header
% 
% ge.img.offset  - the byte offset to the image data
% 
% This function assumes the image file is big endian.  It can
% automatically determine LX2 format (SGI) and 5.X format (Sun);
% other formats fail.
% 
% See also ge_vol_read
% 


% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Souheil J. Inati  <souheil.inati@nyu.edu> at 03/2003
% Dartmouth College, May 2000
% 
% Darren.Weber_at_radiology.ucsf.edu, March 2003
% - Substantially redesigned file handling and function
%   call structures for integration with mri_toolbox at 
%   http://eeg.sf.net
% - Requested permission to distribute code under GPL licence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DLW, found this note in an FSL list message:
%This is further confounded by the GE image format I.001 or I.MR
%formats where the data following the header is inverted (top of the
%image for bottom).



fprintf('\nGE_HDR_READ\n'); tic

fprintf('...reading %s\n',imageFileName);

% Open the Image File
[fid,message] = fopen(imageFileName,'r','b'); % 'b' = Big-Endian format
if fid == -1, error(message); end

% Check for the magic number at the first word
magic = fread(fid,1,'int32');
if magic == 1229801286,
    % This is a 5.X format image
    fprintf('...identified Signa 5.X format\n');
    byte_align = 0;
    ge.hdr.pixel.offset = 0;
    ge.hdr.suite.size   = 114;
    ge.hdr.exam.size    = 1024;
    ge.hdr.series.size  = 1020;
    ge.hdr.image.size   = 1022;
else
    % Magic no. not at the 1st word, try at the LX2 position
    fseek(fid,3228,-1);
    magic = fread(fid,1,'int32');
    if magic == 1229801286,
        % This is an LX2 format image
        fprintf('...identified LX2 format\n');
        byte_align = 1;
        ge.hdr.pixel.offset = 3228;
        ge.hdr.suite.size   = 116;
        ge.hdr.exam.size    = 1040;
        ge.hdr.series.size  = 1028;
        ge.hdr.image.size   = 1044;
    else
        msg = 'This is not a 5.X or LX2 format image.  No Magic Number.';
        error(msg);
    end
end

% Load the pixel header
fseek(fid,ge.hdr.pixel.offset,-1);
ge = ge_readHeaderPixel(fid, byte_align, ge); % see subfunc below


% Check for epirecon images
if (ge.hdr.pixel.img_l_dbHdr == 0 & byte_align == 0),
    error('This is an epirecon image. No header!');
end

% Compute the offsets
ge.hdr.suite.offset  = ge.hdr.pixel.img_p_dbHdr;
ge.hdr.exam.offset   = ge.hdr.suite.offset  + ge.hdr.suite.size;
ge.hdr.series.offset = ge.hdr.exam.offset   + ge.hdr.exam.size;
ge.hdr.image.offset  = ge.hdr.series.offset + ge.hdr.series.size;
ge.img.offset        = ge.hdr.pixel.offset  + ge.hdr.pixel.img_hdr_length;



% Load the suite header
fseek(fid,ge.hdr.suite.offset,-1);
ge.hdr.suite = ge_readHeaderSuite(fid, byte_align); % see subfunc below

% Load the exam header
fseek(fid,ge.hdr.exam.offset,-1);
ge.hdr.exam = ge_readHeaderExam(fid, byte_align); % see subfunc below

% Load the series header
fseek(fid,ge.hdr.series.offset,-1);
ge.hdr.series = ge_readHeaderSeries(fid, byte_align); % see subfunc below

% Load the image header
fseek(fid,ge.hdr.image.offset,-1);
ge.hdr.image = ge_readHeaderImage(fid, byte_align); % see subfunc below

% Close the file
fclose(fid);

t=toc; fprintf('...done (%5.2f sec).\n',t);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function su_hdr = ge_readHeaderSuite(fid, byte_align)

% Loads the suite header from an image file (with fid)
% and returns it as a structure.
% 
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)

% define the structure and read
su_hdr = struct('su_id', fread(fid,4,'char'));                   % Suite ID
su_hdr = setfield(su_hdr, 'su_uniq', fread(fid,1,'int16'));      % The Make-Unique Flag
su_hdr = setfield(su_hdr, 'su_diskid', fread(fid,1,'char'));     % Disk ID
su_hdr = setfield(su_hdr, 'prodid', fread(fid,13,'char'));       % Product ID
su_hdr = setfield(su_hdr, 'su_verscre', fread(fid,2,'char'));    % Genesis Version
su_hdr = setfield(su_hdr, 'su_verscur', fread(fid,2,'char'));    % Genesis Version
su_hdr = setfield(su_hdr, 'su_checksum', fread(fid,1,'uint32')); % Suite Record Checksum
su_hdr = setfield(su_hdr, 'su_padding', fread(fid,85,'char'));   % Spare Space
fseek(fid,1,0); % 16-bit alignment
if byte_align, fseek(fid,2,0); end % 32-bit alignment

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ex_hdr = ge_readHeaderExam(fid, byte_align)

% Loads the exam header from the image file
% and returns it as a structure.
% 
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)

% define the structure and read ifseek(fid,1,0);  % byte align the data
% to overcome the byte alignment problems
% break up the assignment into pieces using the setfield function
ex_hdr = struct('ex_suid', fread(fid,4,'uchar')); %Suite ID for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_uniq', fread(fid,1,'int16'));        %The Make-Unique Flag%
ex_hdr = setfield(ex_hdr, 'ex_diskid', fread(fid,1,'uchar'));      %Disk ID for this Exam%
fseek(fid,1,0); % 16-bit alignment
ex_hdr = setfield(ex_hdr, 'ex_no', fread(fid,1,'uint16'));         %Exam Number%
ex_hdr = setfield(ex_hdr, 'hospname', fread(fid,33,'uchar'));      %Hospital Name%
fseek(fid,1,0); % 16-bit alignment
ex_hdr = setfield(ex_hdr, 'detect', fread(fid,1,'int16'));         %Detector Type%
if byte_align; fseek(fid,2,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'numcells', fread(fid,1,'int32'));       %Number of cells in det%
ex_hdr = setfield(ex_hdr, 'zerocell', fread(fid,1,'float32'));     %Cell number at theta%
ex_hdr = setfield(ex_hdr, 'cellspace', fread(fid,1,'float32'));    %Cell spacing%
ex_hdr = setfield(ex_hdr, 'srctodet', fread(fid,1,'float32'));     %Distance from source to detector%
ex_hdr = setfield(ex_hdr, 'srctoiso', fread(fid,1,'float32'));     %Distance from source to iso%
ex_hdr = setfield(ex_hdr, 'tubetyp', fread(fid,1,'int16'));        %Tube type%
ex_hdr = setfield(ex_hdr, 'dastyp', fread(fid,1,'int16'));         %DAS type%
ex_hdr = setfield(ex_hdr, 'num_dcnk', fread(fid,1,'int16'));       %Number of Decon Kernals%
ex_hdr = setfield(ex_hdr, 'dcn_len', fread(fid,1,'int16'));        %Number of elements in a Decon Kernal%
ex_hdr = setfield(ex_hdr, 'dcn_density', fread(fid,1,'int16'));    %Decon Kernal density%
ex_hdr = setfield(ex_hdr, 'dcn_stepsize', fread(fid,1,'int16'));   %Decon Kernal stepsize%
ex_hdr = setfield(ex_hdr, 'dcn_shiftcnt', fread(fid,1,'int16'));   %Decon Kernal Shift Count%
if byte_align; fseek(fid,2,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'magstrength', fread(fid,1,'int32'));    %Magnet strength (in gauss)%
ex_hdr = setfield(ex_hdr, 'patid', fread(fid,13,'uchar'));         %Patient ID for this Exam%
ex_hdr = setfield(ex_hdr, 'patname', fread(fid,25,'uchar'));       %Patientsda Name%
ex_hdr = setfield(ex_hdr, 'patage', fread(fid,1,'int16'));         %Patient Age (years, months or days)%
ex_hdr = setfield(ex_hdr, 'patian', fread(fid,1,'int16'));         %Patient Age Notation%
ex_hdr = setfield(ex_hdr, 'patsex', fread(fid,1,'int16'));         %Patient Sex%
ex_hdr = setfield(ex_hdr, 'patweight', fread(fid,1,'int32'));      %Patient Weight%
ex_hdr = setfield(ex_hdr, 'trauma', fread(fid,1,'int16'));         %Trauma Flag%
ex_hdr = setfield(ex_hdr, 'hist', fread(fid,61,'uchar'));          %Patient History%
ex_hdr = setfield(ex_hdr, 'reqnum', fread(fid,13,'uchar'));        %Requisition Number%
ex_hdr = setfield(ex_hdr, 'ex_datetime', fread(fid,1,'int32'));    %Exam date/time stamp%
ex_hdr = setfield(ex_hdr, 'refphy', fread(fid,33,'uchar'));        %Referring Physician%
ex_hdr = setfield(ex_hdr, 'diagrad', fread(fid,33,'uchar'));       %Diagnostician/Radiologist%
ex_hdr = setfield(ex_hdr, 'op', fread(fid,4,'uchar'));             %Operator%
ex_hdr = setfield(ex_hdr, 'ex_desc', fread(fid,23,'uchar'));       %Exam Description%
ex_hdr = setfield(ex_hdr, 'ex_typ', fread(fid,3,'uchar'));         %Exam Type%
ex_hdr = setfield(ex_hdr, 'ex_format', fread(fid,1,'int16'));      %Exam Format%
if byte_align; fseek(fid,6,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'firstaxtime', fread(fid,1,'float64'));  %Start time(secs) of first axial in exam%
ex_hdr = setfield(ex_hdr, 'ex_sysid', fread(fid,9,'uchar'));       %Creator Suite and Host%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
ex_hdr = setfield(ex_hdr, 'ex_lastmod', fread(fid,1,'int32'));     %Date/Time of Last Change%
ex_hdr = setfield(ex_hdr, 'protocolflag', fread(fid,1,'int16'));   %Non-Zero indicates Protocol Exam%
ex_hdr = setfield(ex_hdr, 'ex_alloc_key', fread(fid,13,'uchar'));  %Process that allocated this record%
fseek(fid,1,0); % 16-bit alignment
ex_hdr = setfield(ex_hdr, 'ex_delta_cnt', fread(fid,1,'int32'));   %Indicates number of updates to header%
ex_hdr = setfield(ex_hdr, 'ex_verscre', fread(fid,2,'uchar'));     %Genesis Version - Created%
ex_hdr = setfield(ex_hdr, 'ex_verscur', fread(fid,2,'uchar'));     %Genesis Version - Now%
ex_hdr = setfield(ex_hdr, 'ex_checksum', fread(fid,1,'uint32'));   %Exam Record Checksum%
ex_hdr = setfield(ex_hdr, 'ex_complete', fread(fid,1,'int32'));    %Exam Complete Flag%
ex_hdr = setfield(ex_hdr, 'ex_seriesct', fread(fid,1,'int32'));    %Last Series Number Used%
ex_hdr = setfield(ex_hdr, 'ex_numarch', fread(fid,1,'int32'));     %Number of Series Archived%
ex_hdr = setfield(ex_hdr, 'ex_numseries', fread(fid,1,'int32'));   %Number of Series Existing%
ex_hdr = setfield(ex_hdr, 'ex_series', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32')));  %Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_numunser', fread(fid,1,'int32'));    %Number of Unstored Series%
ex_hdr = setfield(ex_hdr, 'ex_unseries', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32')));    %Unstored Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_toarchcnt', fread(fid,1,'int32'));   %Number of Unarchived Series%
ex_hdr = setfield(ex_hdr, 'ex_toarchive', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32')));     %Unarchived Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_prospcnt', fread(fid,1,'int32'));    %Number of Prospective/Scout Series%
ex_hdr = setfield(ex_hdr, 'ex_prosp', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32'))); %Prospective/Scout Series Keys for this Exam%
ex_hdr = setfield(ex_hdr, 'ex_modelnum', fread(fid,1,'int32'));    %Last Model Number used%
ex_hdr = setfield(ex_hdr, 'ex_modelcnt', fread(fid,1,'int32'));    %Number of ThreeD Models%
ex_hdr = setfield(ex_hdr, 'ex_models', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32')));  %ThreeD Model Keys for Exam%
ex_hdr = setfield(ex_hdr, 'ex_stat', fread(fid,1,'int16'));        %Patient Status%
ex_hdr = setfield(ex_hdr, 'uniq_sys_id', fread(fid,16,'uchar'));   %Unique System ID%
ex_hdr = setfield(ex_hdr, 'service_id', fread(fid,16,'uchar'));    %Unique Service ID%
ex_hdr = setfield(ex_hdr, 'mobile_loc', fread(fid,4,'uchar'));     %Mobile Location Number%
ex_hdr = setfield(ex_hdr, 'study_uid', fread(fid,32,'uchar'));     %Study Entity Unique ID%
ex_hdr = setfield(ex_hdr, 'study_status', fread(fid,1,'int16'));   %indicates if study has complete info(DICOM/genesis)%
ex_hdr = setfield(ex_hdr, 'ex_padding', fread(fid,516,'uchar'));   %Spare Space%
if byte_align; fseek(fid,4,0); end % byte alignment

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function se_hdr = ge_readHeaderSeries(fid, byte_align)

% returns the series header as a structure.
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)

% define the structure and read in the data
% to overcome the byte alignment problems
% break up the assignment into pieces using the setfield function
se_hdr = struct('se_suid', fread(fid,4,'uchar'));                       %Suite ID for this Series%
se_hdr = setfield(se_hdr, 'se_uniq', fread(fid,1,'int16'));            %The Make-Unique Flag%
se_hdr = setfield(se_hdr, 'se_diskid', fread(fid,1,'uchar'));          %Disk ID for this Series%
fseek(fid,1,0); % 16-bit alignment
se_hdr = setfield(se_hdr, 'se_exno', fread(fid,1,'uint16'));            %Exam Number%
se_hdr = setfield(se_hdr, 'se_no ', fread(fid,1,'int16'));              %Series Number%
se_hdr = setfield(se_hdr, 'se_datetime', fread(fid,1,'int32'));        %Allocation Series Data/Time stamp%
se_hdr = setfield(se_hdr, 'se_actual_dt', fread(fid,1,'int32'));       %Actual Series Data/Time stamp%
se_hdr = setfield(se_hdr, 'se_desc', fread(fid,30,'uchar'));       %Series Description%
se_hdr = setfield(se_hdr, 'pr_sysid', fread(fid,9,'uchar'));       %Primary Receiver Suite and Host%
se_hdr = setfield(se_hdr, 'pansysid', fread(fid,9,'uchar'));       %Archiver Suite and Host%
se_hdr = setfield(se_hdr, 'se_typ ', fread(fid,1,'int16'));             %Series Type%
se_hdr = setfield(se_hdr, 'se_source ', fread(fid,1,'int16'));          %Series from which prescribed%
se_hdr = setfield(se_hdr, 'se_plane ', fread(fid,1,'int16'));           %Most-like Plane (for L/S)%
se_hdr = setfield(se_hdr, 'scan_type ', fread(fid,1,'int16'));          %Scout or Axial (for CT)%
se_hdr = setfield(se_hdr, 'position', fread(fid,1,'int32'));           %Patient Position%
se_hdr = setfield(se_hdr, 'entry', fread(fid,1,'int32'));              %Patient Entry%
se_hdr = setfield(se_hdr, 'anref', fread(fid,3,'uchar'));          %Anatomical reference%
fseek(fid,1,0); % 16-bit alignment
se_hdr = setfield(se_hdr, 'lmhor', fread(fid,1,'float32'));              %Horizontal Landmark%
se_hdr = setfield(se_hdr, 'prtcl', fread(fid,25,'uchar'));         %Scan Protocol Name%
fseek(fid,1,0); % 16-bit alignment
se_hdr = setfield(se_hdr, 'se_contrast ', fread(fid,1,'int16'));        %Non-zero if > 0 image used contrast(L/S)%
se_hdr = setfield(se_hdr, 'start_ras', fread(fid,1,'uchar'));          %RAS letter for first scan location (L/S)%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'start_loc', fread(fid,1,'float32'));          %First scan location (L/S)%
se_hdr = setfield(se_hdr, 'end_ras', fread(fid,1,'uchar'));            %RAS letter for last scan location (L/S)%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'end_loc', fread(fid,1,'float32'));            %Last scan location (L/S)%
se_hdr = setfield(se_hdr, 'se_pseq ', fread(fid,1,'int16'));            %Last Pulse Sequence Used (L/S)%
se_hdr = setfield(se_hdr, 'se_sortorder ', fread(fid,1,'int16'));       %Image Sort Order (L/S)%
se_hdr = setfield(se_hdr, 'se_lndmrkcnt', fread(fid,1,'int32'));       %Landmark Counter%
se_hdr = setfield(se_hdr, 'se_nacq ', fread(fid,1,'int16'));            %Number of Acquisitions%
se_hdr = setfield(se_hdr, 'xbasest ', fread(fid,1,'int16'));            %Starting number for baselines%
se_hdr = setfield(se_hdr, 'xbaseend', fread(fid,1,'int16'));            %Ending number for baselines%
se_hdr = setfield(se_hdr, 'xenhst', fread(fid,1,'int16'));             %Starting number for enhanced scans%
se_hdr = setfield(se_hdr, 'xenhend', fread(fid,1,'int16'));            %Ending number for enhanced scans%
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'se_lastmod', fread(fid,1,'int32'));         %Date/Time of Last Change%
se_hdr = setfield(se_hdr, 'se_alloc_key', fread(fid,13,'uchar'));  %Process that allocated this record%
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid,2,0); end % 32-bit alignment
se_hdr = setfield(se_hdr, 'se_delta_cnt', fread(fid,1,'int32'));       %Indicates number of updates to header%
se_hdr = setfield(se_hdr, 'se_verscre', fread(fid,2,'uchar'));     %Genesis Version - Created%
se_hdr = setfield(se_hdr, 'se_verscur', fread(fid,2,'uchar'));     %Genesis Version - Now%
se_hdr = setfield(se_hdr, 'se_pds_a', fread(fid,1,'float32'));           %PixelData size - as stored%
se_hdr = setfield(se_hdr, 'se_pds_c', fread(fid,1,'float32'));           %PixelData size - Compressed%
se_hdr = setfield(se_hdr, 'se_pds_u', fread(fid,1,'float32'));           %PixelData size - UnCompressed%
se_hdr = setfield(se_hdr, 'se_checksum', fread(fid,1,'uint32'));        %Series Record checksum%
se_hdr = setfield(se_hdr, 'se_complete', fread(fid,1,'int32'));        %Series Complete Flag%
se_hdr = setfield(se_hdr, 'se_numarch', fread(fid,1,'int32'));         %Number of Images Archived%
se_hdr = setfield(se_hdr, 'se_imagect', fread(fid,1,'int32'));         %Last Image Number Used%
se_hdr = setfield(se_hdr, 'se_numimages', fread(fid,1,'int32'));       %Number of Images Existing%
se_hdr = setfield(se_hdr, 'se_images', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32'))); %Image Keys for this Series%
se_hdr = setfield(se_hdr, 'se_numunimg', fread(fid,1,'int32'));        %Number of Unstored Images%
se_hdr = setfield(se_hdr, 'se_unimages', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32'))); %Unstored Image Keys for this Series%
se_hdr = setfield(se_hdr, 'se_toarchcnt', fread(fid,1,'int32'));       %Number of Unarchived Images%
se_hdr = setfield(se_hdr, 'se_toarchive', struct('length', fread(fid,1,'uint32'), ...
    'data', fread(fid,1,'uint32'))); %Unarchived Image Keys for this Series%
se_hdr = setfield(se_hdr, 'echo1_alpha', fread(fid,1,'float32'));        %Echo 1 Alpha Value%
se_hdr = setfield(se_hdr, 'echo1_beta', fread(fid,1,'float32'));         %Echo 1 Beta Value%
se_hdr = setfield(se_hdr, 'echo1_window', fread(fid,1,'uint16'));       %Echo 1 Window Value%
se_hdr = setfield(se_hdr, 'echo1_level', fread(fid,1,'int16'));        %Echo 1 Level Value%
se_hdr = setfield(se_hdr, 'echo2_alpha', fread(fid,1,'float32'));        %Echo 2 Alpha Value%
se_hdr = setfield(se_hdr, 'echo2_beta', fread(fid,1,'float32'));         %Echo 2 Beta Value%
se_hdr = setfield(se_hdr, 'echo2_window', fread(fid,1,'uint16'));       %Echo 2 Window Value%
se_hdr = setfield(se_hdr, 'echo2_level', fread(fid,1,'int16'));        %Echo 2 Level Value%
se_hdr = setfield(se_hdr, 'echo3_alpha', fread(fid,1,'float32'));        %Echo 3 Alpha Value%
se_hdr = setfield(se_hdr, 'echo3_beta', fread(fid,1,'float32'));         %Echo 3 Beta Value%
se_hdr = setfield(se_hdr, 'echo3_window', fread(fid,1,'uint16'));       %Echo 3 Window Value%
se_hdr = setfield(se_hdr, 'echo3_level', fread(fid,1,'int16'));        %Echo 3 Level Value%
se_hdr = setfield(se_hdr, 'echo4_alpha', fread(fid,1,'float32'));        %Echo 4 Alpha Value%
se_hdr = setfield(se_hdr, 'echo4_beta', fread(fid,1,'float32'));         %Echo 4 Beta Value%
se_hdr = setfield(se_hdr, 'echo4_window', fread(fid,1,'uint16'));       %Echo 4 Window Value%
se_hdr = setfield(se_hdr, 'echo4_level', fread(fid,1,'int16'));        %Echo 4 Level Value%
se_hdr = setfield(se_hdr, 'echo5_alpha', fread(fid,1,'float32'));        %Echo 5 Alpha Value%
se_hdr = setfield(se_hdr, 'echo5_beta', fread(fid,1,'float32'));         %Echo 5 Beta Value%
se_hdr = setfield(se_hdr, 'echo5_window', fread(fid,1,'uint16'));       %Echo 5 Window Value%
se_hdr = setfield(se_hdr, 'echo5_level', fread(fid,1,'int16'));        %Echo 5 Level Value%
se_hdr = setfield(se_hdr, 'echo6_alpha', fread(fid,1,'float32'));        %Echo 6 Alpha Value%
se_hdr = setfield(se_hdr, 'echo6_beta', fread(fid,1,'float32'));         %Echo 6 Beta Value%
se_hdr = setfield(se_hdr, 'echo6_window', fread(fid,1,'uint16'));       %Echo 6 Window Value%
se_hdr = setfield(se_hdr, 'echo6_level', fread(fid,1,'int16'));        %Echo 6 Level Value%
se_hdr = setfield(se_hdr, 'echo7_alpha', fread(fid,1,'float32'));        %Echo 7 Alpha Value%
se_hdr = setfield(se_hdr, 'echo7_beta', fread(fid,1,'float32'));         %Echo 7 Beta Value%
se_hdr = setfield(se_hdr, 'echo7_window', fread(fid,1,'uint16'));       %Echo 7 Window Value%
se_hdr = setfield(se_hdr, 'echo7_level', fread(fid,1,'int16'));        %Echo 7 Level Value%
se_hdr = setfield(se_hdr, 'echo8_alpha', fread(fid,1,'float32'));        %Echo 8 Alpha Value%
se_hdr = setfield(se_hdr, 'echo8_beta', fread(fid,1,'float32'));         %Echo 8 Beta Value%
se_hdr = setfield(se_hdr, 'echo8_window', fread(fid,1,'uint16'));       %Echo 8 Window Value%
se_hdr = setfield(se_hdr, 'echo8_level', fread(fid,1,'int16'));        %Echo 8 Level Value%
se_hdr = setfield(se_hdr, 'series_uid', fread(fid,32,'uchar'));    %Series Entity Unique ID%
se_hdr = setfield(se_hdr, 'landmark_uid', fread(fid,32,'uchar'));  %Landmark Unique ID%
se_hdr = setfield(se_hdr, 'equipmnt_uid', fread(fid,32,'uchar'));  %Equipment Unique ID%
se_hdr = setfield(se_hdr, 'se_padding', fread(fid,588,'uchar'));   %Spare Space%

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function im_hdr = ge_readHeaderImage(fid, byte_align)

% returns the image header as a structure.
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)

% define the structure and read in the data
% to overcome the byte alignment problems
% break up the assignment into pieces using the setfield function

im_hdr = struct('im_suid', fread(fid,4,'uchar')); %Suite id for this image
im_hdr = setfield(im_hdr, 'im_uniq', fread(fid,1,'int16'));            %The Make-Unique Flag
im_hdr = setfield(im_hdr, 'im_diskid', fread(fid,1,'uchar'));          %Disk ID for this Image
fseek(fid, 1, 0);% 16-bit alignment
im_hdr = setfield(im_hdr, 'im_exno', fread(fid,1,'uint16'));            %Exam number for this image
im_hdr = setfield(im_hdr, 'im_seno', fread(fid,1,'int16'));            %Series Number for this image
im_hdr = setfield(im_hdr, 'im_no', fread(fid,1,'int16'));              %Image Number
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'im_datetime', fread(fid,1,'int32'));        %Allocation Image date/time stamp
im_hdr = setfield(im_hdr, 'im_actual_dt', fread(fid,1,'int32'));       %Actual Image date/time stamp
im_hdr = setfield(im_hdr, 'sctime', fread(fid,1,'float32'));             %Duration of scan
im_hdr = setfield(im_hdr, 'slthick', fread(fid,1,'float32'));            %Slice Thickness (mm)
im_hdr = setfield(im_hdr, 'imatrix_X', fread(fid,1,'int16'));          %Image matrix size - X
im_hdr = setfield(im_hdr, 'imatrix_Y', fread(fid,1,'int16'));          %Image matrix size - Y
im_hdr = setfield(im_hdr, 'dfov', fread(fid,1,'float32'));               %Display field of view - X (mm)
im_hdr = setfield(im_hdr, 'dfov_rect', fread(fid,1,'float32'));          %Display field of view - Y (if different)
im_hdr = setfield(im_hdr, 'dim_X', fread(fid,1,'float32'));              %Image dimension - X
im_hdr = setfield(im_hdr, 'dim_Y', fread(fid,1,'float32'));              %Image dimension - Y
im_hdr = setfield(im_hdr, 'pixsize_X', fread(fid,1,'float32'));          %Image pixel size - X
im_hdr = setfield(im_hdr, 'pixsize_Y', fread(fid,1,'float32'));          %Image pixel size - Y
im_hdr = setfield(im_hdr, 'pdid', fread(fid,14,'uchar'));          %Pixel Data ID
im_hdr = setfield(im_hdr, 'contrastIV', fread(fid,17,'uchar'));    %IV Contrast Agent
im_hdr = setfield(im_hdr, 'contrastOral', fread(fid,17,'uchar'));  %Oral Contrast Agent
im_hdr = setfield(im_hdr, 'contmode', fread(fid,1,'int16'));           %Image Contrast Mode
im_hdr = setfield(im_hdr, 'serrx', fread(fid,1,'int16'));              %Series from which prescribed
im_hdr = setfield(im_hdr, 'imgrx', fread(fid,1,'int16'));              %Image from which prescribed
im_hdr = setfield(im_hdr, 'screenformat', fread(fid,1,'int16'));       %Screen Format(8/16 bit)
im_hdr = setfield(im_hdr, 'plane', fread(fid,1,'int16'));              %Plane Type
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'scanspacing', fread(fid,1,'float32'));        %Spacing between scans (mm?)
im_hdr = setfield(im_hdr, 'im_compress', fread(fid,1,'int16'));        %Image compression type for allocation
im_hdr = setfield(im_hdr, 'im_scouttype', fread(fid,1,'int16'));       %Scout Type (AP or lateral)
im_hdr = setfield(im_hdr, 'loc_ras', fread(fid,1,'uchar'));            %RAS letter of image location
fseek(fid, 1, 0); % 16-bit alignment
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'loc', fread(fid,1,'float32'));                %Image location
im_hdr = setfield(im_hdr, 'ctr_R', fread(fid,1,'float32'));              %Center R coord of plane image
im_hdr = setfield(im_hdr, 'ctr_A', fread(fid,1,'float32'));              %Center A coord of plane image
im_hdr = setfield(im_hdr, 'ctr_S', fread(fid,1,'float32'));              %Center S coord of plane image
im_hdr = setfield(im_hdr, 'norm_R', fread(fid,1,'float32'));             %Normal R coord
im_hdr = setfield(im_hdr, 'norm_A', fread(fid,1,'float32'));             %Normal A coord
im_hdr = setfield(im_hdr, 'norm_S', fread(fid,1,'float32'));             %Normal S coord
im_hdr = setfield(im_hdr, 'tlhc_R', fread(fid,1,'float32'));             %R Coord of Top Left Hand Corner
im_hdr = setfield(im_hdr, 'tlhc_A', fread(fid,1,'float32'));             %A Coord of Top Left Hand Corner
im_hdr = setfield(im_hdr, 'tlhc_S', fread(fid,1,'float32'));             %S Coord of Top Left Hand Corner
im_hdr = setfield(im_hdr, 'trhc_R', fread(fid,1,'float32'));             %R Coord of Top Right Hand Corner
im_hdr = setfield(im_hdr, 'trhc_A', fread(fid,1,'float32'));             %A Coord of Top Right Hand Corner
im_hdr = setfield(im_hdr, 'trhc_S', fread(fid,1,'float32'));             %S Coord of Top Right Hand Corner
im_hdr = setfield(im_hdr, 'brhc_R', fread(fid,1,'float32'));             %R Coord of Bottom Right Hand Corner
im_hdr = setfield(im_hdr, 'brhc_A', fread(fid,1,'float32'));             %A Coord of Bottom Right Hand Corner
im_hdr = setfield(im_hdr, 'brhc_S', fread(fid,1,'float32'));             %S Coord of Bottom Right Hand Corner
im_hdr = setfield(im_hdr, 'forimgrev', fread(fid,4,'uchar'));      %Foreign Image Revision
im_hdr = setfield(im_hdr, 'tr', fread(fid,1,'int32'));                 %Pulse repetition time(usec)
im_hdr = setfield(im_hdr, 'ti', fread(fid,1,'int32'));                 %Pulse inversion time(usec)
im_hdr = setfield(im_hdr, 'te', fread(fid,1,'int32'));                 %Pulse echo time(usec)
im_hdr = setfield(im_hdr, 'te2', fread(fid,1,'int32'));                %Second echo echo (usec)
im_hdr = setfield(im_hdr, 'numecho', fread(fid,1,'int16'));            %Number of echoes
im_hdr = setfield(im_hdr, 'echonum', fread(fid,1,'int16'));            %Echo Number
im_hdr = setfield(im_hdr, 'tbldlta', fread(fid,1,'float32'));            %Table Delta
im_hdr = setfield(im_hdr, 'nex', fread(fid,1,'float32'));                %Number of Excitations
im_hdr = setfield(im_hdr, 'contig', fread(fid,1,'int16'));             %Continuous Slices Flag
im_hdr = setfield(im_hdr, 'hrtrate', fread(fid,1,'int16'));            %Cardiac Heart Rate (bpm)
im_hdr = setfield(im_hdr, 'tdel', fread(fid,1,'int32'));               %Delay time after trigger (msec)
im_hdr = setfield(im_hdr, 'saravg', fread(fid,1,'float32'));             %Average SAR
im_hdr = setfield(im_hdr, 'sarpeak', fread(fid,1,'float32'));            %Peak SAR
im_hdr = setfield(im_hdr, 'monsar', fread(fid,1,'int16'));             %Monitor SAR flag
im_hdr = setfield(im_hdr, 'trgwindow', fread(fid,1,'int16'));          %Trigger window (% of R-R interval)
im_hdr = setfield(im_hdr, 'reptime', fread(fid,1,'float32'));            %Cardiac repetition time
im_hdr = setfield(im_hdr, 'imgpcyc', fread(fid,1,'int16'));            %Images per cardiac cycle
im_hdr = setfield(im_hdr, 'xmtgain', fread(fid,1,'int16'));            %Actual Transmit Gain (.1 db)
im_hdr = setfield(im_hdr, 'rcvgain1', fread(fid,1,'int16'));           %Actual Receive Gain Analog (.1 db)
im_hdr = setfield(im_hdr, 'rcvgain2', fread(fid,1,'int16'));           %Actual Receive Gain Digital (.1 db)
im_hdr = setfield(im_hdr, 'mr_flip', fread(fid,1,'int16'));            %Flip Angle for GRASS scans (deg.)
if byte_align; fseek(fid, 2, 0); end % byte alignment
im_hdr = setfield(im_hdr, 'mindat', fread(fid,1,'int32'));             %Minimum Delay after Trigger (uSec)
im_hdr = setfield(im_hdr, 'cphase', fread(fid,1,'int16'));             %Total Cardiac Phase prescribed
im_hdr = setfield(im_hdr, 'swappf', fread(fid,1,'int16'));             %Swap Phase/Frequency Axis
im_hdr = setfield(im_hdr, 'pauseint', fread(fid,1,'int16'));           %Pause Interval (slices)
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'pausetime', fread(fid,1,'float32'));          %Pause Time
im_hdr = setfield(im_hdr, 'obplane', fread(fid,1,'int32'));            %Oblique Plane
im_hdr = setfield(im_hdr, 'slocfov', fread(fid,1,'int32'));            %Slice Offsets on Freq axis
im_hdr = setfield(im_hdr, 'xmtfreq', fread(fid,1,'int32'));            %Center Frequency (0.1 Hz)
im_hdr = setfield(im_hdr, 'autoxmtfreq', fread(fid,1,'int32'));        %Auto Center Frequency (0.1 Hz)
im_hdr = setfield(im_hdr, 'autoxmtgain', fread(fid,1,'int16'));        %Auto Transmit Gain (0.1 dB)
im_hdr = setfield(im_hdr, 'prescan_r1', fread(fid,1,'int16'));         %PreScan R1 - Analog
im_hdr = setfield(im_hdr, 'prescan_r2', fread(fid,1,'int16'));         %PreScan R2 - Digital
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'user_bitmap', fread(fid,1,'int32'));        %Bitmap defining user CVs
im_hdr = setfield(im_hdr, 'cenfreq', fread(fid,1,'int16'));            %Center Frequency Method
im_hdr = setfield(im_hdr, 'imode', fread(fid,1,'int16'));              %Imaging Mode
im_hdr = setfield(im_hdr, 'iopt', fread(fid,1,'int32'));               %Imaging Options
im_hdr = setfield(im_hdr, 'pseq', fread(fid,1,'int16'));               %Pulse Sequence
im_hdr = setfield(im_hdr, 'pseqmode', fread(fid,1,'int16'));           %Pulse Sequence Mode
im_hdr = setfield(im_hdr, 'psdname', fread(fid,33,'uchar'));       %Pulse Sequence Name
fseek(fid,1,0); % 16-bit alignment
if byte_align; fseek(fid, 2, 0); end % byte alignment
im_hdr = setfield(im_hdr, 'psd_datetime', fread(fid,1,'int32'));       %PSD Creation Date and Time
im_hdr = setfield(im_hdr, 'psd_iname', fread(fid,13,'uchar'));     %PSD name from inside PSD
fseek(fid, 1, 0); % 16-bit alignment
im_hdr = setfield(im_hdr, 'ctyp', fread(fid,1,'int16'));               %Coil Type
im_hdr = setfield(im_hdr, 'cname', fread(fid,17,'uchar'));         %Coil Name
fseek(fid, 1, 0); % 16-bit alignment
im_hdr = setfield(im_hdr, 'surfctyp', fread(fid,1,'int16'));           %Surface Coil Type
im_hdr = setfield(im_hdr, 'surfcext', fread(fid,1,'int16'));           %Extremity Coil Flag
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'rawrunnum', fread(fid,1,'int32'));          %RawData Run Number
im_hdr = setfield(im_hdr, 'cal_fldstr', fread(fid,1,'uint32'));         %Calibrated Field Strength (x10 uGauss)
im_hdr = setfield(im_hdr, 'supp_tech', fread(fid,1,'int16'));          %SAT fat/water/none
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'vbw', fread(fid,1,'float32'));                %Variable Bandwidth (Hz)
im_hdr = setfield(im_hdr, 'slquant', fread(fid,1,'int16'));            %Number of slices in this scan group
im_hdr = setfield(im_hdr, 'gpre', fread(fid,1,'int16'));               %Graphically prescribed
im_hdr = setfield(im_hdr, 'intr_del', fread(fid,1,'int32'));           %Interimage/interloc delay (uSec)
im_hdr = setfield(im_hdr, 'user0', fread(fid,1,'float32'));              %User Variable 0
im_hdr = setfield(im_hdr, 'user1', fread(fid,1,'float32'));              %User Variable 1
im_hdr = setfield(im_hdr, 'user2', fread(fid,1,'float32'));              %User Variable 2
im_hdr = setfield(im_hdr, 'user3', fread(fid,1,'float32'));              %User Variable 3
im_hdr = setfield(im_hdr, 'user4', fread(fid,1,'float32'));              %User Variable 4
im_hdr = setfield(im_hdr, 'user5', fread(fid,1,'float32'));              %User Variable 5
im_hdr = setfield(im_hdr, 'user6', fread(fid,1,'float32'));              %User Variable 6
im_hdr = setfield(im_hdr, 'user7', fread(fid,1,'float32'));              %User Variable 7
im_hdr = setfield(im_hdr, 'user8', fread(fid,1,'float32'));              %User Variable 8
im_hdr = setfield(im_hdr, 'user9', fread(fid,1,'float32'));              %User Variable 9
im_hdr = setfield(im_hdr, 'user10', fread(fid,1,'float32'));             %User Variable 10
im_hdr = setfield(im_hdr, 'user11', fread(fid,1,'float32'));             %User Variable 11
im_hdr = setfield(im_hdr, 'user12', fread(fid,1,'float32'));             %User Variable 12
im_hdr = setfield(im_hdr, 'user13', fread(fid,1,'float32'));             %User Variable 13
im_hdr = setfield(im_hdr, 'user14', fread(fid,1,'float32'));             %User Variable 14
im_hdr = setfield(im_hdr, 'user15', fread(fid,1,'float32'));             %User Variable 15
im_hdr = setfield(im_hdr, 'user16', fread(fid,1,'float32'));             %User Variable 16
im_hdr = setfield(im_hdr, 'user17', fread(fid,1,'float32'));             %User Variable 17
im_hdr = setfield(im_hdr, 'user18', fread(fid,1,'float32'));             %User Variable 18
im_hdr = setfield(im_hdr, 'user19', fread(fid,1,'float32'));             %User Variable 19
im_hdr = setfield(im_hdr, 'user20', fread(fid,1,'float32'));             %User Variable 20
im_hdr = setfield(im_hdr, 'user21', fread(fid,1,'float32'));             %User Variable 21
im_hdr = setfield(im_hdr, 'user22', fread(fid,1,'float32'));             %User Variable 22
im_hdr = setfield(im_hdr, 'user23', fread(fid,1,'float32'));             %Projection Angle
im_hdr = setfield(im_hdr, 'user24', fread(fid,1,'float32'));             %Concat Sat Type Flag
im_hdr = setfield(im_hdr, 'im_alloc_key', fread(fid,13,'uchar'));
fseek(fid, 1, 0); % 16-bit alignment
if byte_align; fseek(fid, 2, 0); end % 32-bit alignment
im_hdr = setfield(im_hdr, 'im_lastmod', fread(fid,1,'int32'));         %Date/Time of Last Change
im_hdr = setfield(im_hdr, 'im_verscre', fread(fid,2,'uchar'));     %Genesis Version - Created
im_hdr = setfield(im_hdr, 'im_verscur', fread(fid,2,'uchar'));     %Genesis Version - Now
im_hdr = setfield(im_hdr, 'im_pds_a', fread(fid,1,'int32'));           %PixelData size - as stored
im_hdr = setfield(im_hdr, 'im_pds_c', fread(fid,1,'int32'));           %PixelData size - Compressed
im_hdr = setfield(im_hdr, 'im_pds_u', fread(fid,1,'int32'));           %PixelData size - UnCompressed
im_hdr = setfield(im_hdr, 'im_checksum', fread(fid,1,'int32'));        %AcqRecon record checksum
im_hdr = setfield(im_hdr, 'im_archived', fread(fid,1,'int32'));        %Image Archive Flag
im_hdr = setfield(im_hdr, 'im_complete', fread(fid,1,'int32'));        %Image Complete Flag
im_hdr = setfield(im_hdr, 'satbits', fread(fid,1,'int16'));            %Bitmap of SAT selections
im_hdr = setfield(im_hdr, 'scic', fread(fid,1,'int16'));               %Surface Coil Intensity Correction Flag
im_hdr = setfield(im_hdr, 'satxloc1', fread(fid,1,'int16'));           %R-side SAT pulse loc rel to lndmrk
im_hdr = setfield(im_hdr, 'satxloc2', fread(fid,1,'int16'));           %L-side SAT pulse loc rel to lndmrk
im_hdr = setfield(im_hdr, 'satyloc1', fread(fid,1,'int16'));           %A-side SAT pulse loc rel to lndmrk
im_hdr = setfield(im_hdr, 'satyloc2', fread(fid,1,'int16'));           %P-side SAT pulse loc rel to lndmrk
im_hdr = setfield(im_hdr, 'satzloc1', fread(fid,1,'int16'));           %S-side SAT pulse loc rel to lndmrk
im_hdr = setfield(im_hdr, 'satzloc2', fread(fid,1,'int16'));           %I-side SAT pulse loc rel to lndmrk
im_hdr = setfield(im_hdr, 'satxthick', fread(fid,1,'int16'));          %Thickness of X-axis SAT pulse
im_hdr = setfield(im_hdr, 'satythick', fread(fid,1,'int16'));          %Thickness of Y-axis SAT pulse
im_hdr = setfield(im_hdr, 'satzthick', fread(fid,1,'int16'));          %Thickness of Z-axis SAT pulse
im_hdr = setfield(im_hdr, 'flax', fread(fid,1,'int16'));               %Phase contrast flow axis
im_hdr = setfield(im_hdr, 'venc', fread(fid,1,'int16'));               %Phase contrast velocity encoding
im_hdr = setfield(im_hdr, 'thk_disclmr', fread(fid,1,'int16'));        %Slice Thickness
im_hdr = setfield(im_hdr, 'ps_flag', fread(fid,1,'int16'));            %Auto/Manual Prescan flag
im_hdr = setfield(im_hdr, 'ps_status', fread(fid,1,'int16'));          %Bitmap of changed values
im_hdr = setfield(im_hdr, 'image_type', fread(fid,1,'int16'));         %Magnitude, Phase, Imaginary, or Real
im_hdr = setfield(im_hdr, 'vas_collapse', fread(fid,1,'int16'));       %Collapse Image
im_hdr = setfield(im_hdr, 'user23n', fread(fid,1,'float32'));            %User Variable 23
im_hdr = setfield(im_hdr, 'user24n', fread(fid,1,'float32'));            %User Variable 24
im_hdr = setfield(im_hdr, 'proj_alg', fread(fid,1,'int16'));           %Projection Algorithm
im_hdr = setfield(im_hdr, 'proj_name', fread(fid,13,'uchar'));     %Projection Algorithm Name
fseek(fid, 1, 0); % 16-bit alignment
im_hdr = setfield(im_hdr, 'x_axis_rot', fread(fid,1,'float32'));         %X Axis Rotation
im_hdr = setfield(im_hdr, 'y_axis_rot', fread(fid,1,'float32'));         %Y Axis Rotation
im_hdr = setfield(im_hdr, 'z_axis_rot', fread(fid,1,'float32'));         %Z Axis Rotation
im_hdr = setfield(im_hdr, 'thresh_min1', fread(fid,1,'int32'));        %Lower Range of Pixels 1
im_hdr = setfield(im_hdr, 'thresh_max1', fread(fid,1,'int32'));        %Upper Range of Pixels 1
im_hdr = setfield(im_hdr, 'thresh_min2', fread(fid,1,'int32'));        %Lower Range of Pixels 2
im_hdr = setfield(im_hdr, 'thresh_max2', fread(fid,1,'int32'));        %Upper Range of Pixels 2
im_hdr = setfield(im_hdr, 'echo_trn_len', fread(fid,1,'int16'));       %Echo Train Length for Fast Spin Echo
im_hdr = setfield(im_hdr, 'frac_echo', fread(fid,1,'int16'));          %Fractional Echo - Effective TE Flag
im_hdr = setfield(im_hdr, 'prep_pulse', fread(fid,1,'int16'));         %Preporatory Pulse Option
im_hdr = setfield(im_hdr, 'cphasenum', fread(fid,1,'int16'));          %Cardiac Phase Number
im_hdr = setfield(im_hdr, 'var_echo', fread(fid,1,'int16'));           %Variable Echo Flag
im_hdr = setfield(im_hdr, 'ref_img', fread(fid,1,'uchar'));            %Reference Image Field
im_hdr = setfield(im_hdr, 'sum_img', fread(fid,1,'uchar'));            %Summary Image Field
im_hdr = setfield(im_hdr, 'img_window', fread(fid,1,'uint16'));         %Window Value
im_hdr = setfield(im_hdr, 'img_level', fread(fid,1,'int16'));          %Level Value
im_hdr = setfield(im_hdr, 'slop_int_1', fread(fid,1,'int32'));         %Number of 3D Slabs
im_hdr = setfield(im_hdr, 'slop_int_2', fread(fid,1,'int32'));         %Slice Locs Per 3D Slab
im_hdr = setfield(im_hdr, 'slop_int_3', fread(fid,1,'int32'));         %# of Slice Locs on Each Slab Which Overlap Neighbors
im_hdr = setfield(im_hdr, 'slop_int_4', fread(fid,1,'int32'));         %Image Filtering 0.5/0.2T
im_hdr = setfield(im_hdr, 'slop_int_5', fread(fid,1,'int32'));         %Integer Slop Field 5
im_hdr = setfield(im_hdr, 'slop_float_1', fread(fid,1,'float32'));       %Float Slop Field 1
im_hdr = setfield(im_hdr, 'slop_float_2', fread(fid,1,'float32'));       %Float Slop Field 2
im_hdr = setfield(im_hdr, 'slop_float_3', fread(fid,1,'float32'));       %Float Slop Field 3
im_hdr = setfield(im_hdr, 'slop_float_4', fread(fid,1,'float32'));       %Float Slop Field 4
im_hdr = setfield(im_hdr, 'slop_float_5', fread(fid,1,'float32'));       %Float Slop Field 5
im_hdr = setfield(im_hdr, 'slop_str_1', fread(fid,16,'uchar'));    %String Slop Field 1
im_hdr = setfield(im_hdr, 'slop_str_2', fread(fid,16,'uchar'));    %String Slop Field 2
im_hdr = setfield(im_hdr, 'scanactno', fread(fid,1,'int16'));          %Scan Acquisition Number
im_hdr = setfield(im_hdr, 'vasflags', fread(fid,1,'int16'));           %Magnitude Weighting Flag
im_hdr = setfield(im_hdr, 'vencscale', fread(fid,1,'float32'));          %Scale Weighted Venc
im_hdr = setfield(im_hdr, 'integrity', fread(fid,1,'int16'));          %GE Image Integrity
if byte_align; fseek(fid, 2, 0); end % byte alignment
im_hdr = setfield(im_hdr, 'fphase', fread(fid,1,'int32'));             %Number Of Phases
im_hdr = setfield(im_hdr, 'freq_dir', fread(fid,1,'int16'));           %Frequency Direction
im_hdr = setfield(im_hdr, 'vas_mode', fread(fid,1,'int16'));           %Vascular Mode
im_hdr = setfield(im_hdr, 'image_uid', fread(fid,32,'uchar'));     %Image Unique ID
im_hdr = setfield(im_hdr, 'sop_uid', fread(fid,32,'uchar'));       %Service Obj Class Unique ID
im_hdr = setfield(im_hdr, 'dont_use_1', fread(fid,1,'int16'));         %This field is not used
im_hdr = setfield(im_hdr, 'dont_use_2', fread(fid,1,'int16'));         %This field is not used
im_hdr = setfield(im_hdr, 'dont_use_3', fread(fid,1,'int16'));         %This field is not used
im_hdr = setfield(im_hdr, 'pscopts', fread(fid,1,'int16'));            %bitmap of prescan options
im_hdr = setfield(im_hdr, 'asoffsetx', fread(fid,1,'int16'));          %gradient offset in X-direction
im_hdr = setfield(im_hdr, 'asoffsety', fread(fid,1,'int16'));          %gradient offset in Y-direction
im_hdr = setfield(im_hdr, 'asoffsetz', fread(fid,1,'int16'));          %gradient offset in Z-direction
im_hdr = setfield(im_hdr, 'unoriginal', fread(fid,1,'int16'));         %identifies image as original or unoriginal
im_hdr = setfield(im_hdr, 'interleaves', fread(fid,1,'int16'));        %number of EPI shots
im_hdr = setfield(im_hdr, 'effechospace', fread(fid,1,'int16'));       %effective echo spacing for EPI
im_hdr = setfield(im_hdr, 'viewsperseg', fread(fid,1,'int16'));        %views per segment
im_hdr = setfield(im_hdr, 'rbpm', fread(fid,1,'int16'));               %respiratory rate, breaths per min
im_hdr = setfield(im_hdr, 'rtpoint', fread(fid,1,'int16'));            %respiratory trigger point as percent of max.
im_hdr = setfield(im_hdr, 'rcvrtype', fread(fid,1,'int16'));           %type of receiver used
im_hdr = setfield(im_hdr, 'dbdt', fread(fid,1,'float32'));               %peak rate of change of gradient field, tesla/sec
im_hdr = setfield(im_hdr, 'dbdtper', fread(fid,1,'float32'));            %limit in units of percent of theoretical curve
im_hdr = setfield(im_hdr, 'estdbdtper', fread(fid,1,'float32'));         %PSD estimated limit in units of percent
im_hdr = setfield(im_hdr, 'estdbdtts', fread(fid,1,'float32'));          %PSD estimated limit in Teslas/sec
im_hdr = setfield(im_hdr, 'saravghead', fread(fid,1,'float32'));         %Avg head SAR
im_hdr = setfield(im_hdr, 'neg_scanspacing', fread(fid,1,'float32'));    %Negative scan spacing for overlap slices
im_hdr = setfield(im_hdr, 'offsetfreq', fread(fid,1,'int32'));         %Offset Frequency - Mag.Transfer
im_hdr = setfield(im_hdr, 'user_usage_tag', fread(fid,1,'uint32'));     %Defines how following user CVs are to be filled in
%Default value = 0x00000000
%GE range = 0x00000001 - 0x7fffffff
%Research = 0x80000000 - 0xffffffff
im_hdr = setfield(im_hdr, 'user_fill_mapMSW', fread(fid,1,'uint32'));   %Define what process fills in the user CVs, ifcc or TIR
im_hdr = setfield(im_hdr, 'user_fill_mapLSW', fread(fid,1,'uint32'));   %Define what process fills in the user CVs, ifcc or TIR
im_hdr = setfield(im_hdr, 'user25', fread(fid,1,'float32'));             %User Variable 25
im_hdr = setfield(im_hdr, 'user26', fread(fid,1,'float32'));             %User Variable 26
im_hdr = setfield(im_hdr, 'user27', fread(fid,1,'float32'));             %User Variable 27
im_hdr = setfield(im_hdr, 'user28', fread(fid,1,'float32'));             %User Variable 28
im_hdr = setfield(im_hdr, 'user29', fread(fid,1,'float32'));             %User Variable 29
im_hdr = setfield(im_hdr, 'user30', fread(fid,1,'float32'));             %User Variable 30
im_hdr = setfield(im_hdr, 'user31', fread(fid,1,'float32'));             %User Variable 31
im_hdr = setfield(im_hdr, 'user32', fread(fid,1,'float32'));             %User Variable 32
im_hdr = setfield(im_hdr, 'user33', fread(fid,1,'float32'));             %User Variable 33
im_hdr = setfield(im_hdr, 'user34', fread(fid,1,'float32'));             %User Variable 34
im_hdr = setfield(im_hdr, 'user35', fread(fid,1,'float32'));             %User Variable 35
im_hdr = setfield(im_hdr, 'user36', fread(fid,1,'float32'));             %User Variable 36
im_hdr = setfield(im_hdr, 'user37', fread(fid,1,'float32'));             %User Variable 37
im_hdr = setfield(im_hdr, 'user38', fread(fid,1,'float32'));             %User Variable 38
im_hdr = setfield(im_hdr, 'user39', fread(fid,1,'float32'));             %User Variable 39
im_hdr = setfield(im_hdr, 'user40', fread(fid,1,'float32'));             %User Variable 40
im_hdr = setfield(im_hdr, 'user41', fread(fid,1,'float32'));             %User Variable 41
im_hdr = setfield(im_hdr, 'user42', fread(fid,1,'float32'));             %User Variable 42
im_hdr = setfield(im_hdr, 'user43', fread(fid,1,'float32'));             %User Variable 43
im_hdr = setfield(im_hdr, 'user44', fread(fid,1,'float32'));             %User Variable 44
im_hdr = setfield(im_hdr, 'user45', fread(fid,1,'float32'));             %User Variable 45
im_hdr = setfield(im_hdr, 'user46', fread(fid,1,'float32'));             %User Variable 46
im_hdr = setfield(im_hdr, 'user47', fread(fid,1,'float32'));             %User Variable 47
im_hdr = setfield(im_hdr, 'user48', fread(fid,1,'float32'));             %User Variable 48
im_hdr = setfield(im_hdr, 'slop_int_6', fread(fid,1,'int32'));         %Integer Slop Field 6
im_hdr = setfield(im_hdr, 'slop_int_7', fread(fid,1,'int32'));         %Integer Slop Field 7
im_hdr = setfield(im_hdr, 'slop_int_8', fread(fid,1,'int32'));         %Integer Slop Field 8
im_hdr = setfield(im_hdr, 'slop_int_9', fread(fid,1,'int32'));         %Integer Slop Field 9
im_hdr = setfield(im_hdr, 'mr_padding', fread(fid,32,'uchar'));    %Spare Space

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ge = ge_readHeaderPixel(fid, byte_align, ge)

% returns the pixel header as a structure.
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)

% read magic number
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_magic', fread(fid,1,'int32'));
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_hdr_length', fread(fid,1,'int32'));% length of total header in bytes and
% a byte displacement to the 'pixel data area'  
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_width', fread(fid,1,'int32'));  % width (pixels) of image  
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_height', fread(fid,1,'int32')); % height (pixels) of image 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_depth', fread(fid,1,'int32'));  % depth (1, 8, 16, or 24 bits) of pixel 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_compress', fread(fid,1,'int32'));   % type of compression; see IC_* below 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_dwindow', fread(fid,1,'int32'));    % default window setting 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_dlevel', fread(fid,1,'int32')); % default level setting 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_bgshade', fread(fid,1,'int32'));    % background shade to use for non-image 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_ovrflow', fread(fid,1,'int32'));    % overflow value 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_undflow', fread(fid,1,'int32'));    % underflow value 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_top_offset', fread(fid,1,'int32')); % number of blank lines at image top 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_bot_offset', fread(fid,1,'int32')); % number of blank lines at image bottom 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_version', fread(fid,1,'int16'));    % version of the header structure
if byte_align, fseek(fid,2,0); end % 32-bit alignment
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_checksum', fread(fid,1,'uint16')); % 16 bit end_around_carry sum of pixels 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_id', fread(fid,1,'int32'));   % a byte disp to unique image identifier 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_id', fread(fid,1,'int32'));   % byte length of unique image identifier 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_unpack', fread(fid,1,'int32'));   % a byte disp to 'unpack control' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_unpack', fread(fid,1,'int32'));   % byte length of 'unpack control' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_compress', fread(fid,1,'int32')); % a byte disp to 'compression control' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_compress', fread(fid,1,'int32')); % byte length of 'compression control' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_histo', fread(fid,1,'int32'));    % a byte disp to 'histogram control' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_histo', fread(fid,1,'int32'));    % byte length of 'histogram control' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_text', fread(fid,1,'int32')); % a byte disp to 'text plane data' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_text', fread(fid,1,'int32')); % byte length of 'text plane data' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_graphics', fread(fid,1,'int32')); % a byte disp to 'graphics plane data' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_graphics', fread(fid,1,'int32')); % byte length of 'graphics plane data' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_dbHdr', fread(fid,1,'int32'));    % a byte disp to 'data base header data' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_dbHdr', fread(fid,1,'int32'));    % byte length of 'data base header data' 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_levelOffset', fread(fid,1,'int32'));% value to add to stored Pixel Data values
% to get the correct presentation value 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_p_user', fread(fid,1,'int32')); % byte displacement to user defined data 
ge.hdr.pixel = setfield(ge.hdr.pixel, 'img_l_user', fread(fid,1,'int32')); % byte length of user defined data 

return
