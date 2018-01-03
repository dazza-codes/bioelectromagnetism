function [ avw ] = avw_hdr_check_datatype(avw,verbose)

% avw_hdr_check_datatype - read Analyze format data header (*.hdr)
%
% avw = avw_hdr_check_datatype(avw,verbose)
%
% attempts to set the datatype based on the bits per pixel, although
% this really needs to be done the other way around.  The following
% table indicates the values recognised for 
% avw.hdr.dime.datatype and avw.hdr.dime.bitpix
%
% short int datatype      /* Datatype for this image set */ 
% /*Acceptable values for datatype are*/ 
% #define DT_NONE             0
% #define DT_UNKNOWN          0    /*Unknown data type*/ 
% #define DT_BINARY           1    /*Binary             ( 1 bit per voxel)*/ 
% #define DT_UNSIGNED_CHAR    2    /*Unsigned character ( 8 bits per voxel)*/ 
% #define DT_SIGNED_SHORT     4    /*Signed short       (16 bits per voxel)*/ 
% #define DT_SIGNED_INT       8    /*Signed integer     (32 bits per voxel)*/ 
% #define DT_FLOAT           16    /*Floating point     (32 bits per voxel)*/ 
% #define DT_COMPLEX         32    /*Complex (64 bits per voxel; 2 floating point numbers)/* 
% #define DT_DOUBLE          64    /*Double precision   (64 bits per voxel)*/ 
% #define DT_RGB            128    /*A Red-Green-Blue datatype*/
% #define DT_ALL            255    /*Undocumented*/
% 
% short int bitpix;    /* Number of bits per pixel; 1, 8, 16, 32, or 64. */ 
% 
% verbose - the default is to output processing information to the command
%           window.  If verbose = 0, this will not happen.
%
% See also avw_hdr_write, avw_hdr_make, avw_view_hdr, avw_view
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze format and c code below is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('verbose','var'), verbose = 1; end

if verbose,
    version = '[$Revision: 1.1 $]';
    fprintf('\nAVW_HDR_CHECK_DATATYPE [v%s]\n',version(12:16));  tic;
end

if ~exist('avw','var'),
    error('...no input avw');
end


switch avw.hdr.dime.bitpix,
  
  case 0,
    error('avw.hdr.dime.bitpix = 0, unknown datatype');
    
  case 1,
    avw.hdr.dime.datatype = 1;
    
  case 8,
    avw.hdr.dime.datatype = 2;
    
  case 16,
    avw.hdr.dime.datatype = 4;
    
  case 32,
    warning('bitpix = 32, assuming datatype is float (rather than signed int)');
    avw.hdr.dime.datatype = 16;
    
  case 64,
    warning('bitpix = 64, assuming datatype is double (rather than complex)');
    avw.hdr.dime.datatype = 64;
    
  case 128,
    avw.hdr.dime.datatype = 128;
    
  otherwise
    error('unknown bitpix and datatype');
    
end

if verbose,
    t = toc; fprintf('...done ( %6.2f sec)\n\n',t);
end

return
