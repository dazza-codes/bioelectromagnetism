function [ avw, machine ] = avw_read(fileprefix,IMGorient,machine,verbose)

% avw_read - read Analyze format data image (*.img)
% 
% [ avw, machine ] = avw_read([fileprefix], [orient], [machine], [verbose])
% 
% This function calls avw_hdr_read and avw_img_read, providing
% just a slightly easier command line function call.
% 
% fileprefix - a string, the filename without the .img extension.
%              A gui prompt appears if this argument is missing.
% 
% orient - read a specified orientation, integer values:
% 
%          '', use header history orient field
%          0,  transverse unflipped (LAS*)
%          1,  coronal unflipped (LA*S)
%          2,  sagittal unflipped (L*AS)
%          3,  transverse flipped (LPS*)
%          4,  coronal flipped (LA*I)
%          5,  sagittal flipped (L*AI)
% 
% where * follows the slice dimension and letters indicate +XYZ
% orientations (L left, R right, A anterior, P posterior,
% I inferior, & S superior).
% 
% It is rare, but a dataset might be stored in the 3-5 
% orientations. For more information about orientation, see the
% documentation in the mri_toolbox doc directory and the
% extensive comments in avw_img_read and avw_hdr_read.  
% See also the avw_flip function for any orthogonal 
% reorientation (although this should not be necessary
% and must be done with great care, try not to invalidate
% the Analyze orientation specifications).
% 
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-le' but the routine
%           will automatically switch between little and big
%           endian to read any such Analyze header.  It
%           reports the appropriate machine format and can
%           return the machine value.
% 
% verbose - the default is to output processing information to the command
%           window.  If verbose = 0, this will not happen.
%
% Returned values:
% 
% avw.hdr - a struct with image data parameters.
% avw.img - a 3D matrix of image data (double precision).
% 
% All going well, the returned 3D matrix in avw.img will 
% correspond with the default ANALYZE coordinate system, 
% which is a Left-handed coordinate system (radiological
% orientation):
% 
% X-Y plane is Transverse/Axial
% X-Z plane is Coronal
% Y-Z plane is Sagittal
% 
% X axis runs from patient Right (low X) to patient Left (high X)
% Y axis runs from Posterior (low Y) to Anterior (high Y)
% Z axis runs from Inferior (low Z) to Superior (high Z)
% 
% See also: avw_hdr_read, avw_img_read (both called by this function), 
%           avw_write, avw_img_write, avw_hdr_write, 
%           avw_view, avw_flip
% 


% $Revision: 1.2 $ $Date: 2006/06/17 02:47:56 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2003, Darren.Weber_at_radiology.ucsf.edu
%                    The Analyze format is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%                    - created this wrapper for avw_img_read, see
%                      that function for extensive comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~exist('fileprefix','var'),
    [fileprefix, pathname, filterindex] = uigetfile('*.hdr','locate an Analyze .hdr file');
    if pathname, cd(pathname); end
    if ~fileprefix,
      error('no .hdr file specified');
    end
  end
  if isempty(fileprefix),
    [fileprefix, pathname, filterindex] = uigetfile('*.hdr','locate an Analyze .hdr file');
    if pathname, cd(pathname); end
    if ~fileprefix,
      error('no .hdr file specified');
    end
  end

  if ~exist('IMGorient','var'), IMGorient = ''; end
  if ~exist('machine','var'), machine = 'ieee-le'; end
  if ~exist('verbose','var'), verbose = 1; end

  if isempty(IMGorient), IMGorient = ''; end
  if isempty(machine), machine = 'ieee-le'; end
  if isempty(verbose), verbose = 1; end

  fileprefix = strrep(fileprefix,'.hdr','');
  fileprefix = strrep(fileprefix,'.img','');

  % MAIN

  fid = fopen(sprintf('%s.img',fileprefix),'r',machine);
  if fid < 0,
    msg = sprintf('...cannot open file %s.img\n\n',fileprefix);
    error(msg);
  else
    fclose(fid);
    % avw_img_read will call avw_hdr_read also
    avw = avw_img_read(fileprefix,IMGorient,machine,verbose);
  end

  return
