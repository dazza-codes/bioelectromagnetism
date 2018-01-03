function [XYZm] = avw_metric(hdr,XYZimg)

% AVW_METRIC - Convert image location from CRS to meters
%
% Useage:  XYZm = avw_metric(hdr,XYZimg)
%
% hdr    - avw.hdr from avw_hdr_read
% XYZimg - Nx3 matrix of image coordinates
% XYZm   - Nx3 matrix of image meter coordinates (not mm)
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  06/2002, Darren.Weber@flinders.edu.au
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(XYZimg), error('XYZimg is empty'); end;
if isempty(hdr),    error('hdr is empty'); end;

if size(XYZimg,2) ~= 3,
    msg = sprintf('AVW_METRIC: XYZimg must be Nx3 matrix\n');
    error(msg);
end

XYZpixdim = double(hdr.dime.pixdim(2:4));

if findstr(hdr.dime.vox_units,'mm'),
    fprintf('AVW_METRIC: voxel units: mm\n');
    XYZpixdim = XYZpixdim ./ 1000;
end
if findstr(hdr.dime.vox_units,'cm'),
    fprintf('AVW_METRIC: voxel units: cm\n');
    XYZpixdim = XYZpixdim ./ 100;
end

XYZpixdim = repmat(XYZpixdim,size(XYZimg,1),1);


XYZm = XYZimg .* XYZpixdim;

return
