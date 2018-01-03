function [reg] = emse_read_reg(file)

% emse_read_reg - Read EMSE/MRVU coregistration matrices
%
% [reg] = emse_read_reg(file)
%
% reg is a struct with the following fields:
%
% reg.translation - the translation in meters along the
%                   x, y and z axes respectively, from
%                   the MRI image frame to head/elec frame.
%
% reg.rotation - The rotation vector contains the angles
%                (in radians) about the x, y and z axes,
%                also from the MRI image frame to the
%                head/elec frame.
%
% reg.elec2mri - 'HeadToImageMatrix' is the 4 x 4 matrix
%                containing the electrode to MRI translation and
%                rotation transformations in homogeneous coordinates:
%                * the upper left 3 x 3 submatrix is rotations
%                  around z, y, x in that order;
%                * the rightmost 3 x 1 column is a projection
%                  vector (all zeros here);
%                * the bottom 1 x 3 row is a translation vector,
%                  equal to -1 * reg.translation here; and
%                * the bottom right (1 x 1) scalar is the
%                  homogenous scale unit, usually 1
%
% reg.mri2elec - 'ImageToHeadMatrix' is the inverse of elec2mri,
%                ie, reg.mri2elec = inv(reg.elec2mri).
%
% This function also reads the fiducial points and the electrode
% coordinates from the registration file, they are returned into:
% reg.RPA, reg.LPA, reg.NAS, reg.Helec, and reg.Melec.  Each of
% the fiducial structs (RPA,LPA,NAS) contains the electrode
% fiducials in the head space (Hh) and the MRI space (Hm), plus the
% MRI fiducials in the head space (Mh) and the MRI space (Mm).
%
% The transformation matrices (T) multiply a column vector, so that
% [x', y', z', 1] = [x, y, z, 1] * T;
% where x',y',z' are in the other coordinate system. For example,
% MRI coordinates into head space:
% tmp = [ reg.Melec ones(size(reg.Melec,1),1) ] * reg.mri2elec;
% Note reg.Helec ~= tmp(:,1:3) due to floating point rounding only.
% Similarly, head space (electrodes) into MRI coordinates:
% tmp = [ reg.Helec ones(size(reg.Helec,1),1) ] * reg.elec2mri;
% Note reg.Melec ~= tmp(:,1:3) due to floating point rounding only.
%
% EMSE Note: The origin in the head frame is at or near the center of
% the skull, while the origin in the image frame is located at the
% bottom right front corner of the bounding box (and so would be
% located at the upper left corner of the first axial slice as
% displayed by MR Viewer).
%
% A useful chapter on homogeneous coordinates, among other things,
% may be found in Mortenson, M. (1985, Chpt. 8), Geometric Modelling,
% New York: John Wiley & Sons.
%


% $Revision: 1.2 $ $Date: 2007/10/27 01:10:04 $

% Licence:  GNU GPL, no express or implied warranties
% History:  06/2002, Darren.Weber@flinders.edu.au
%           09/2002, Darren.Weber@flinders.edu.au
%                    - transposed HeadToImageMatrix so it
%                      can be used as described above
%                    - added reading of most other fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('file','var'),
    fprintf('No input file - see help open_emse_reg\n');
    return;
end

[path,name,ext] = fileparts(file);
file = fullfile(path,[name ext]);

[fid,msg] = fopen(file,'r');
if ~isempty(msg), error(msg); end

fprintf('EMSE_READ_REG: Reading registration data...');
tic
fid = fopen(file,'r','ieee-le');
reg = read_reg(fid);
fclose(fid);
t = toc;
fprintf('done (%6.2f sec).\n',t);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reg] = read_reg(fid)

while 1,
    text = fgetl(fid);
    if ~ischar(text), break, end

    if strmatch('Offset',text),
        % Offset is the translation in meters along the x, y and z axes
        % respectively, from the MRI image frame to head/elec frame.
        text = strrep(text,sscanf(text,'%c',8),'');
        text = strrep(text,']','');
        text = strrep(text,',','');
        reg.translation = sscanf(text,'%f')';
    end
    if strmatch('Rotation',text),
        % The Rotation vector contains the angles (in radians) about
        % the x, y and z axes, also from the MRI image frame to the
        % head/elec frame.
        text = strrep(text,sscanf(text,'%c',10),'');
        text = strrep(text,']','');
        text = strrep(text,',','');
        reg.rotation = sscanf(text,'%f')';
    end
    if strmatch('HeadToImageMatrix',text),
        reg.elec2mri = zeros(4,4);
        for i=1:4,
            text = fgetl(fid);
            reg.elec2mri(i,:) = sscanf(text,'%f')';
        end
        % The emse matrix requires transposition
        reg.elec2mri = reg.elec2mri';
        % It is more accurate to do this:
        reg.mri2elec = inv(reg.elec2mri);
    end
    % See inverse calculation above to get reg.mri2elec
    %if strmatch('ImageToHeadMatrix',text),
    %    reg.mri2elec = zeros(4,4);
    %    for i=1:4,
    %        text = fgetl(fid);
    %        reg.mri2elec(i,:) = sscanf(text,'%f')';
    %    end
    %    % The emse matrix requires transposition
    %    reg.mri2elec = reg.mri2elec';
    %end

    % The coordinates of the three fiducials are given in both frames.
    % For example, Head lists the fiducial coordinates (taken from the
    % electrode data) in the head frame, while Head' are the fiducial
    % coordinates from the image data expressed in the head frame.
    % Similarly, Image lists the fiducial coordinates from the image
    % data in the image frame while Image' lists those from the electrode
    % data in the image frame. The two sets of numbers should be close but
    % not identical.

    if strmatch('RPA',text,'exact'),
        format = '%7c %f %f %f';
        % Read the Right Preauricular coordinates
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.RPA.Hh = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.RPA.Mh = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.RPA.Mm = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.RPA.Hm = tmp(8:10);
    end

    if strmatch('LPA',text,'exact'),
        format = '%7c %f %f %f';
        % Read the Left Preauricular coordinates
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.LPA.Hh = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.LPA.Mh = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.LPA.Mm = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.LPA.Hm = tmp(8:10);
    end

    if strmatch('Nasion',text,'exact'),
        format = '%7c %f %f %f';
        % Read the Nasion coordinates
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.NAS.Hh = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.NAS.Mh = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.NAS.Mm = tmp(8:10);
        text = fgetl(fid);
        tmp = sscanf(text,format)';
        reg.NAS.Hm = tmp(8:10);
    end

    % The Electrode Positions block lists the coordinates (x, y, and z)
    % first in the head frame and then in the image frame.
    if strmatch('Electrode Positions',text),
        reg.Helec = zeros(1,3);
        reg.Melec = zeros(1,3);
        n = 1;
        while n < 400,
            % Read the Head space coordinates
            text = fgetl(fid);
            if isempty(text), break; end
            tmp = sscanf(text,'%f : %f %f')';
            reg.Helec(n,1:3) = tmp(2:4);
            % Read the MRI space coordinates
            text = fgetl(fid);
            tmp = sscanf(text,'%s %f %f %f')';
            reg.Melec(n,1:3) = tmp(2:4);
            n = n + 1;
        end
    end

end

% Create essential fiducial marker matrices
% The order of these points in the matrices is very
% important if they are used for coregistration
reg.fiducials.head = [ reg.NAS.Hh; reg.RPA.Hh; reg.LPA.Hh ];
reg.fiducials.mri  = [ reg.NAS.Mm; reg.RPA.Mm; reg.LPA.Mm ];


return
