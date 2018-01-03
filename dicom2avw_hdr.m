function [avw, orient] = dicom2avw_hdr(info)

% This function converts the dicom info header into an analyze7.5 header
% for the UCSF China Basin 3T GE scanner

avw = avw_hdr_make;
avw.hdr.hist.descrip = info.StudyDescription;
avw.hdr.hist.exp_date = [info.AcquisitionDate(7:8),'-',info.AcquisitionDate(5:6),'-',info.AcquisitionDate(1:4)];
avw.hdr.hist.exp_time = info.AcquisitionTime;

% to comply with HIPPA (in USA), do not allocate this
%avw.hdr.hist.patient_id = info.PatientID;

avw.hdr.dime.bitpix = info.BitsAllocated;
avw = avw_hdr_check_datatype(avw);

% Orientation information (from spm_dicom_convert.m)
%-------------------------------------------------------------------
% Axial Analyze voxel co-ordinate system (left handed):
% x increases     right to left
% y increases posterior to anterior
% z increases  inferior to superior

% DICOM patient co-ordinate system (right handed):
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior

% The bottom corner of the image is given by
% info.ImagePositionPatient
% and the orientation of the image by
% info.ImageOrientationPatient

%info.ImagePositionPatient
%(x,y,z)
%centre of top-left pixel 

%info.ImageOrientationPatient 
%(rx, ry, rz, cx, cy, cz) 
% image row and column directions

rowOrient = info.ImageOrientationPatient(1:3)';
colOrient = info.ImageOrientationPatient(4:6)';

orient = '';
if isequal(rowOrient,[1 0 0]) && isequal(colOrient,[0 1 0]),
    % This is an axial slice plane, with 
    % rows in x (left to right) and 
    % columns in y (posterior to anterior)
    orient = 'axial';
end

if isequal(rowOrient,[1 0 0]) && isequal(colOrient,[0 0 1]),
    % This is a coronal slice plane, with 
    % rows in x (left to right) and 
    % columns in z (inferior to superior)
    orient = 'coronal';
end

if isequal(rowOrient,[0 1 0]) && isequal(colOrient,[0 0 1]),
    % This is a sagittal slice plane, with 
    % rows in y (anterior to posterior) and 
    % columns in z (inferior to superior)
    orient = 'sagittal';
end

% The China Basin 3T GE scanner puts the slices / volume into this:
if isfield(info,'Private_0021_104f'),
    Nslices = info.Private_0021_104f;
else
    Nslices = 1;
end


% this depends on image orientation.
switch orient,
    case {'axial','axial-flipped'},
        avw.hdr.dime.dim(2) = info.Rows;    % x dim
        avw.hdr.dime.dim(3) = info.Columns; % y dim
        avw.hdr.dime.dim(4) = Nslices;      % z dim
        avw.hdr.dime.pixdim(2) = info.PixelSpacing(1);
        avw.hdr.dime.pixdim(3) = info.PixelSpacing(2);
        avw.hdr.dime.pixdim(4) = info.SliceThickness;
    case {'coronal','coronal-flipped'},
        avw.hdr.dime.dim(2) = info.Rows;    % x dim
        avw.hdr.dime.dim(3) = Nslices;      % y dim
        avw.hdr.dime.dim(4) = info.Columns; % z dim
        avw.hdr.dime.pixdim(2) = info.PixelSpacing(1);
        avw.hdr.dime.pixdim(3) = info.SliceThickness;
        avw.hdr.dime.pixdim(4) = info.PixelSpacing(2);
    case {'sagittal','sagittal-flipped',}
        avw.hdr.dime.dim(2) = Nslices;      % x dim
        avw.hdr.dime.dim(3) = info.Rows;    % y dim
        avw.hdr.dime.dim(4) = info.Columns; % z dim
        avw.hdr.dime.pixdim(2) = info.SliceThickness;
        avw.hdr.dime.pixdim(3) = info.PixelSpacing(1);
        avw.hdr.dime.pixdim(4) = info.PixelSpacing(2);
    otherwise,
        error('unknown image orientation');
end

avw.hdr.dime.pixdim(5) = info.RepetitionTime;


return




% 	% Orientation information (from spm_dicom_convert.m)
% 	%-------------------------------------------------------------------
% 	% Axial Analyze voxel co-ordinate system:
% 	% x increases     right to left
% 	% y increases posterior to anterior
% 	% z increases  inferior to superior
% 
% 	% DICOM patient co-ordinate system:
% 	% x increases     right to left
% 	% y increases  anterior to posterior
% 	% z increases  inferior to superior
% 
% 	% T&T co-ordinate system:
% 	% x increases      left to right
% 	% y increases posterior to anterior
% 	% z increases  inferior to superior
% 
% 	analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];
% 
% 	vox    = [info.PixelSpacing; info.SpacingBetweenSlices];
% 	pos    = info.ImagePositionPatient';
% 	orient = reshape(info.ImageOrientationPatient,[3 2]);
% 	orient(:,3) = null(orient');
% 	if det(orient)<0, orient(:,3) = -orient(:,3); end;
% 
% 	% The image position vector is not correct. In dicom this vector points to
% 	% the upper left corner of the image. Perhaps it is unlucky that this is
% 	% calculated in the syngo software from the vector pointing to the center of
% 	% the slice (keep in mind: upper left slice) with the enlarged FoV.
% 	dicom_to_patient = [orient*diag(vox) pos ; 0 0 0 1];
% 	truepos          = dicom_to_patient *[(size(mosaic)-dim(1:2))/2 0 1]';
% 	dicom_to_patient = [orient*diag(vox) truepos(1:3) ; 0 0 0 1];
% 	patient_to_tal   = diag([-1 -1 1 1]);
% 	mat              = patient_to_tal * dicom_to_patient * analyze_to_dicom;

