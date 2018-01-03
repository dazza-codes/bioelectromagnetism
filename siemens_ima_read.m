function [ ima, machine ] = siemens_ima_read(fileprefix,IMGorient)

% SIEMENS_IMA_READ - read Siemens .ima data file (*.ima)
% 
% [ ima, machine ] = siemens_ima_read(fileprefix, orient, machine)
%
% fileprefix - a string, the filename without the .ima extension
% 
% orient - force reading of specified orientation, integer values:
%
%          '', read header history orient field
%          0,  transverse/axial unflipped
%          1,  coronal unflipped
%          2,  sagittal unflipped
%          3,  transverse/axial flipped
%          4,  coronal flipped
%          5,  sagittal flipped
% 
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-be' but the routine
%           will automatically switch between little and big
%           endian.  It reports the appropriate machine 
%           format and can return the machine value.
% 
% Returned values:
% 
% ima.hdr - a struct with image data parameters.
% ima.img - a 3D matrix of image data (double precision).
% 
% See also: IMA_HDR_READ (called by this function)
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Licence:  GNU GPL, no express or implied warranties
% History:  10/2002, Darren.Weber@flinders.edu.au
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\nUnder development - not functional yet...sorry.\n\n');
return


if ~exist('fileprefix','var'),
    fprintf('SIEMENS_IMA_READ: No input fileprefix - see help siemens_ima_read\n');
    return;
end
if ~exist('IMGorient','var'), IMGorient = ''; end
if ~exist('machine','var'), machine = 'ieee-le'; end


% MAIN

fid = fopen(sprintf('%s.ima',fileprefix),'r',machine);
if fid < 0,
    msg = sprintf('Cannot open file %s.ima\n',fileprefix);
    error(msg);
else
    ima = read_image(fid,fileprefix,IMGorient,machine);
    ima.fileprefix = fileprefix;
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ima ] = read_image(fid,fileprefix,IMGorient,machine)
    
    % Read the file header
    [ ima, machine ] = ima_hdr_read(fileprefix,machine);
    
    if ~isequal(ima.hdr.hk.sizeof_hdr,6144),
        fprintf('SIEMENS_IMA_READ: Failed reading %s Analyze image\n',fileprefix);
        ima.img = [];
        fclose(fid);
        return
    end
    
    
    % read the whole .img file into matlab (faster)
    fprintf('SIEMENS_IMA_READ: Reading Siemens image data.\n');
    fseek(fid,6144,'bof');
    precision = 'double';
    tmp = fread(fid,inf,sprintf('%s=>double',precision));
    fclose(fid);
    
    % Now partition the img data into xyz
    
    
    
    PixelDim = double(ima.hdr.dime.dim(2));
    RowDim   = double(ima.hdr.dime.dim(3));
    SliceDim = double(ima.hdr.dime.dim(4));
    
    
    
    if ~isempty(IMGorient), ima.hdr.hist.orient = IMGorient; end,
    
    
    switch double(ima.hdr.hist.orient),
    
    case 0, % transverse/axial unflipped
        
        % For the 'transverse unflipped' type, the voxels are stored with
        % Pixels in 'x' axis (varies fastest) - from patient right to left
        % Rows in   'y' axis                  - from patient posterior to anterior
        % Slices in 'z' axis                  - from patient inferior to superior
        
        fprintf('SIEMENS_IMA_READ: Image orient appears to be axial unflipped\n');
        
        ima.img = zeros(PixelDim,RowDim,SliceDim);
        
        n = 1;
        for z = 1:SliceDim,
            for y = 1:RowDim,
                x = 1:PixelDim;
                ima.img(x,y,z) = tmp(n:n+(PixelDim-1));
                n = n + PixelDim;
            end
        end
        
    end
return


function [ ima, machine ] = ima_hdr_read(fileprefix,machine)



% Siemens Magnetom Vision Native Format
% 
% The exact details of the format are not known, but 
% a little guess work has determined what follows. The 
% data types are Sun, hence the byte order is big-endian 
% and the all the floats I have found are doubles. Offsets 
% here are in bytes from the start of the header. The 
% uncompressed image data starts at offset 6144.
% 
%         0         u_int      SiemensStudyDateYYYY 
%         4         u_int      SiemensStudyDateMM 
%         8         u_int      SiemensStudyDateDD 
%         12        u_int      AcquisitionDateYYYY 
%         16        u_int      AcquisitionDateMM 
%         20        u_int      AcquisitionDateDD 
%         24        u_int      ImageDateYYYY             
%         28        u_int      ImageDateMM             
%         32        u_int      ImageDateDD             
%         36        u_int      SiemensStudyTimeHH 
%         40        u_int      SiemensStudyTimeMM 
%         44        u_int      SiemensStudyTimeSS
%         52        u_int      AcquisitionTimeHH 
%         56        u_int      AcquisitionTimeMM 
%         60        u_int      AcquisitionTimeSS 
%         68        u_int      ImageTimeHH             
%         72        u_int      ImageTimeMM             
%         76        u_int      ImageTimeSS             
%         96        char[7]    Manufacturer
%         105       char[25]   InstitutionName 
%         186       char[4]    Annotation             
%         281       char[15]   ModelName 
%         412       u_int      LastMoveDateYYYY 
%         416       u_int      LastMoveDateMM 
%         420       u_int      LastMoveDateDD 
%         424       u_int      LastMoveTimeHH 
%         428       u_int      LastMoveTimeMM 
%         432       u_int      LastMoveTimeSS 
%         768       char[25]   PatientName
%         795       char[12]   PatientID
%         808       u_int      DOBYYYY             
%         812       u_int      DOBMM                   
%         816       u_int      DOBDD                   
%         851       char[3]    PatientAge             
%         854       char       PatientAgeUnits      ('Y'=years) 
%         1052      u_int      RegistrationDateYYYY 
%         1056      u_int      RegistrationDateMM 
%         1060      u_int      RegistrationDateDD 
%         1064      u_int      RegistrationTimeHH 
%         1068      u_int      RegistrationTimeMM 
%         1072      u_int      RegistrationTimeSS 
%         1544      double     SliceThickness 
%         1560      double     RepetitionTime
%         1568      double     EchoTime
%         1592      double     FrequencyMHz
%         1639      char[5]    Station             
%         1712      u_int      CalibrationDateYYYY 
%         1716      u_int      CalibrationDateMM 
%         1720      u_int      CalibrationDateDD 
%         1724      u_int      CalibrationTimeHH 
%         1728      u_int      CalibrationTimeMM 
%         1732      u_int      CalibrationTimeSS 
%         1767      char[16]   ReceivingCoil
%         1828      char[4]    ImagedNucleus
%         2112      double     FlipAngle
%         2560      double     MagneticFieldStrength
%         2864      u_int      DisplayMatrixSize 
%         2944      char[65]   SequencePrgName
%         3009      char[65]   SequenceWkcName 
%         3074      char[9]    SequenceAuthor 
%         3083      char[8]    SequenceType
%         3744      double     FOVRow             
%         3752      double     FOVColumn             
%         3768      double     CenterPointX             
%         3776      double     CenterPointY             
%         3784      double     CenterPointZ             
%         3792      double     NormalVectorX             
%         3800      double     NormalVectorY             
%         3808      double     NormalVectorZ             
%         3816      double     DistanceFromIsocenter 
%         3832      double     RowVectorX             
%         3840      double     RowVectorY             
%         3848      double     RowVectorZ             
%         3856      double     ColumnVectorX             
%         3864      double     ColumnVectorY             
%         3872      double     ColumnVectorZ             
%         3880      char[3]    OrientationSet1Top 
%         3884      char[3]    OrientationSet1Left 
%         3888      char[3]    OrientationSet1Back 
%         3892      char[3]    OrientationSet2Down 
%         3896      char[3]    OrientationSet2Right 
%         3900      char[3]    OrientationSet2Front 
%         3904      char[32]   SequenceName             
% 
%         4136      double     GapRatio (see below)
% 
%         5000      double     PixelSizeRow             
%         5008      double     PixelSizeColumn 
% 
%         5504      char[12]   TextPatientID             
%         5517      char       TextPatientSex 
%         5518      char[3]    TextPatientAge 
%         5521      char       TextPatientAgeUnits       ('Y'=years) 
%         5529      char[7]    TextPatientPosition 
%         5541      char[5]    TextImageNumberFlag       ('IMAGE'=image) 
%         5546      char[3]    TextImageNumber
%         5559      char[2]    TextDateDD             
%         5562      char[3]    TextDateMM             
%         5566      char[4]    TextDateYYYY             
%         5571      char[2]    TextTimeHH             
%         5574      char[2]    TextTimeMM             
%         5577      char[2]    TextAcquisitionTimeFlag   ('TA'=acquisition time) 
%         5583      char[2]    TextAcquisitionTimeMM             
%         5586      char[2]    TextAcquisitionTimeSS             
%         5601      char[4]    TextAnnotation
%         5655      char[25]   TextOrganization             
%         5682      char[5]    TextStation                   
%         5695      char[3]    TextAcquisitionMatrixPhase
%         5698      char       TextAcquisitionMatrixPhaseAxis  ('h'=horizontal,' '=vertical) 
%         5700      char[3]    TextAcquisitionMatrixFreq 
%         5703      char       TextAcquisitionMatrixFreqO      ('o'=o,' '=blank) 
%         5704      char       TextAcquisitionMatrixFreqS      ('s'=s,' '=blank) 
%         5706      char[8]    TextSequence                   
%         5714      char[3]    TextFlipAngle                   
%         5718      char[4]    TextScanNumberFlag        ('SCAN'=scan) 
%         5723      char[3]    TextScanNumberA             
%         5726      char[3]    TextScanNumberB             
%         5730      char[2]    TextRepetitionTimeFlag    ('TR'=tr) 
%         5734      char[7]    TextRepetitionTime             
%         5742      char[2]    TextEchoTimeFlag          ('TE'=te) 
%         5746      char[5]    TextEchoTime                   
%         5752      char       TextEchoNumber             
%         5790      char[2]    TextSliceThicknessFlag    ('SL'=slice thickness) 
%         5794      char[7]    TextSliceThickness
%         5802      char[2]    TextSlicePositionFlag     ('SP'=slice position) 
%         5806      char[7]    TextSlicePosition
%         5814      char[3]    TextAngleFlag1            ('Sag'=sagittal,'Cor'=coronal,'Tra'=transverse) 
%         5817      char       TextAngleFlag2            ('>'=gt,'<'=lt) 
%         5818      char[3]    TextAngleFlag3            ('Sag'=sagittal,'Cor'=coronal,'Tra'=transverse) 
%         5821      char[4]    TextAngle                   
%         5838      char[3]    TextFOVFlag               ('FoV'=field of view) 
%         5842      char[3]    TextFOVH                   
%         5846      char[3]    TextFOVV                   
%         5874      char[2]    TextTablePositionFlag     ('TP'=table position) 
%         5878      char[7]    TextTablePosition             
%         5938      char[5]    TextStudyNumberFlag       ('STUDY'=study) 
%         5943      char[2]    TextStudyNumber
%         5956      char[2]    TextDOBDD                   
%         5959      char[3]    TextDOBMM                   
%         5963      char[4]    TextDOBYYYY                   
%         5992      char[3]    TextStudyNumberFlag2      ('STU'=study) 
%         5996      char[3]    TextImageNumberFlag2      ('IMA'=study) 
%         5999      char[2]    TextStudyNumber2             
%         6002      char[2]    TextImageNumber2             
%         6013      char[5]    TextStudyImageNumber3             
%         6031      char[15]   TextModelName                   
%         6058      char[25]   TextPatientName             
%         6085      char[2]    TextScanStartTimeHH             
%         6088      char[2]    TextScanStartTimeMM             
%         6091      char[2]    TextScanStartTimeSS        


% I have asked David Clunie to update his FAQ to include information for
% finding size of the gap given a Siemens Vision format image. For the record
% the answer is:  GapRatio is at byte 4136 (measured from the beginning of
% the file) and it is of type double (8-byte float, Big Endian). For example,
% if SliceThickness is 8mm and the GapRatio is 0.25, then the Gap is
% (8x0.25=) 2mm, and the interslice distance is (8+2=) 10mm.


return
