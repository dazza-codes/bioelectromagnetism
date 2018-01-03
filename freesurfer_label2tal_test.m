
file = which('freesurfer_read_label');
if file,
    [path, file, ext] = fileparts(file);
    cd(path)
else
    error('cannot find "freesurfer_read_label" on the matlab path');
    return
end

% I read a subject rh.pial into tksurfer and selected a small patch of
% vertices to save them into a small .label file.  These two files were
% copied into:
% freesurfer_label2tal_test.pial
% freesurfer_label2tal_test.label

surfFile = 'freesurfer_label2tal_test.pial';
pial = freesurfer_read_surf(surfFile);
surfFile = 'freesurfer_label2tal_test.orig';
orig = freesurfer_read_surf(surfFile);

lname = 'freesurfer_label2tal_test';
label = freesurfer_read_label([], lname);
LABELi = label(:,1);
LABELv = label(:,2:4);

xfmFile = 'freesurfer_label2tal_test.talairach.xfm';
TalairachXFM = freesurfer_read_talxfm(xfmFile);

[MNIv,TALv] = freesurfer_label2tal(LABELv,TalairachXFM);

% These are the values noted from csurf/tksurfer for three vertices in the
% .label area:
% vertex   Vertex RAS          Vertex Talairach     MNI Tal              Curvature
%  79497   58.81 -4.05 16.10   69.44 -14.69 -1.09   70.14 -15.07 -2.17   -0.175448
%  81708   58.42 -2.60 15.52   69.04 -13.30 -1.81   69.74 -13.60 -2.94   -0.169868
%  82883   58.10 -1.81 15.14   68.72 -12.55 -2.26   69.41 -12.80 -3.43   -0.145019
% These vertex indices are from tksurfer, which are 0 indexed, so we need
% to add 1 for matlab compatibility and then check the results
tksurferIndex = [79497, 81708, 82883] + 1;
checkLABELv = zeros(length(tksurferIndex),3);
checkMNIv = checkLABELv;
checkTALv = checkLABELv;
for i = 1:length(tksurferIndex),
    labelIndex = find(label(:,1) == tksurferIndex(i));
    checkLABELv(i,:) = LABELv(labelIndex,:);
    checkMNIv(i,:) = MNIv(labelIndex,:);
    checkTALv(i,:) = TALv(labelIndex,:);
end

% These are the final values from these operations:
[checkLABELv, checkTALv, checkMNIv]
% 
% ans =
% 
%    58.8150   -4.0450   16.1020   69.4361  -14.6923   -1.0908   70.1375  -15.0708   -2.1713
%    58.4180   -2.5990   15.5180   69.0438  -13.3031   -1.8091   69.7412  -13.6034   -2.9427
%    58.1010   -1.8110   15.1420   68.7161  -12.5491   -2.2592   69.4102  -12.8039   -3.4329


% However, note that .label files contain .orig values:
% >> pial(LABELi(1:5),:)
% 
% ans =
% 
%    61.9822   -5.1003   17.5164
%    61.5910   -5.5861   16.0358
%    61.9058   -3.7053   18.3539
%    61.5532   -4.1727   16.7315
%    61.4371   -3.5764   16.1880
% 
% >> orig(LABELi(1:5),:)
% 
% ans =
% 
%    59.4524   -5.5179   17.1148
%    59.3218   -5.5571   16.2510
%    59.3762   -4.5750   17.4985
%    59.1861   -4.5583   16.4948
%    58.8150   -4.0450   16.1021
% 
% >> LABELv(1:5,:)
% 
% ans =
% 
%    59.4520   -5.5180   17.1150
%    59.3220   -5.5570   16.2510
%    59.3760   -4.5750   17.4980
%    59.1860   -4.5580   16.4950
%    58.8150   -4.0450   16.1020

