function [summary_clusters, summary_peaks]=stat_summary(input_file, fwhm, ...
   df, mask_file, mask_thresh, input_thresh, flip, nconj, nvar);

%STAT_SUMMARY produces SPM-style summary analyses of T or F statistic images
%
% [SUMMARY_CLUSTERS SUMMARY_PEAKS] = STAT_SUMMARY( INPUT_FILE [, FWHM [, DF [,
% MASK_FILE [, MASK_THRESH [, INPUT_THRESH [, FLIP [, NCONJ [, NVAR ]]]]]]]])
%
% Produces an SPM-style glass brain and summary analysis of a T or F statistic
% image. P-values for local maxima and cluster sizes are based on non-isotropic
% random field theory if an FWHM image is provided (or a Bonferroni correction, 
% if smaller). The random field theory is based on the assumption that
% the search region is a sphere (in isotropic space), which is a very tight 
% lower bound for any non-spherical region, unless you supply all the
% resels. It also produces a volume of clusters labelled by their index
% (as in the printout below) in INPUT_FILE with '_cluster' before the
% extension, handy for identifying the clusters in 'register'.
%
% INPUT_FILE: Finds its local maxima and clusters above INPUT_THRESH.
% Clusters are voxels connected in any of the 2*D directions. A local maximum 
% is a voxel which is greater than or equal to all its 2*D neighbours,
% and strictly greater than at least one of them. If empty, just prints out
% average FWHM (see below).
%
% FWHM is the fwhm of a smoothing kernel applied to the data, either as a
% fwhm file from fmrilm, multistat or glim_image, or as a scalar. Default
% is 0.0, i.e. no smoothing, which is roughly the case for raw fMRI data.
% For motion corrected fMRI data, use at least 6mm;
% for PET data, this would be the scanner fwhm of about 6mm.
% If FWHM is a vector, these are treated as resels of the mask.
%
% DF=[DF1 DF2; DFW1 DFW2] is a 2 x 2 matrix of degrees of freedom.
% If DF2 is 0, then DF1 is the df of the T statistic image.
% If DF1=Inf then it calculates thresholds for the Gaussian image. 
% If DF2>0 then DF1 and DF2 are the df's of the F statistic image.
% DFW1 and DFW2 are the numerator and denominator df of the FWHM image. 
% If DF=[DF1 DF2] (row) then DFW1=DFW2=Inf, i.e. FWHM is fixed.
% If DF=[DF1; DFW1] (column) then DF2=0 and DFW2=DFW1, i.e. t statistic. 
% If DF=DF1 (scalar) then DF2=0 and DFW1=DFW2=Inf, i.e. t stat, fixed FWHM.
% If any component of DF >= 1000 then it is set to Inf. Default is Inf. 
% If FWHM is estimated by FMRILM, set DFW1=DFW2=DF outputted by FMRILM;
% if FWHM is estimated by MULTISAT, set DFW1=DF_RESID, DFW2=DF outputted 
% by MULTISTAT.
%
% MASK_FILE is a mask file. If empty, it is ignored. Default is [].
%
% MASK_THRESH defines the search volume as the first frame of MASK_FILE 
% > MASK_THRESH. If MASK_THRESH is a vector [a b], a<=b, then mask 
% is a < MASK_FILE <= b. If empty (default), calls fmri_mask_thresh. 
%
% INPUT_THRESH: If <= 1 then the second element is taken as a probability and 
% the threshold is chosen so that the uncorrected P-value is this probability.
% If INPUT_THRESH is a scalar, the second element is set equal to the first.
% The default is 0.001, i.e. the threshold satisfies P=0.001 (uncorrected).
%
% FLIP: INPUT_FILE is multiplied by FLIP before processing. For T statistic
% images, FLIP = -1 will look for negative peaks and clusters. Default is 1.
%
% NCONJ is the number of conjunctions. If NCONJ > 1, calculates P-values
% for peaks (but not clusters) of the minimum of NCONJ independent 
% SPM's - see Friston, K.J., Holmes, A.P., Price, C.J., Buchel, C.,
% Worsley, K.J. (1999). Multi-subject fMRI studies and conjunction analyses.
% NeuroImage, 10:385-396. Default is NCONJ = 1 (no conjunctions). 
%
% NVAR is the number of variables for multivariate equivalents of T and F 
% statistics, found by maximizing T^2 or F over all linear combinations of 
% variables, i.e. Hotelling's T^2 for DF1=1, Roy's maximum root for DF1>1. 
% Default is 1, i.e. univariate statistics.
%
% SUMMARY_CLUSTERS is a matrix with 6 columns:
% Col 1: index of cluster, in descending order of cluster resels. 
% Col 2: volume of cluster in mm^3.
% Col 3: resels of cluster.
% Col 4: P-value of cluster extent.
% Col 5: P-value if the cluster was chosen in advance, e.g. nearest to an ROI.
%
% SUMMARY_PEAKS is a matrix with 11 columns. 
% Col 1: index of cluster. 
% Col 2: values of local maxima, sorted in descending order withihn cluster.
% Col 3: P-value of local maxima.
% Col 4: P-value if the peak was chosen in advance, e.g. nearest to an ROI.
% Col 5: Q-value or false discovery rate ~ probability that voxel is not signal.
% Cols 6-8: i,j,k coords of local maxima in voxels, starting at 0, as in 'register'.
% Cols 9-11: x,y,z coords of local maxima in world coordinates (mm).

%############################################################################
% COPYRIGHT:   Copyright 2002 K.J. Worsley, 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that this copyright
%              notice appears in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:

if nargin < 2
   fwhm=0
end
if nargin < 3
   df=Inf
end
if nargin < 4
   mask_file=[];
end
if nargin < 5
   mask_thresh=[];
end
if nargin < 6
   input_thresh=0.001
end
if nargin < 7
   flip=1
end
if nargin<8;  nconj=1;  end
if nargin<9;  nvar=1;  end

if length(input_thresh)==1
   input_thresh=[input_thresh input_thresh];
end
if input_thresh(1) > 1
   input_thresh=input_thresh(1)
else
   x=stat_threshold(1,1,0,df,input_thresh,0.001,0.05,nconj,nvar);
   input_thresh=x(2)
end

if ~isempty(mask_file) & isempty(mask_thresh)
   mask_thresh=fmri_mask_thresh(mask_file);
end

if ~isempty(input_file)
   base=input_file(1:(length(input_file)-4));
   ext=input_file((length(input_file)-2):length(input_file));
   cluster_file=[base '_cluster.' ext];
   if exist(cluster_file);  delete(cluster_file); end
   lm=locmax(input_file,input_thresh,mask_file,mask_thresh,fwhm,flip,cluster_file);
   if isempty(lm)
      summary_clusters=[];
      summary_peaks=[];
      return
   end
   [search_volume, num_voxels]= ...
      glass_brain(input_file,input_thresh,mask_file,mask_thresh,flip);
   colormap(spectral);
end;

if isstr(fwhm)
   d=fmris_read_image(fwhm,0,0);
   numslices=d.dim(3);
   if ~isempty(mask_file)
      mask_thresh1=mask_thresh(1);
      if length(mask_thresh)>=2
         mask_thresh2=mask_thresh(2);
      else
         mask_thresh2=Inf;
      end
      d=fmris_read_image(mask_file,1:numslices,1);
      mask=d.data;
      mask= (mask>mask_thresh1 & mask<=mask_thresh2);
   else
      mask=ones(d.dim(1),d.dim(2),d.dim(3));
   end
   d=fmris_read_image(fwhm,0,0);
   if d.dim(4)>=2
      d=fmris_read_image(fwhm,1:numslices,2);
   else
      d=fmris_read_image(fwhm,1:numslices,1);
      d.data=abs(prod(d.vox(1:3)))./(d.data+(d.data<=0)).^3.*(d.data>0);
   end
   search_resels=sum(sum(sum(mask.*d.data)))
   num_voxels=sum(sum(sum(mask)))
   search_volume=num_voxels*abs(prod(d.vox(1:3)))
   average_fwhm=(search_volume/search_resels)^(1/3)
   if isempty(input_file)
      return
   end
   [p_peak, p_cluster, p_peak1, p_cluster1]= ...
      stat_threshold(search_resels, num_voxels, 1, ... 
      df, [10; lm(:,1)], input_thresh, [10; lm(:,7)], nconj, nvar);
elseif length(fwhm)==1
   [p_peak, p_cluster, p_peak1, p_cluster1]= ...
      stat_threshold(search_volume, num_voxels, fwhm, ... 
      df, [10; lm(:,1)], input_thresh, [10; lm(:,6)], nconj, nvar);
else
   [p_peak, p_cluster, p_peak1, p_cluster1]= ...
      stat_threshold(fwhm, num_voxels, 1, ... 
      df, [10; lm(:,1)], input_thresh, [10; lm(:,6)], nconj, nvar);   
end
p_peak=p_peak(2:length(p_peak));
p_cluster=p_cluster(2:length(p_cluster));
p_peak1=p_peak1(2:length(p_peak1));
p_cluster1=p_cluster1(2:length(p_cluster1));

q_value = fdr_threshold( input_file, input_thresh, ... 
   mask_file, mask_thresh, df, lm(:,1), flip, nconj, nvar);

if isnan(p_cluster(1))
   summary_clusters=unique(lm(:,5:7),'rows');
   n=size(summary_clusters,1);
   ['clus    vol  resel']
   [repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,1))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,2))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,3)*100)/100) ]
else
   summary_clusters=unique([lm(:,5:7) p_cluster p_cluster1],'rows');
   n=size(summary_clusters,1);
   ['clus    vol  resel  p-val  (one)']
   [repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,1))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,2))) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,3)*100)/100) repmat('  ',n,1) ...
         num2str(round(summary_clusters(:,4)*1000)/1000) repmat(' (',n,1)  ... 
         num2str(round(summary_clusters(:,5)*1000)/1000) repmat(')  ',n,1)]
   n_clus=sum(summary_clusters(:,4)<=0.05);
   if n_clus>0
      p=get(gca,'Position');
      axes('position',[p(1)+p(3)*7/6+0.05 p(2) 0.9-p(1)-p(3)*7/6 p(4)]);
      blob_brain(input_file,input_thresh*flip,cluster_file,[0.5 n_clus+0.5]);
      title('Clusters, P<0.05, index=');
   end
end

ext=input_file((length(input_file)-2):length(input_file));
isanalyze= all(ext=='img')
if isanalyze
   d=fmris_read_image(input_file,0,0);
   coord=lm(:,2:4).*(ones(size(lm,1),1)*d.vox)+ones(size(lm,1),1)*d.origin;
else
   h=openimage(input_file);
   coord=voxeltoworld(h,lm(:,2:4)','xyzorder zerobase noflip')';
   closeimage(h);
end
summary_peaks=[lm(:,5) lm(:,1) p_peak p_peak1 q_value' lm(:,2:4) coord ];
summary_peaks=flipud(sortrows(summary_peaks,2));
n=size(lm,1);
['clus   peak   p-val  (one)   q-val   (i   j   k)  (   x      y      z )']
[repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,1))) repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,2)*100)/100) repmat('  ',n,1)  ...  
   num2str(round(summary_peaks(:,3)*1000)/1000) repmat(' (',n,1)  ...  
   num2str(round(summary_peaks(:,4)*1000)/1000) repmat(')  ',n,1)  ... 
   num2str(round(summary_peaks(:,5)*1000)/1000) repmat('  (',n,1)  ... 
   num2str(round(summary_peaks(:,6:8))) repmat(')  (',n,1)  ... 
   num2str(round(summary_peaks(:,9)*10)/10)  repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,10)*10)/10) repmat('  ',n,1) ...
   num2str(round(summary_peaks(:,11)*10)/10) repmat(')',n,1)]

return


