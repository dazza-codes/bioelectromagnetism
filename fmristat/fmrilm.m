function [df, p, spatial_av]= ...
   fmrilm(input_file, output_file_base, X_cache, ...
   contrast, exclude, which_stats, fwhm_rho, n_trends, confounds, ...
   contrast_is_delay, num_hrf_bases, basis_type, numlags, df_limit)

%FMRILM fits a linear model to fMRI time series data.
%
% The method is based on linear models with correlated AR(p) errors:
% Y = hrf*X b + e, e_t=a_1 e_(t-1) + ... + a_p e_(t-p) + white noise_t. 
% 
% [DF,P,SPATIAL_AV] = 
%    FMRILM( INPUT_FILE, OUTPUT_FILE_BASE, X_CACHE, CONTRAST  
%    [, EXCLUDE [, WHICH_STATS [, FWHM_RHO [, N_TRENDS [, CONFOUNDS
%    [, CONTRAST_IS_DELAY [,NUM_HRF_BASES [, BASIS_TYPE [,NUMLAGS 
%    [, DF_LIMIT ]]]]]]]]]] )
% 
% INPUT_FILE (Y) is the name of a single 4D fMRI image file with multiple  
% frames, or a matrix of image file names, each with a single 3D frame,
% either ANALYZE (.img) or MINC (.mnc) format. Extra blanks are ignored. 
% File separator can be / or \ on Windows. Gzipped files are gunzipped on     
% unix. Frames must be equally spaced in time, unless FWHM_RHO=Inf.
% 
% OUTPUT_FILE_BASE: matrix whose rows are the base for output statistics,
% one base for each row of CONTRAST, padded with extra blanks if 
% necessary; they will be removed before use.
%
% X_CACHE: A structure usually supplied by FMRIDESIGN.
% X_CACHE.TR: TR, average time between frames (secs), only used for 
% calculating the number of temporal drift terms (see N_TEMPORAL below).  
% X_CACHE.X: A cache of the design matrices (hrf*X) stored as a 4D array. 
% Dim 1: frames; Dim 2: response variables; Dim 3: 4 values, corresponding 
% to the stimuli convolved with: hrf, derivative of hrf, first and second 
% spectral basis functions over the range in SHIFT; Dim 4: slices. 
% X_CACHE.W: A 3D array of coefficients of the basis functions in X_CACHE.X.
% Dim 1: frames; Dim 2: response variables; Dim 3: 5 values: 
% coefficients of the hrf and its derivative, coefficients of the first and 
% second spectral basis functions, shift values.
%
% CONTRAST is a matrix whose rows are contrasts for the
% response variables, and columns are the contrast. Extra columns can 
% be added to estimate contrasts in the temporal and spatial trends, and 
% the confounds (in that order - see N_TRENDS and CONFOUNDS below). 
% 
% EXCLUDE is a list of frames that should be excluded from the
% analysis. This must be used with Siemens EPI scans to remove the
% first few frames, which do not represent steady-state images.
% If NUMLAGS=1, the excluded frames can be arbitrary, otherwise they 
% should be from the beginning and/or end. Default is [1].
% 
% WHICH_STATS: logical matrix indicating output statistics by 1.
% If WHICH_STATS to a row vector of length 9, then it is used for all
% contrasts; if the number of columns is less than 9 it is padded
% with zeros. If empty (default), only SPATIAL_AV is returned.
% Rows correspond to rows of CONTRAST, columns correspond to:
%   1: T statistic image, OUTPUT_FILE_BASE_mag_t.mnc or .img (for contrasts 
%      in the delays, mag is replaced by del in the extension).  
%      The degrees of freedom is DF. Note that tstat=effect/sdeffect. 
%      To avoid large values in 'register', T is truncated to 100.
%   2: effect (b) image, OUTPUT_FILE_BASE_mag_ef.mnc or .img.
%   3: standard deviation of the effect, OUTPUT_FILE_BASE_mag_sd.mnc or .img. 
%   4: F-statistic image, OUTPUT_FILE_BASE_mag_F.mnc or .img, for testing
%      all rows of CONTRAST that deal with magnitudes (see CONTRAST_IS_DELAY
%      below). The degrees of freedom are [DF, P].
%      To avoid large values in 'register', F is truncated to 10000.
%   5: the temporal autocorrelation(s), OUTPUT_FILE_BASE_rho.mnc or .img.
%   6: the residuals from the model, OUTPUT_FILE_BASE_resid.mnc or .img,
%      only for the non-excluded frames (warning: uses lots of memory). 
%   7: the whitened residuals from the model normalized by dividing
%      by their root sum of squares, OUTPUT_FILE_BASE_wresid.mnc or .img, 
%      only for the non-excluded frames (warning: uses lots of memory).
%   8: the AR parameter(s) a_1 ... a_p, in OUTPUT_FILE_BASE_A.mnc or .img.
%   9: the estimated FWHM in the first frame of OUTPUT_FILE_BASE_FWHM.mnc 
%      or .img, and the resels in second frame - uses more time and memory.
% Estimates the effective FWHM in mm of the residuals from the linear model,
% as if the residuals were white noise smoothed with a Gaussian filter 
% whose fwhm was FWHM. Applies a first-order correction for small degrees 
% of freedom and small FWHM relative to the voxel size. Great care has been
% taken to make sure that FWHM is unbiased, particularly for low df, so that
% if it is smoothed spatially then it remains unbiased. The bias in FWHM
% is less than 1/10 of the voxel size for FWHM > 2*voxel size. However the 
% continuity correction for small FWHM breaks down if FWHM < voxel size,
% in which case FWHM is generally too large. If FWHM > 50, FWHM = 50.
%   
% FWHM_RHO is the fwhm in mm of a 3D Gaussian kernel needed to smooth the
% autocorrelation of residuals with 100 df. This is NOT the fwhm of the  
% image data as used by stat_threshold, but should be fairly large. The 
% actual fwhm used is FWHM_RHO*(100/DF)^(1/3), where DF is the residual 
% df of the linear model (outputed - see below). Setting FWHM_RHO to Inf
% smooths the autocorrelation to 0, i.e. it assumes the frames
% are uncorrelated (useful for TR>10 seconds). Default is 10. 
% If FWHM_RHO is a character string, then this file is used for the
% autocorrelations, e.g. the _rho.mnc or .img file created by a previous 
% run - this saves execution time.
% 
% N_TRENDS=[N_TEMPORAL N_SPATIAL PCNT] is the number of trends to remove:
% N_TEMPORAL: number of cubic spline temporal trends to be removed per 6 
%  minutes of scanner time (so it is backwards compatible). Temporal  
%  trends are modeled by cubic splines, so for a 6 minute run, N_TEMPORAL
%  <=3 will model a polynomial trend of degree N_TEMPORAL in frame times, 
%  and N_TEMPORAL>3 will add (N_TEMPORAL-3) equally spaced knots.
%  N_TEMPORAL=0 will model just the constant level and no temporal trends.
%  N_TEMPORAL=-1 will not remove anything, in which case the design matrix 
%  is completely determined by X_CACHE.X. 
% N_SPATIAL: order of the polynomial in the spatial average (SPATIAL_AV)  
%  weighted by first non-excluded frame; 0 will remove no spatial trends. 
%  [-1 0 0] will not remove anything, in which case the design matrix is
%  completely determined by X_CACHE.X and CONFOUNDS. 
% PCNT: if PCNT=1, then the data is converted to percentages before
%  analysis by dividing each frame by its spatial average, * 100%.
% Default is [3 1 1], i.e. 3 temporal trends per 6 minutes, 
% linear spatial trend, and conversion of data to percent of whole brain.
% For backwards compatibility N_TRENDS padded with defaults if too short.
% If you've previously calculated SPATIAL_AV, add it to the end of N_TRENDS. 
%
% CONFOUNDS: A matrix or array of extra columns for the design matrix
% that are not convolved with the HRF, e.g. movement artifacts. 
% If a matrix, the same columns are used for every slice; if an array,
% the first two dimensions are the matrix, the third is the slice.
% For functional connectivity with a single voxel, use fmri_interp
% to resample the reference data at different slice times. 
% Default is [], i.e. no confounds.
% 
% CONTRAST_IS_DELAY is a logical vector indicating if the row of CONTRAST
% is a contrast in the delays (1) or magnitudes (0). Delays are shifts
% of the time origin of the HRF, measured in seconds. Note that you cannot
% estimate delays of the polynomial terms or confounds. Note that
% F statistics are not yet available with this option. If the length
% of CONTRAST_IS_DELAY is less then the number of contrasts, it is padded
% with zeros. Default is 0, i.e. all contrasts refer to magnitudes.
%
% NUM_HRF_BASES is a row vector indicating the number of basis functions
% for the hrf for each response, either 1 or 2 at the moment. At least  
% one basis functions is needed to estimate the magnitude, but two basis 
% functions are needed to estimate the delay. If empty (default), then 
% NUM_HRF_BASES is 2 for each response where CONTRAST is non-zero and  
% CONTRAST_IS_DELAY = 1, otherwise NUM_HRF_BASES = 1. By setting 
% NUM_HRF_BASES = 2 you can allow for an unknown delay, without  
% actually estimating it.  Example:   
%     CONTRAST=[1 -1 0 0; 1 0 -1 0];    CONTRAST_IS_DELAY=[0 1];  
% The first contrast is the magnitude of response_1 - response_2;
% the second contrast is the delay of response_1 - response_3. 
% The default setting of NUM_HRF_BASES is [2 1 2 1]. By setting        
% NUM_HRF_BASES=[2 2 2 1] unknown delays are allowed for the second 
% response but not actually estimated.
%
% BASIS_TYPE selects the basis functions for the hrf used for delay
% estimation, or whenever NUM_HRF_BASES = 2. These are convolved with
% the stimulus to give the responses in Dim 3 of X_CACHE.X:
% 'taylor' - use hrf and its first derivative (components 1 and 2), or 
% 'spectral' - use first two spectral bases (components 3 and 4 of Dim 3).
% Ignored if NUM_HRF_BASES = 1, in which case it always uses component 1,  
% i.e. the hrf is convolved with the stimulus. Default is 'spectral'. 
%
% NUMLAGS is the order (p) of the autoregressive model. Default is 1.
%
% DF_LIMIT controls which method is used for estimating FWHM. If DF > 
% DF_LIMIT, then the FWHM is calculated assuming the Gaussian filter is 
% arbitrary. However if DF is small, this gives inaccurate results, so
% if DF <= DF_LIMIT, the FWHM is calculated assuming that the axes of
% the Gaussian filter are aligned with the x, y and z axes of the data. 
% Default is 4. 
%
% DF is the (residual) degrees of freedom of the statistics.
%
% P is the numerator degrees of freedom of the F statistic.
%
% SPATIAL_AV is the column vector of the spatial average (SPATIAL_AV) 
% of the frames weighted by the first non-excluded frame.

%############################################################################
% COPYRIGHT:   Copyright 2002 K.J. Worsley
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca, liao@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that the above copyright
%              notice appear in all copies.  The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:

if nargin < 5
   exclude=[1]
end
if nargin < 6
   which_stats=[]
end
if nargin < 7
   fwhm_rho=10
end
if nargin < 8
   n_trends=[3 1 1];
end
if nargin < 9
   confounds=[]
end
if nargin < 10
   contrast_is_delay=[]
end
if nargin < 11
   num_hrf_bases=[];
end
if nargin < 12
   basis_type='spectral'
end
if nargin < 13
   numlags=1
end
if nargin < 14
   df_limit=4
end

if isempty(contrast)
   contrast=[1]
end
if isempty(contrast_is_delay)
   contrast_is_delay=[0]
end
if isempty(fwhm_rho)
   fwhm_rho=10
end
if isempty(n_trends)
   n_trends=[3 1 1]
end
n_temporal=round(n_trends(1))
if length(n_trends)>=2
   n_spatial=n_trends(2)
else
   n_spatial=1
end
if length(n_trends)>=3
   ispcnt=(n_trends(3)==1)
else
   ispcnt=1
end

if isempty(basis_type)
   basis_type='spectral'
end   
switch lower(basis_type)
case 'taylor',    
   basis1=1;
   basis2=2;
case 'spectral',    
   basis1=3;
   basis2=4;
otherwise, 
   disp('Unknown basis_type.'); 
   return
end

% Open image:

numfiles=size(input_file,1)
d=fmris_read_image(deblank(input_file(1,:)),0,0);
d.dim
numframes=max(d.dim(4),numfiles);
numslices=d.dim(3);
numys=d.dim(2);
numxs=d.dim(1);
numpix=numxs*numys;
Steps=d.vox

% Keep time points that are not excluded:

allpts = 1:numframes;
allpts(exclude) = zeros(1,length(exclude));
keep = allpts( find( allpts ) );
n=length(keep)

% Create spatial average, weighted by the first frame:

if (ispcnt | n_spatial>=1)
   if length(n_trends)<=3
      spatial_av=zeros(numframes,1);
      if numfiles==1
         d1=fmris_read_image(input_file,1:numslices,keep(1));
      else
         d1=fmris_read_image(deblank(input_file(keep(1),:)),1:numslices,1);
      end
      tot=sum(sum(sum(d1.data)));
      for i=1:numframes
         if numfiles==1
            d=fmris_read_image(input_file,1:numslices,i);
         else
            d=fmris_read_image(deblank(input_file(i,:)),1:numslices,1);
         end
         spatial_av(i)=sum(sum(sum(d.data.*d1.data)))/tot;
      end
      clear d1;
   else
      spatial_av=n_trends((1:numframes)+3)';
   end
end
   
if isempty(which_stats)
   df=[]
   p=[]
   return
end

% Create temporal trends:

n_spline=round(n_temporal*X_cache.TR*n/360)
if n_spline>=0 
   trend=((2*keep-(max(keep)+min(keep)))./(max(keep)-min(keep)))';
   if n_temporal<=3
      temporal_trend=(trend*ones(1,n_spline+1)).^(ones(n,1)*(0:n_spline));
   else
      temporal_trend=(trend*ones(1,4)).^(ones(n,1)*(0:3));
      knot=(1:(n_spline-3))/(n_spline-2)*(max(keep)-min(keep))+min(keep);
      for k=1:length(knot)
         cut=keep'-knot(k);
         temporal_trend=[temporal_trend (cut>0).*(cut./max(cut)).^3];
      end
   end
else
   temporal_trend=[];
end 

% Create spatial trends:

if n_spatial>=1 
   trend=spatial_av(keep)-mean(spatial_av(keep));
   spatial_trend=(trend*ones(1,n_spatial)).^(ones(n,1)*(1:n_spatial));
else
   spatial_trend=[];
end 

trend=[temporal_trend spatial_trend];

% Add confounds:

numtrends=size(trend,2)+size(confounds,2)
Trend=zeros(n,numtrends,numslices);
for slice=1:numslices
   if isempty(confounds)
      Trend(:,:,slice)=trend;
   else  
      if length(size(confounds))==2
         Trend(:,:,slice)=[trend confounds(keep,:)];
      else
         Trend(:,:,slice)=[trend confounds(keep,:,slice)];
      end
   end
end

if ~isempty(X_cache)
   numresponses=size(X_cache.X,2)
else
   numresponses=0
end

% Make full contrasts:

numcontrasts=size(contrast,1)
contrast_is_delay=[contrast_is_delay ...
      zeros(1,numcontrasts-length(contrast_is_delay))]
if isempty(num_hrf_bases)
   num_hrf_bases=ones(1,numresponses);
end
for k=find(contrast_is_delay)
   num_hrf_bases(find(contrast(k,:)~=0))=2;
end
num_hrf_bases
contrasts=[contrast zeros(numcontrasts,numresponses+numtrends-size(contrast,2))]
p=rank(contrast(~contrast_is_delay,:))

% Check for estimability:

tolerance=0.0000001;
for slice=1:numslices
   if ~isempty(X_cache)
      X=[squeeze(X_cache.X(keep,:,1,slice)) Trend(:,:,slice)];
   else
      X=Trend(:,:,slice);
   end
   dfs(slice)=n-rank(X);
   NullSpaceX=null(X);
   Cmhalf=diag(1./sqrt(diag(contrasts*contrasts')));
   ContrastNullSpaceX=Cmhalf*contrasts*NullSpaceX;
   nonest=sum(abs(ContrastNullSpaceX)>tolerance,2);
   if sum(nonest)>0
      fprintf(['Error: the following contrasts are nonestimable in slice ' ...
            num2str(slice) ':']);
      RowsNonEstContrast=find(nonest>0)
      NonEstContrast=contrasts(RowsNonEstContrast,:)
      NullSpaceX
      ContrastNullSpaceX
      return
   end
end

% Calculate df (needed for the fwhm for smoothing the AR coeffs):

dfs=zeros(1,numslices);
for slice=1:numslices
   if ~isempty(X_cache)
      X=[squeeze(X_cache.X(keep,num_hrf_bases==1,1,slice)) ... 
            squeeze(X_cache.X(keep,num_hrf_bases==2,basis1,slice)) ...
            squeeze(X_cache.X(keep,num_hrf_bases==2,basis2,slice)) ...
            Trend(:,:,slice)];
   else
      X=Trend(:,:,slice);
   end
   dfs(slice)=n-rank(X);
end
dfs
df=round(mean(dfs))

% Open files for output:

if size(which_stats,1)==1
   which_stats=repmat(which_stats,numcontrasts,1);
end
which_stats=[which_stats zeros(numcontrasts,9-size(which_stats,2))]

parent_file=deblank(input_file(1,:));
ext=parent_file(min(findstr('.',parent_file))+(1:3));

typestat=['_t. '; '_ef.'; '_sd.'];
for i=1:numcontrasts
   if contrast_is_delay(i)
      param='_del';
   else
      param='_mag';
   end
   for stat=1:3
      if which_stats(i,stat) 
         out(i,stat).file_name=[deblank(output_file_base(i,:)) param deblank(typestat(stat,:)) ext];
         out(i,stat).dim=[numxs numys numslices 0];
         out(i,stat).parent_file=parent_file;
      end
   end
end

if which_stats(1,4)
   out_Fstat.file_name=[deblank(output_file_base(1,:)) '_mag_F.' ext];
   out_Fstat.dim=[numxs numys numslices 0];
   out_Fstat.parent_file=parent_file;
end
if which_stats(1,6)
   out_resid.file_name=[deblank(output_file_base(1,:)) '_resid.' ext];
   out_resid.dim=[numxs numys numslices n];
   out_resid.parent_file=parent_file;
end
if which_stats(1,7)
   out_wresid.file_name=[deblank(output_file_base(1,:)) '_wresid.' ext];
   out_wresid.dim=[numxs numys numslices n];
   out_wresid.parent_file=parent_file;
end
if which_stats(1,8)
   out_A.file_name=[deblank(output_file_base(1,:)) '_A.' ext];
   out_A.dim=[numxs numys numslices numlags];
   out_A.parent_file=parent_file;
end
if which_stats(1,9)
   out_fwhm.file_name=[deblank(output_file_base(1,:)) '_fwhm.' ext];
   out_fwhm.dim=[numxs numys numslices 2];
   out_fwhm.parent_file=parent_file;
end

% Setup for finding rho:

rho_slice=zeros(numxs,numys);
rho_vol=squeeze(zeros(numpix, numslices, numlags));

indk1=((keep(2:n)-keep(1:n-1))==1);
k1=find(indk1)+1;

% If fwhm_rho is a file name, then read in rho, otherwise estimte it:

if isnumeric(fwhm_rho) & fwhm_rho < Inf 
   
   if which_stats(1,5)
      out_rho.file_name=[deblank(output_file_base(1,:)) '_rho.' ext];
      out_rho.dim=[numxs numys numslices numlags-(numlags==1)];
      out_rho.parent_file=parent_file;
   end
   
   % Smoothing rho in slice is done using conv2 with a kernel ker_x and ker_y.
   
   if fwhm_rho>0
      fwhm_x=fwhm_rho/abs(Steps(1))*(100/df)^(1/3);
      ker_x=exp(-(-ceil(fwhm_x):ceil(fwhm_x)).^2*4*log(2)/fwhm_x^2);
      ker_x=ker_x/sum(ker_x);
      fwhm_y=fwhm_rho/abs(Steps(2))*(100/df)^(1/3);
      ker_y=exp(-(-ceil(fwhm_y):ceil(fwhm_y)).^2*4*log(2)/fwhm_y^2);
      ker_y=ker_y/sum(ker_y);
   else
      ker_x=1;
      ker_y=1;
   end
   
   % First loop over slices, then pixels, to get the AR parameter:
   Diag1=diag(indk1,1)+diag(indk1,-1);
   Y=zeros(n,numpix);

   for slice=1:numslices
      First_pass_slice=slice
      if ~isempty(X_cache)
         X=[squeeze(X_cache.X(keep,num_hrf_bases==1,1,slice)) ... 
               squeeze(X_cache.X(keep,num_hrf_bases==2,basis1,slice)) ...
               squeeze(X_cache.X(keep,num_hrf_bases==2,basis2,slice)) ...
               Trend(:,:,slice)];
      else
         X=Trend(:,:,slice);
      end
      pinvX=pinv(X);
      
      % Preliminary calculations for unbiased estimates of autocorrelation:
      
      R=eye(n)-X*pinvX;
      if numlags==1
         M11=trace(R);
         M12=trace(R*Diag1);
         M21=M12/2;
         M22=trace(R*Diag1*R*Diag1)/2;
         M=[M11 M12; M21 M22];
      else
         M=zeros(numlags+1);
         for i=1:(numlags+1)
            for j=1:(numlags+1)
               Di=(diag(ones(1,n-i+1),i-1)+diag(ones(1,n-i+1),-i+1))/(1+(i==1));
               Dj=(diag(ones(1,n-j+1),j-1)+diag(ones(1,n-j+1),-j+1))/(1+(j==1));
               M(i,j)=trace(R*Di*R*Dj)/(1+(i>1));
            end
         end
      end
      invM=inv(M);
      %invM=eye(numlags+1); %this will undo the correction.
      
      % Read in data:
      if numfiles==1
         d=fmris_read_image(input_file,slice,keep);
         Y=reshape(d.data,numpix,n)';
      else
         for i=1:n
            d=fmris_read_image(input_file(keep(i),:),slice,1);
            Y(i,:)=reshape(d.data,1,numpix);
         end
      end
      
      % Convert to percent:
      
      if ispcnt
         for i=1:n
            Y(i,:)=Y(i,:)*(100/spatial_av(keep(i)));
         end
      end
      
      % Least squares:
      
      betahat_ls=pinvX*Y;
      resid=Y-X*betahat_ls;
      if numlags==1
         Cov0=sum(resid.*resid,1);
         Cov1=sum(resid(k1,:).*resid(k1-1,:),1);
         Covadj=invM*[Cov0; Cov1];
         rho_slice(:)=(Covadj(2,:)./ ...
            (Covadj(1,:)+(Covadj(1,:)<=0)).*(Covadj(1,:)>0))';
         rho_vol(:,slice)=reshape(conv2(ker_x,ker_y,rho_slice,'same'),numpix,1);   
      else
         for lag=0:numlags
            Cov(lag+1,:)=sum(resid(1:(n-lag),:).*resid((lag+1):n,:));
         end
         Covadj=invM*Cov;
         Coradj_slice= ( Covadj(2:(numlags+1),:) ...
            .*( ones(numlags,1)*((Covadj(1,:)>0)./ ...
            (Covadj(1,:)+(Covadj(1,:)<=0)))) )';
         for lag=1:numlags
            Cor_xy=reshape(Coradj_slice(:,lag),numxs,numys);
            rho_vol(:,slice,lag)=reshape(conv2(ker_x,ker_y,Cor_xy,'same'),numpix,1);   
         end
      end 
   end
   
   % Smoothing rho betwen slices is done by straight matrix multiplication
   % by a toeplitz matrix K normalized so that the column sums are 1.
   
   if fwhm_rho>0
      fwhm_z=fwhm_rho/abs(Steps(3))*(100/df)^(1/3);
      ker_z=exp(-(0:(numslices-1)).^2*4*log(2)/fwhm_z^2);
      K=toeplitz(ker_z);
      K=K./(ones(numslices)*K);
      for lag=1:numlags
         rho_vol(:,:,lag)=squeeze(rho_vol(:,:,lag))*K;
      end
   end
   
   if which_stats(1,5)         
      out_rho.data=squeeze(reshape(rho_vol,numxs,numys,numslices,numlags));
      fmris_write_image(out_rho);
   end
   
end

if ~isnumeric(fwhm_rho)
   out_rho=fmris_read_image(fwhm_rho);
   rho_vol=squeeze(reshape(out_rho.data,numpix,numslices,numlags));
end

% Make space for results. Last dim of stat_slice: 1=Tstat, 2=Effect, 3=Sd:

stat_slice=zeros(numpix,numcontrasts,3);
if which_stats(1,4)
   Fstat_slice=zeros(numpix,1);
end
if which_stats(1,6)
   resid_slice=zeros(numpix,n);
end
if which_stats(1,7) | which_stats(1,9)
   wresid_slice=zeros(numpix,n);
end
if which_stats(1,8)
   A_slice=zeros(numpix,numlags);
end
if which_stats(1,9)   
   % setup for estimating the FWHM:
   I=numxs;
   J=numys;
   IJ=I*J;
   Im=I-1;
   Jm=J-1;
   nxy=conv2(ones(Im,Jm),ones(2));
   f=zeros(I,J);
   r=zeros(I,J);
   ip=[0 1 0 1];
   jp=[0 0 1 1];
   is=[1 -1  1 -1];
   js=[1  1 -1 -1];
   D=2+(numslices>1);
   alphaf=-1/(2*D);
   alphar=1/2;
   Step=abs(prod(Steps(1:D)))^(1/D);
end

% Set up for second loop:

X_type=[ones(1,sum(num_hrf_bases==1))*1 ones(1,sum(num_hrf_bases==2))*2 ...
      ones(1,sum(num_hrf_bases==2))*3 ones(1,numtrends)*4 ]
find_X_is_u1=find(X_type==2)         
find_X_is_u2=find(X_type==3)         
find_X_is_u12=[find_X_is_u1 find_X_is_u2]         
find_X_is_mag=[find(X_type==1) find(X_type==2) find(X_type==4)]

find_contrast_is_mag=find(~contrast_is_delay)
find_response_is_mag=[find(num_hrf_bases==1) find(num_hrf_bases==2) ...
      numresponses+(1:numtrends)]
contr_mag=contrasts(find_contrast_is_mag,find_response_is_mag)

find_contrast_is_delay=find(contrast_is_delay)
find_response_is_delay=find(num_hrf_bases==2)
contr_delay=contrast(find_contrast_is_delay,find_response_is_delay)
numcontr_delay=length(find_contrast_is_delay)
contr_delay2=repmat(contr_delay,1,2)';
contr_delay_is_1_col=zeros(numcontr_delay,1);
find_delay_is_1_col=ones(numcontr_delay,1);
for i=1:numcontr_delay
   pos=find(contr_delay(i,:)~=0);
   if length(pos)==1
      contr_delay_is_1_col(i)=1;
      find_delay_is_1_col(i)=pos;
   end
end
contr_delay_is_1_col 
find_delay_is_1_col

% Fit a tangent function to the basis coefficients, W:
cv=zeros(numresponses,1);
dd0v=zeros(numresponses,1);
for k=1:numresponses
   delta=X_cache.W(:,k,5);
   R=X_cache.W(:,k,basis2)./X_cache.W(:,k,basis1);
   ddelta=gradient(delta)./gradient(R);
   dd0=ddelta(delta==0);
   c=max(delta)/(pi/2);
   deltahat=atan(R/c*dd0)*c;
   for niter=1:5
      c=pinv(deltahat/c-cos(deltahat/c).^2.*R/c*dd0)*(delta-deltahat)+c;
      deltahat=atan(R/c*dd0)*c;
   end
   cv(k)=c;
   dd0v(k)=dd0;
end
C=cv(find_response_is_delay);
C*pi/2
Dd0=dd0v(find_response_is_delay)

% Second loop over voxels to get statistics:

drho=0.01;
Y=zeros(n,numpix);
for slice=1:numslices
   Second_pass_slice=slice
   if ~isempty(X_cache)
      X=[squeeze(X_cache.X(keep,num_hrf_bases==1,1,slice)) ... 
            squeeze(X_cache.X(keep,num_hrf_bases==2,basis1,slice)) ...
            squeeze(X_cache.X(keep,num_hrf_bases==2,basis2,slice)) ...
            Trend(:,:,slice)];
   else
      X=Trend(:,:,slice);
   end
   Xstar=X;
   Df=dfs(slice);
   
   % bias and continuity corrections for estimating the FWHM:
   
   if which_stats(1,9)
      df_resid=Df;
      dr=df_resid/Df;
      dv=df_resid-dr-(0:D-1);
      if df_resid>df_limit
         % constants for arbitrary filter method:
         biasf=exp(sum(gammaln(dv/2+alphaf)-gammaln(dv/2)) ...
            +gammaln(Df/2-D*alphaf)-gammaln(Df/2))*dr^(-D*alphaf);
         biasr=exp(sum(gammaln(dv/2+alphar)-gammaln(dv/2)) ...
            +gammaln(Df/2-D*alphar)-gammaln(Df/2))*dr^(-D*alphar);
         correctf=1/4;
         correctr=1/4;
      else
         % constants for filter aligned with axes method:
         biasf=exp((gammaln(dv(1)/2+alphaf)-gammaln(dv(1)/2))*D+ ...
            +gammaln(Df/2-D*alphaf)-gammaln(Df/2))*dr^(-D*alphaf);
         biasr=exp((gammaln(dv(1)/2+alphar)-gammaln(dv(1)/2))*D+ ...
            +gammaln(Df/2-D*alphar)-gammaln(Df/2))*dr^(-D*alphar);
         correctf=(df_resid+((4*D-2)*alphaf-3)*dr)/(df_resid-dr+2*alphaf)/4;
         correctr=(df_resid+((4*D-2)*alphar-3)*dr)/(df_resid-dr+2*alphar)/4;
      end
      consf=(4*log(2))^(-D*alphaf)/biasf*Step;
      consr=(4*log(2))^(-D*alphar)/biasr;
   end
   
   % read in data:
   
   if numfiles==1
      d=fmris_read_image(input_file,slice,keep);
      Y=reshape(d.data,numpix,n)';
   else
      for i=1:n
         d=fmris_read_image(input_file(keep(i),:),slice,1);
         Y(i,:)=reshape(d.data,1,numpix);
      end
   end
   
   % Convert to percent:
   
   if ispcnt
      for i=1:n
         Y(i,:)=Y(i,:)*(100/spatial_av(keep(i)));
      end
   end
   
   if numlags==1
      % bin rho to intervals of length drho, avoiding -1 and 1:
      irho=round(rho_vol(:,slice)/drho)*drho;
      irho=min(irho,1-drho);
      irho=max(irho,-1+drho);
   else
      % use dummy unique values so every pixel is analysed seperately:
      irho=(1:numpix)';
   end
   
   for rho=unique(irho)'
      pix=find(irho==rho);
      numrhos=length(pix);
      Ystar=Y(:,pix);
      if numlags==1
         factor=1./sqrt(1-rho^2);
         Ystar(k1,:)=(Y(k1,pix)-rho*Y(k1-1,pix))*factor;
         Xstar(k1,:)=(X(k1,:)-rho*X(k1-1,:))*factor;
      else
         Coradj_pix=squeeze(rho_vol(pix,slice,:));
         [Ainvt posdef]=chol(toeplitz([1 Coradj_pix']));
         nl=size(Ainvt,1);
         A=inv(Ainvt');
         if which_stats(1,8)
            A_slice(pix,1:(nl-1))=-A(nl,(nl-1):-1:1)/A(nl,nl);
         end
         B=ones(n-nl,1)*A(nl,:);
         Vmhalf=spdiags(B,1:nl,n-nl,n);
         Ystar=zeros(n,1);
         Ystar(1:nl)=A*Y(1:nl,pix);
         Ystar((nl+1):n)=Vmhalf*Y(:,pix);
         Xstar(1:nl,:)=A*X(1:nl,:);
         Xstar((nl+1):n,:)=Vmhalf*X;
      end
      pinvXstar=pinv(Xstar);
      betahat=pinvXstar*Ystar;
      resid=Ystar-Xstar*betahat;
      if which_stats(1,6)
         resid_slice(pix,:)=(Y(:,pix)-X*betahat)';
      end
      SSE=sum(resid.^2,1);
      sd=sqrt(SSE/Df);
      if which_stats(1,7) | which_stats(1,9)
         sdd=(sd>0)./(sd+(sd<=0))/sqrt(Df);
         wresid_slice(pix,:)=(resid.*repmat(sdd,n,1))';
      end
      V=pinvXstar*pinvXstar';
      sdbetahat=sqrt(diag(V))*sd;
      T0=betahat./(sdbetahat+(sdbetahat<=0)).*(sdbetahat>0);
      
      if any(~contrast_is_delay)
         
         % estimate magnitudes:
         
         mag_ef=contr_mag*betahat(find_X_is_mag,:);
         mag_sd=sqrt(diag(contr_mag*V(find_X_is_mag,find_X_is_mag)*contr_mag'))*sd;
         stat_slice(pix,find_contrast_is_mag,2)=mag_ef';
         stat_slice(pix,find_contrast_is_mag,3)=mag_sd';
         stat_slice(pix,find_contrast_is_mag,1)= ...
            (mag_ef./(mag_sd+(mag_sd<=0)).*(mag_sd>0))';
         
         if which_stats(1,4)
            cVcinv=pinv(contr_mag*V(find_X_is_mag,find_X_is_mag)*contr_mag');
            SST=sum((cVcinv*mag_ef).*mag_ef,1);
            Fstat_slice(pix,1)=(SST./(SSE+(SSE<=0)).*(SSE>0)/p*Df)';
         end      
      end
      
      if any(contrast_is_delay)
         
         % estimate delays:
         
         betaw1=betahat(find_X_is_u1,:);
         betaw2=betahat(find_X_is_u2,:);
         inv_betaw1=(betaw1~=0)./(betaw1+(betaw1==0));
         rhat=betaw2.*inv_betaw1;
         T0_2=T0(find_X_is_u1,:).^2;
         c0=T0_2./(T0_2+1);
         rhat_T=rhat.*c0;
         delay=atan(rhat_T.*repmat(Dd0./C,1,numrhos)).*repmat(C,1,numrhos);
         del_ef=contr_delay*delay;
         
         % estimates of sd of delay:
         
         drs=cos(delay./repmat(C,1,numrhos)).^2.*repmat(Dd0,1,numrhos);
         gdot1=c0./(T0_2+1).*(1-T0_2).*rhat.*inv_betaw1.*drs;
         gdot2=c0.*inv_betaw1.*drs;
         gdot=[gdot1; gdot2]; 
         gdotc=kron(gdot,ones(1,numcontr_delay)).*repmat(contr_delay2,1,numrhos);
         Vcontr=sum((V(find_X_is_u12,find_X_is_u12)*gdotc).*gdotc,1);
         del_sd=reshape(sqrt(Vcontr),numcontr_delay,numrhos) ...
            .*repmat(sd,numcontr_delay,1);
         
         % write out delays:
         
         stat_slice(pix,find_contrast_is_delay,2)=del_ef';
         stat_slice(pix,find_contrast_is_delay,3)=del_sd';
         contr_delay_is_1_cols=repmat(contr_delay_is_1_col,1,numrhos);
         sdbetaw2=sdbetahat(find_X_is_u2,:);
         T1=betaw2./(sdbetaw2+(sdbetaw2<=0)).*(sdbetaw2>0);
         stat_slice(pix,find_contrast_is_delay,1)= ...
            (contr_delay_is_1_cols.*T1(find_delay_is_1_col,:)+ ...
            ~contr_delay_is_1_cols.*del_ef./(del_sd+(del_sd<=0)).*(del_sd>0))';
      end
   end
   
   % Write out results; make sure T statistics don't exceed 100:
   
   stat_slice(:,:,1)=min(stat_slice(:,:,1),100);
   
   for i=1:numcontrasts
      for stat=1:3
         if which_stats(i,stat) 
            out(i,stat).data=reshape(stat_slice(:,i,stat),numxs,numys);
            fmris_write_image(out(i,stat),slice,1);
         end
      end
   end
   if which_stats(1,4)
      Fstat_slice=min(Fstat_slice,10000);
      out_Fstat.data=reshape(Fstat_slice,numxs,numys);
      fmris_write_image(out_Fstat,slice,1);
   end  
   if which_stats(1,6)
      out_resid.data=reshape(resid_slice,numxs,numys,n);
      fmris_write_image(out_resid,slice,1:n);
   end  
   if which_stats(1,7)
      out_wresid.data=reshape(wresid_slice,numxs,numys,n);
      fmris_write_image(out_wresid,slice,1:n);
   end  
   if which_stats(1,8)
      out_A.data=reshape(A_slice,numxs,numys,numlags);
      fmris_write_image(out_A,slice,1:numlags);
   end 
   if which_stats(1,9)
      
      % Finds an estimate of the fwhm for each of the 8 cube corners surrounding
      % a voxel, then averages. 
      
      if slice==1
         u=reshape(wresid_slice,I,J,n);
         ux=diff(u,1,1);
         uy=diff(u,1,2);
         Axx=sum(ux.^2,3);
         Ayy=sum(uy.^2,3);
         if D==2
            for index=1:4
               i=(1:Im)+ip(index);
               j=(1:Jm)+jp(index);
               axx=Axx(:,j);
               ayy=Ayy(i,:);
               if df_resid>df_limit
                  axy=sum(ux(:,j,:).*uy(i,:,:),3)*is(index)*js(index);
                  detlam=(axx.*ayy-axy.^2);
                  % Continuity correction:
                  detlamtrlaminvlam2=(axx+2*axy+ayy).*axx.*ayy- ...
                     4*(axx-axy+ayy).*axy.^2;
                  detlamf=detlam+correctf*detlamtrlaminvlam2;
                  detlamr=detlam+correctr*detlamtrlaminvlam2;
               else
                  detlamf=axx.*ayy.*(1+correctf*(axx+ayy));
                  detlamr=axx.*ayy.*(1+correctr*(axx+ayy));
               end
               f(i,j)=f(i,j)+(detlamf>0).*(detlamf+(detlamf<=0)).^alphaf;
               r(i,j)=r(i,j)+(detlamr>0).*(detlamr+(detlamr<=0)).^alphar;
            end
         end
      else 
         uz=reshape(wresid_slice,I,J,n)-u;
         Azz=sum(uz.^2,3);
         % The 4 upper cube corners:
         for index=1:4
            i=(1:Im)+ip(index);
            j=(1:Jm)+jp(index);
            axx=Axx(:,j);
            ayy=Ayy(i,:);
            azz=Azz(i,j);
            if Df>df_limit
               axy=sum(ux(:,j,:).*uy(i,:,:),3)*is(index)*js(index);
               axz=sum(ux(:,j,:).*uz(i,j,:),3)*is(index);
               ayz=sum(uy(i,:,:).*uz(i,j,:),3)*js(index);
               detlam=(axx.*ayy-axy.^2).*azz-(axz.*ayy-2*axy.*ayz).*axz-axx.*ayz.^2;
               s1=axx+ayy+azz;
               s2=axy+ayz+axz;
               % Continuity correction:
               detlamtrlaminvlam2=(s1+2*s2).*axx.*ayy.*azz+ ...
                  4*(2*s1-s2).*axz.*axy.*ayz- ...
                  (axx+4*(ayy-ayz+azz)).*axx.*ayz.^2- ...
                  (ayy+4*(axx-axz+azz)).*ayy.*axz.^2- ...
                  (azz+4*(axx-axy+ayy)).*azz.*axy.^2- ...
                  2*(axx.*ayz.*(axz.*ayy+axy.*azz)+ayy.*azz.*axy.*axz);
               detlamf=detlam+correctf*detlamtrlaminvlam2;
               detlamr=detlam+correctr*detlamtrlaminvlam2;
            else
               detlamf=axx.*ayy.*azz.*(1+correctf*(axx+ayy+azz));
               detlamr=axx.*ayy.*azz.*(1+correctr*(axx+ayy+azz));
            end
            f(i,j)=f(i,j)+(detlamf>0).*(detlamf+(detlamf<=0)).^alphaf;
            r(i,j)=r(i,j)+(detlamr>0).*(detlamr+(detlamr<=0)).^alphar;
         end
         f=consf/((slice>2)+1)*f./nxy;
         r=consr/((slice>2)+1)*r./nxy;
         out_fwhm.data=f.*(f<50)+50*(f>=50);
         fmris_write_image(out_fwhm,slice-1,1);
         out_fwhm.data=r;
         fmris_write_image(out_fwhm,slice-1,2);
         
         f=zeros(I,J);
         r=zeros(I,J);
         u=reshape(wresid_slice,I,J,n);
         ux=diff(u,1,1);
         uy=diff(u,1,2);
         Axx=sum(ux.^2,3);
         Ayy=sum(uy.^2,3);
         % The 4 lower cube corners:
         for index=1:4
            i=(1:Im)+ip(index);
            j=(1:Jm)+jp(index);
            axx=Axx(:,j);
            ayy=Ayy(i,:);
            azz=Azz(i,j);
            if Df>df_limit
               axy=sum(ux(:,j,:).*uy(i,:,:),3)*is(index)*js(index);
               axz=-sum(ux(:,j,:).*uz(i,j,:),3)*is(index);
               ayz=-sum(uy(i,:,:).*uz(i,j,:),3)*js(index);
               detlam=(axx.*ayy-axy.^2).*azz-(axz.*ayy-2*axy.*ayz).*axz-axx.*ayz.^2;
               s1=axx+ayy+azz;
               s2=axy+ayz+axz;
               % Continuity correction:
               detlamtrlaminvlam2=(s1+2*s2).*axx.*ayy.*azz+ ...
                  4*(2*s1-s2).*axz.*axy.*ayz- ...
                  (axx+4*(ayy-ayz+azz)).*axx.*ayz.^2- ...
                  (ayy+4*(axx-axz+azz)).*ayy.*axz.^2- ...
                  (azz+4*(axx-axy+ayy)).*azz.*axy.^2- ...
                  2*(axx.*ayz.*(axz.*ayy+axy.*azz)+ayy.*azz.*axy.*axz);
               detlamf=detlam+correctf*detlamtrlaminvlam2;
               detlamr=detlam+correctr*detlamtrlaminvlam2;
            else
               detlamf=axx.*ayy.*azz.*(1+correctf*(axx+ayy+azz));
               detlamr=axx.*ayy.*azz.*(1+correctr*(axx+ayy+azz));
            end
            f(i,j)=f(i,j)+(detlamf>0).*(detlamf+(detlamf<=0)).^alphaf;
            r(i,j)=r(i,j)+(detlamr>0).*(detlamr+(detlamr<=0)).^alphar;
         end
         if slice==numslices
            f=consf*f./nxy;
            r=consr*r./nxy;
            out_fwhm.data=f.*(f<50)+50*(f>=50);
            fmris_write_image(out_fwhm,slice,1);
            out_fwhm.data=r;
            fmris_write_image(out_fwhm,slice,2);
         end
      end
   end
end

return

