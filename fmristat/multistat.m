function [df, df_resid] = multistat(input_files_Y,input_files_sd, ...
   input_files_df,input_files_fwhm,X,contrast,output_file_base, ...
   which_stats,fwhm_varatio,niter,df_limit)

%MULTISTAT fits a mixed effects linear model.
%
% Combines effects (E) and their standard errors (S) using a linear mixed 
% effects model:     E = X b + e_fixed + e_random,     where
%    b is a vector of unknown coefficients,
%    e_fixed  is normal with mean zero, standard deviation S,
%    e_random is normal with mean zero, standard deviation sigma (unknown).
% The model is fitted by REML using the EM algorithm with NITER iterations. 
%
% Gives the conjunction (minimum) of the random effects T stats for the data.
%
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
% WARNING: Multistat is very slow if the number of columns of X is more than 
% 1 and less than the number of input files, and INPUT_FILES_SD is not empty,
% since it loops over voxels, rather than doing calculations in parallel.
% 
% [DF DF_RESID] = MULTISTAT( INPUT_FILES_EFFECT , INPUT_FILES_SDEFFECT , 
%           INPUT_FILES_DF [, INPUT_FILES_FWHM [, X [, 
%           [, CONTRAST [, OUTPUT_FILE_BASE [, WHICH_STATS 
%           [, FWHM_VARATIO [, NITER [, DF_LIMIT ]]]]]]]]] )
% 
% INPUT_FILES_EFFECT is the input fmri effect files, the dependent variables,
% usually the _ef.img or _ef.mnc files, padded with extra blanks if necessary; 
% they will be removed before use.
% 
% INPUT_FILES_SDEFFECT is the input fmri sdeffect files, the standard
% deviations of the dependent variables, usually the _sd.img or _sd.mnc files,
% padded with extra blanks if necessary; they will be removed before use.
% If INPUT_FILES_SDEFFECT=[], then INPUT_FILES_SDEFFECT is assumed to be
% zero for all voxels, INPUT_FILES_DF is set to Inf, and FWHM_VARATIO now
% smoothes the voxel sd. This allows multistat to duplicate DOT for 
% analysing PET data (with smoothed voxel sd, rather than pooled sd).
% 
% INPUT_FILES_DF is the row vector of degrees of freedom of the input files
% as printed out by fmrilm.m. If it is a scalar, it is repeated to fill.
%
% INPUT_FILES_FWHM is the fwhm in mm of the original fMRI data. It is only  
% used to calculate the degrees of freedom, printed out at the beginning. 
% Default is 6.
%
% X is the design matrix, whose rows are the files, and columns
% are the explanatory (independent) variables of interest. 
% Default is X=[1; 1; 1; ..1] which just averages the files. If the
% rank of X equals the number of files, e.g. if X is square, then 
% the random effects cannot be estinmated, but the fixed effects
% sd's can be used for the standard error. This is done very quickly.
% 
% CONTRAST is a matrix whose rows are contrasts for the statistic images.
% Default is [1 0 ... 0], i.e. it picks out the first column of X.
% 
% OUTPUT_FILE_BASE: Base for output statistics. Default is the first 
% INPUT_FILES_EFFECT name minus the extension.
%
% WHICH_STATS: logical matrix indicating output statistics by 1.
% Rows are contrasts, columns correspond to:
%   1: T statistic image, OUTPUT_FILE_BASE_t.mnc or .img. The degrees 
%      of freedom is DF. Note that tstat=effect/sdeffect. 
%   2: effect image, OUTPUT_FILE_BASE_ef.mnc or .img.
%   3: standard deviation of the effect, OUTPUT_FILE_BASE_sd.mnc or .img. 
%   4: ratio of random effects standard deviation to fixed effects 
%      standard deviation, OUTPUT_FILE_BASE_sdratio.mnc or .img.
%   5: conjunction (minimum) of the T statistics for the data, i.e. 
%      min(INPUT_FILES_EFFECT/sd), in OUTPUT_FILE_BASE_conj.mnc or .img.
%   6: the residuals from the model, OUTPUT_FILE_BASE_resid.mnc or .img,
%      only for the non-excluded frames (warning: uses lots of memory). 
%   7: the whitened residuals from the model normalized by dividing
%      by ( estimated standard deviation * sqrt(DF_RESID) ), 
%      OUTPUT_FILE_BASE_wresid.mnc or .img (warning: uses lots of memory).
%   8: Not used yet.
%   9: the estimated FWHM in OUTPUT_FILE_BASE_FWHM.mnc or .img - note this
%      uses more time and memory.
% If the number of columns is less than 9 it is padded with zeros. 
% If there is only one row, it is repeated for all contrasts.
% Default is 1. 
% 
% FWHM_VARATIO is the fwhm in mm of the Gaussian filter used to smooth the
% ratio of the random effects variance divided by the fixed effects variance.
%  - 0 will do no smoothing, and give a purely random effects analysis;
%  - Inf will do complete smoothing to a global ratio of one, giving a 
%    purely fixed effects analysis. 
% The higher the FWHM_VARATIO, the higher the ultimate degrees of
% freedom DF of the tstat image (printed out at the beginning of the run), 
% and the more sensitive the test. However too much smoothing will
% bias the results. The program prints and returns DF as its value.
% Alternatively, if FWHM_VARATIO is negative, it is taken as the desired
% df, and the fwhm is chosen to get as close to this as possible (if fwhm>50,  
% fwhm=Inf). Default is -100, i.e. the fwhm is chosen to achieve 100 df. 
%
% NITER is the number of iterations of the EM algorithm. Default is 10.
%
% DF_LIMIT controls which method is used for estimating FWHM. If DF_RESID > 
% DF_LIMIT, then the FWHM is calculated assuming the Gaussian filter is 
% arbitrary. However if DF is small, this gives inaccurate results, so
% if DF_RESID <= DF_LIMIT, the FWHM is calculated assuming that the axes of
% the Gaussian filter are aligned with the x, y and z axes of the data. 
% Default is 4. 
%
% DF is the approximate (residual) degrees of freedom of the T statistics.
%
% DF_RESID is the residual degrees of freedom of a random effects analysis.

%############################################################################
% COPYRIGHT:   Copyright 2000 K.J. Worsley and C. Liao, 
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

if nargin<4; input_files_fwhm=[]; end
if nargin<5; X=[]; end
if nargin<6; contrast=[]; end
if nargin<7; output_file_base=[]; end
if nargin<8; which_stats=[]; end
if nargin<9; fwhm_varatio=[]; end
if nargin<10; niter=[]; end
if nargin<11; df_limit=[]; end

parent_file=deblank(input_files_Y(1,:));
if isempty(input_files_fwhm); input_files_fwhm=6; end
if isempty(X); X=ones(size(input_files_Y,1),1); end
if isempty(contrast); contrast=[1 zeros(1,size(X,2)-1)]; end
if isempty(output_file_base); output_file_base=parent_file(1:(min(findstr('.',parent_file))-1)); end
if isempty(which_stats); which_stats=1; end
if isempty(fwhm_varatio); fwhm_varatio=-100; end
if isempty(niter); niter=10; end
if isempty(df_limit); df_limit=4; end

% Open images:

n=size(input_files_Y,1)
d=fmris_read_image(parent_file);
numslices=d.dim(3)
numys=d.dim(2)
numxs=d.dim(1)
numpix=numys*numxs;
Steps=d.vox
numcolX=size(contrast,2)
D=2+(numslices>1)

% Open files for writing:

ext=parent_file(min(findstr('.',parent_file))+(1:3));
out.parent_file=parent_file;
numcontrasts=size(contrast,1)
which_stats=[which_stats zeros(size(which_stats,1),9-size(which_stats,2))];
if size(which_stats,1)==1
   which_stats=repmat(which_stats,numcontrasts,1);
end
which_stats

X2=X'.^2;
p=rank(X)
df_resid=n-p
if df_resid>0
   
   % Degrees of freedom is greater than zero, so do mixed effects analysis:
   
   if isempty(input_files_sd)
      input_files_df=Inf
   end
   if length(input_files_df)==1
      input_files_df=ones(1,n)*input_files_df
   end
   df_fixed=sum(input_files_df)
   
   % Set-up for loop over voxels:
   
   Y=zeros(n,numpix);
   S=ones(n,numpix);
   varfix=ones(1,numpix);
   varatio_vol=zeros(numpix, numslices);
   
   if fwhm_varatio<0
      
      % find fwhm to achieve target df:
      
      df_target=-fwhm_varatio;
      if df_target<=df_resid
         fwhm_varatio=0
      elseif df_target>=df_fixed
         fwhm_varatio=Inf
      else
         ff=0:1:100;
         dfs=[];
         for f=ff;
            df=regularized_df(f, ...
               D,Steps,numslices,df_resid,df_fixed,input_files_fwhm);
            dfs=[dfs df];
         end
         fwhm_varatio=interp1(dfs,ff,df_target);
         if isnan(fwhm_varatio)
            fwhm_varatio=Inf;
         end
         if fwhm_varatio>50
            fwhm_varatio=Inf;
         end
         fwhm_varatio
      end
   end
   
   if fwhm_varatio<Inf
      
      % Smoothing varatio in slice is done using conv2 with a kernel ker_xy.
      % Smoothing betwen slices is done by straight matrix multiplication
      % by a toeplitz matrix K normalized so that the column sums are 1.
      
      [df, df_indep, df_correl, ker_x, ker_y, K]=regularized_df(fwhm_varatio, ...
         D,Steps,numslices,df_resid,df_fixed,input_files_fwhm);
      df=round(df)
      
      % Now loop over voxels to get variance ratio, varatio:
      
      Sreduction=0.99
      for slice=1:numslices
         slice
         for ifile=1:n
            d=fmris_read_image(deblank(input_files_Y(ifile,:)),slice,1);
            Y(ifile,:)=reshape(d.data,1,numpix);
         end
         sigma2=sum((Y-X*pinv(X)*Y).^2,1)/df_resid;
         if ~isempty(input_files_sd)
            for ifile=1:n
               d=fmris_read_image(deblank(input_files_sd(ifile,:)),slice,1);
               S(ifile,:)=reshape(d.data,1,numpix).^2;
            end
            minS=min(S)*Sreduction;
            Sm=S-ones(n,1)*minS;
            if size(X,2)==1
               % When X is a vector, calculations can be done in parallel:
               for iter=1:niter
                  Sms=Sm+ones(n,1)*sigma2;
                  W=(Sms>0)./(Sms+(Sms<=0));
                  X2W=X2*W;
                  XWXinv=(X2W>0)./(X2W+(X2W<=0));
                  betahat=XWXinv.*(X'*(W.*Y));
                  R=W.*(Y-X*betahat);
                  ptrS=p+sum(Sm.*W,1)-(X2*(Sm.*W.^2)).*XWXinv;
                  sigma2=(sigma2.*ptrS+(sigma2.^2).*sum(R.^2,1))/n; 
               end
               sigma2=sigma2-minS;
            else
               % Otherwise when X is a matrix, we have to loop over voxels:
               for pix=1:numpix
                  sigma2_pix=sigma2(pix);
                  Sm_pix=Sm(:,pix);
                  Y_pix=Y(:,pix);
                  for iter=1:niter
                     Sms=Sm_pix+sigma2_pix;
                     W=(Sms>0)./(Sms+(Sms<=0));
                     Whalf=diag(sqrt(W));
                     WhalfX=Whalf*X;
                     pinvX=pinv(WhalfX);
                     WhalfY=Whalf*Y_pix;
                     betahat=pinvX*WhalfY;
                     R=WhalfY-WhalfX*betahat;
                     SW=diag(Sm_pix.*W);
                     ptrS=p+sum(Sm_pix.*W,1)-sum(sum((SW*WhalfX).*pinvX'));   
                     sigma2_pix=(sigma2_pix.*ptrS+ ...
                        (sigma2_pix.^2).*sum(W.*(R.^2),1))/n; 
                  end
                  sigma2(pix)=sigma2_pix-minS(pix);
               end         
            end
            varfix=input_files_df*S/df_fixed;   
         end
         varatio=sigma2./(varfix+(varfix<=0)).*(varfix>0);
         varatio_slice=reshape(varatio,numxs,numys);
         varatio_slice=conv2(ker_x,ker_y,varatio_slice,'same');
         varatio_vol(:,slice)=reshape(varatio_slice,numpix,1);        
      end
      
      % Smooth varatio in z direction:
      
      if fwhm_varatio>0 
         varatio_vol=varatio_vol*K;
      end
      
   else
      df=round(df_fixed)
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
      alphaf=-1/(2*D);
      alphar=1/2;
      dr=df_resid/df;
      dv=df_resid-dr-(0:D-1);
      if df_resid>df_limit
         % constants for arbitrary filter method:
         biasf=exp(sum(gammaln(dv/2+alphaf)-gammaln(dv/2)) ...
            +gammaln(df/2-D*alphaf)-gammaln(df/2))*dr^(-D*alphaf);
         biasr=exp(sum(gammaln(dv/2+alphar)-gammaln(dv/2)) ...
            +gammaln(df/2-D*alphar)-gammaln(df/2))*dr^(-D*alphar);
         correctf=1/4;
         correctr=1/4;
      else
         % constants for filter aligned with axes method:
         biasf=exp((gammaln(dv(1)/2+alphaf)-gammaln(dv(1)/2))*D+ ...
            +gammaln(df/2-D*alphaf)-gammaln(df/2))*dr^(-D*alphaf);
         biasr=exp((gammaln(dv(1)/2+alphar)-gammaln(dv(1)/2))*D+ ...
            +gammaln(df/2-D*alphar)-gammaln(df/2))*dr^(-D*alphar);
         correctf=(df_resid+((4*D-2)*alphaf-3)*dr)/(df_resid-dr+2*alphaf)/4;
         correctr=(df_resid+((4*D-2)*alphar-3)*dr)/(df_resid-dr+2*alphar)/4;
      end
      Step=abs(prod(Steps(1:D)))^(1/D);
      consf=(4*log(2))^(-D*alphaf)/biasf*Step;
      consr=(4*log(2))^(-D*alphar)/biasr;
   end
   
   pinvX=pinv(X);
   ncpinvX=sqrt(sum((contrast*pinvX).^2,2));
   
   % Second loop over slices to get statistics:
   
   for slice=1:numslices
      slice            
      for ifile=1:n
         d=fmris_read_image(deblank(input_files_Y(ifile,:)),slice,1);
         Y(ifile,:)=reshape(d.data,1,numpix);
      end
      if ~isempty(input_files_sd)
         for ifile=1:n
            d=fmris_read_image(deblank(input_files_sd(ifile,:)),slice,1);
            S(ifile,:)=reshape(d.data,1,numpix).^2;
         end
         varfix=input_files_df*S/df_fixed;   
         sigma2=varatio_vol(:,slice)'.*varfix;
         Sigma=S+ones(n,1)*sigma2;
         % To ensure that Sigma is always positive:
         if fwhm_varatio>0
            Sigma=max(Sigma,S*0.25);
         end
         W=(Sigma>0)./(Sigma+(Sigma<=0));
         if size(X,2)==1
            % Do calculations in parallel:
            X2W=X2*W;
            XWXinv=(X2W>0)./(X2W+(X2W<=0));
            betahat=XWXinv.*(X'*(W.*Y));
            sdeffect_slice=contrast*sqrt(XWXinv);
         else
            % Loop over pixels:
            betahat=zeros(size(X,2),numpix);
            sdeffect_slice=zeros(numcontrasts,numpix);
            for pix=1:numpix
               Whalf=diag(sqrt(W(:,pix)));
               pinvX=pinv(Whalf*X);
               betahat(:,pix)=pinvX*Whalf*Y(:,pix);
               sdeffect_slice(:,pix)=sqrt(sum((contrast*pinvX).^2,2));
            end  
         end
      else
         betahat=pinvX*Y;
         sigma2=varatio_vol(:,slice)';            
         sdeffect_slice=ncpinvX*sqrt(sigma2);
         Sigma=ones(n,1)*sigma2;
         W=(Sigma>0)./(Sigma+(Sigma<=0));
      end
      effect_slice=contrast*betahat;
      tstat_slice=effect_slice./(sdeffect_slice+(sdeffect_slice<=0)) ...
         .*(sdeffect_slice>0);
      sdratio_slice=sqrt(varatio_vol(:,slice)+1);
      
      % output:
      
      for k=1:numcontrasts
         if which_stats(k,1)
            out.file_name=[deblank(output_file_base(k,:)) '_t.' ext];
            out.dim=[numxs numys numslices 0];
            out.data=reshape(tstat_slice(k,:),numxs,numys);
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,2)
            out.file_name=[deblank(output_file_base(k,:)) '_ef.' ext];
            out.dim=[numxs numys numslices 0];
            out.data=reshape(effect_slice(k,:),numxs,numys);
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,3)
            out.file_name=[deblank(output_file_base(k,:)) '_sd.' ext];
            out.dim=[numxs numys numslices 0];
            out.data=reshape(sdeffect_slice(k,:),numxs,numys);
            fmris_write_image(out,slice,1);
         end
      end
      if which_stats(1,4)
         out.file_name=[deblank(output_file_base(1,:)) '_sdratio.' ext];
         out.dim=[numxs numys numslices 0];
         out.data=reshape(sdratio_slice,numxs,numys);
         fmris_write_image(out,slice,1);
      end
      if which_stats(1,5)
         out.file_name=[deblank(output_file_base(1,:)) '_conj.' ext];
         out.dim=[numxs numys numslices 0];
         out.data=reshape(min(sqrt(W).*Y,[],1),numxs,numys);
         fmris_write_image(out,slice,1);
      end
      if which_stats(1,6)
         out.file_name=[deblank(output_file_base(1,:)) '_resid.' ext];
         out.dim=[numxs numys numslices n];
         out.data=reshape((Y-X*betahat)',numxs,numys,n);
         fmris_write_image(out,slice,1:n);
      end  
      if which_stats(1,7) | which_stats(9)
         wresid_slice=(sqrt(W/df_resid).*(Y-X*betahat))';
      end
      if which_stats(1,7) 
         out.file_name=[deblank(output_file_base(1,:)) '_wresid.' ext];
         out.dim=[numxs numys numslices n];
         out.data=reshape(wresid_slice,numxs,numys,n);
         fmris_write_image(out,slice,1:n);
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
               if df_resid>df_limit
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
            out.file_name=[deblank(output_file_base(1,:)) '_fwhm.' ext];
            out.dim=[numxs numys numslices 2];
            out.data=f.*(f<50)+50*(f>=50);
            fmris_write_image(out,slice-1,1);
            out.data=r;
            fmris_write_image(out,slice-1,2);
            
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
               if df_resid>df_limit
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
         end
         if slice==numslices
            f=consf*f./nxy;
            r=consr*r./nxy;
            out.file_name=[deblank(output_file_base(1,:)) '_fwhm.' ext];
            out.dim=[numxs numys numslices 2];
            out.data=f.*(f<50)+50*(f>=50);
            fmris_write_image(out_fwhm,slice,1);
            out.data=r;
            fmris_write_image(out_fwhm,slice,2);
         end
      end
   end
else
   % If degrees of freedom is zero, estimate effects by least squares,
   % and use the standard errors to estimate the sdeffect.
   
   Y=zeros(n,numpix);
   S=ones(n,numpix);
   cXinv=contrast*pinv(X);
   cXinv2=cXinv.^2;    
   for slice=1:numslices
      slice            
      for ifile=1:n
         d=fmris_read_image(deblank(input_files_Y(ifile,:)),slice,1);
         Y(ifile,:)=reshape(d.data,1,numpix);
      end
      if ~isempty(input_files_sd)
         for ifile=1:n
            d=fmris_read_image(deblank(input_files_sd(ifile,:)),slice,1);
            S(ifile,:)=reshape(d.data,1,numpix).^2;
         end
      end
      effect_slice=cXinv*Y;
      sdeffect_slice=sqrt(cXinv2*S);
      tstat_slice=effect_slice./(sdeffect_slice+(sdeffect_slice<=0)) ...
         .*(sdeffect_slice>0);
      for k=1:numcontrasts
         if which_stats(k,1)
            out.file_name=[deblank(output_file_base(k,:)) '_t.' ext];
            out.dim=[numxs numys numslices 0];
            out.data=reshape(tstat_slice(k,:),numxs,numys);
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,2)
            out.file_name=[deblank(output_file_base(k,:)) '_ef.' ext];
            out.dim=[numxs numys numslices 0];
            out.data=reshape(effect_slice(k,:),numxs,numys);
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,3)
            out.file_name=[deblank(output_file_base(k,:)) '_sd.' ext];
            out.dim=[numxs numys numslices 0];
            out.data=reshape(sdeffect_slice(k,:),numxs,numys);
            fmris_write_image(out,slice,1);
         end
      end
   end
   fprintf('Note: df_resid=0, so doing a fixed effects analysis.');
   fprintf('The approximate df is:\n');
   if ~isempty(input_files_df)
      df=sum(cXinv2,2).^2./((cXinv2.^2)*(1./input_files_df'))
   else
      df=Inf
   end
   
end

return



function [df, df_indep, df_correl, ker_x, ker_y, K]=regularized_df(fwhm_varatio, ...
   D,Steps,numslices,df_resid,df_fixed,input_files_fwhm);

if fwhm_varatio>0
   fwhm_x=fwhm_varatio/abs(Steps(1));
   ker_x=exp(-(-ceil(fwhm_x):ceil(fwhm_x)).^2*4*log(2)/fwhm_x^2);
   ker_x=ker_x/sum(ker_x);
   fwhm_y=fwhm_varatio/abs(Steps(2));
   ker_y=exp(-(-ceil(fwhm_y):ceil(fwhm_y)).^2*4*log(2)/fwhm_y^2);
   ker_y=ker_y/sum(ker_y);
   fwhm_z=fwhm_varatio/abs(Steps(3));
   ker_z=exp(-(0:(numslices-1)).^2*4*log(2)/fwhm_z^2);
   K=toeplitz(ker_z);
   K=K./(ones(numslices)*K);
else
   ker_x=1;
   ker_y=1;
   K=eye(numslices);
end
df_indep=df_resid/(sum(ker_x.^2)*sum(ker_y.^2)*sum(sum(K.^2))/numslices);
df_correl=df_resid*(2*(fwhm_varatio/input_files_fwhm)^2+1)^(D/2);
df=1/(1/min(df_indep,df_correl)+1/df_fixed);

return

