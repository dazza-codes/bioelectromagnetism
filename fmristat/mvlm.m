function mvlm(input_files, X, contrast, ...
   output_file_base, which_stats, df_limit)

%MVLM Multivariate linear model statistics such as Hotelling's T^2 
%
% MVLM( INPUT_FILES, X, CONTRAST, OUTPUT_FILE_BASE,  
%              WHICH_STATS [, DF_LIMIT])
%
% Fits multivariate linear model Y = X B + E, where Y is n x NVAR.
% At the moment NVAR<=3 only. 
% 
% INPUT_FILES: rows are the 4D input files of data for the dependent 
% variable Y, frames are NVAR variates. If the number of frames=1, 
% and the file name contains _dx, then
% the program assumes NVAR=3 and looks for the other 2 variates in 
% separate files by replacing _dx with _dy and _dz in the file name. 
%
% X is the design matrix, whose rows are the files, and columns
% are the explanatory (independent) variables of interest. 
% 
% CONTRAST is a matrix whose rows are contrasts for the effects.
%
% WHICH_STATS: logical matrix indicating output statistics by 1.
% Rows correspond to contrasts, columns correspond to:
%   1: T^2 statistic image, OUTPUT_FILE_BASE_T2.mnc or .img 
%   2: effect image, OUTPUT_FILE_BASE_ef.mnc or .img, with
%      frames for variates.
%   3: covariance matrix of the effect, OUTPUT_FILE_BASE_cov.mnc or 
%      .img, with frames for the upper triangle only, e.g. 
%      for 3 variates: 11, 12 13, 22, 23, 33.
%   4: Roy's maximum root, analogue of F statistic for all contrasts,
%      OUTPUT_FILE_BASE_Roy.mnc or .img. 
%   9: the estimated FWHM in OUTPUT_FILE_BASE_FWHM.mnc or .img
%      - note this uses more time and memory.
% If the number of columns is less than 9 it is padded with zeros. 
%
% OUTPUT_FILE_BASE: Rows are bases for output statistics for each
% contrast.
%
% DF_LIMIT controls which method is used for estimating FWHM. If DF > 
% DF_LIMIT, then the FWHM is calculated assuming the Gaussian filter is 
% arbitrary. However if DF is small, this gives inaccurate results, so
% if DF <= DF_LIMIT, the FWHM is calculated assuming that the axes of
% the Gaussian filter are aligned with the x, y and z axes of the data. 
% Default is 4. 

if nargin < 6
   df_limit=4;
end

n=size(input_files,1)
p=rank(contrast)
df=n-rank(X)
parent_file=deblank(input_files(1,:));
d=fmris_read_image(parent_file,0,0)
D=2+(d.dim(3)>1)
ext=parent_file(size(parent_file,2)+(-2:0))
out.parent_file=parent_file;

% number of variates
L=d.dim(4); 
if L==1 & ~isempty(strfind(parent_file,'_dx'))
   L=3; 
end
LL=round(L*(L+1)/2);

% Check for estimability:
tolerance=0.0000001;
NullSpaceX=null(X);
Cmhalf=diag(1./sqrt(diag(contrast*contrast')));
ContrastNullSpaceX=Cmhalf*contrast*NullSpaceX;
nonest=sum(abs(ContrastNullSpaceX)>tolerance,2);
if sum(nonest)>0
   fprintf(['Error: the following contrasts are nonestimable:']);
   RowsNonEstContrast=find(nonest>0)
   NonEstContrast=contrast(RowsNonEstContrast,:)
   NullSpaceX
   ContrastNullSpaceX
   return
end

numcontrasts=size(contrast,1)
which_stats=[which_stats zeros(size(which_stats,1),9-size(which_stats,2))];
if size(which_stats,1)==1
   which_stats=repmat(which_stats,numcontrasts,1);
end
which_stats

if which_stats(1,9)   
   wresid_slice=zeros(d.dim(1)*d.dim(2),n,L);
   diags=round(((1:L)-1).*(2*L-(1:L)+2)/2+1);
   I=d.dim(1);
   J=d.dim(2);
   D=3;
   df_resid=df;
   IJ=I*J;
   Im=I-1;
   Jm=J-1;
   nxy=repmat(conv2(ones(Im,Jm),ones(2)),[1 1 L]);
   f=zeros(I,J,L);
   r=zeros(I,J,L);
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
   Step=abs(prod(d.vox(1:D)))^(1/D);
   consf=(4*log(2))^(-D*alphaf)/biasf*Step;
   consr=(4*log(2))^(-D*alphar)/biasr;
end

pinvX=pinv(X);
covef=(contrast*pinvX)*(contrast*pinvX)';
covef_inv=pinv(covef);

numcolX=size(X,2);
Y=zeros(n,d.dim(1)*d.dim(2),L);
betahat=zeros(numcolX,d.dim(1)*d.dim(2),L);
ef=zeros(numcontrasts,d.dim(1)*d.dim(2),L);
cov=zeros(d.dim(1)*d.dim(2),LL);
compext=['_dx'; '_dy'; '_dz']; 
for slice=1:d.dim(3)
   slice
   for i=1:n
      for j=1:L
         if d.dim(4)>=L
            d=fmris_read_image(deblank(input_files(i,:)),slice,j);
         else
            input_file1=strrep(deblank(input_files(i,:)),compext(1,:),compext(j,:));
            d=fmris_read_image(input_file1,slice,1);
         end
         Y(i,:,j)=reshape(d.data,1,d.dim(1)*d.dim(2));
      end
   end
   
   for j=1:L
      betahat(:,:,j)=pinvX*Y(:,:,j);
      ef(:,:,j)=contrast*betahat(:,:,j);
      Y(:,:,j)=Y(:,:,j)-X*betahat(:,:,j);
   end
   
   for k=1:numcontrasts
      if which_stats(k,2)
         out.file_name=[deblank(output_file_base(k,:)) '_ef.' ext];
         out.dim=[d.dim(1:3) L];
         out.data=reshape(ef(k,:,:),d.dim(1),d.dim(2),L);
         fmris_write_image(out,slice,1:L);
      end
   end
   
   l=0;
   for j=1:L
      for k=j:L
         l=l+1;
         cov(:,l)=sum(Y(:,:,j).*Y(:,:,k),1)'/df;
      end
   end
   
   for k=1:numcontrasts
      if which_stats(k,3)
         out.file_name=[deblank(output_file_base(k,:)) '_cov.' ext];
         out.dim=[d.dim(1:3) LL];
         out.data=reshape(cov,d.dim(1),d.dim(2),LL)*covef(k,k);
         fmris_write_image(out,slice,1:LL);
      end
   end
   
   if any(which_stats(:,1)) | which_stats(1,4)
      if L==1
         det=cov(:,1);
         for k=1:numcontrasts
            cov_inv_ef(k,:,1)= ef(k,:,1)';
         end
      if L==2
         det=cov(:,1).*cov(:,3)-cov(:,2).^2;
         for k=1:numcontrasts
            cov_inv_ef(k,:,1)= ef(k,:,1)'.*cov(:,3)-ef(k,:,2)'.*cov(:,2);
            cov_inv_ef(k,:,2)=-ef(k,:,1)'.*cov(:,2)+ef(k,:,2)'.*cov(:,1);
         end
      end
      if L==3
         det=cov(:,1).*(cov(:,4).*cov(:,6)-cov(:,5).^2) ...
            -cov(:,2).*(cov(:,2).*cov(:,6)-cov(:,5).*cov(:,3)) ...
            +cov(:,3).*(cov(:,2).*cov(:,5)-cov(:,4).*cov(:,3));
         for k=1:numcontrasts
            cov_inv_ef(k,:,1)=ef(k,:,1)'.*(cov(:,4).*cov(:,6)-cov(:,5).^2) ...
               -ef(k,:,2)'.*(cov(:,2).*cov(:,6)-cov(:,5).*cov(:,3)) ...
               +ef(k,:,3)'.*(cov(:,2).*cov(:,5)-cov(:,4).*cov(:,3));
            cov_inv_ef(k,:,2)=-ef(k,:,1)'.*(cov(:,2).*cov(:,6)-cov(:,3).*cov(:,5)) ...
               +ef(k,:,2)'.*(cov(:,1).*cov(:,6)-cov(:,3).^2) ...
               -ef(k,:,3)'.*(cov(:,1).*cov(:,5)-cov(:,2).*cov(:,3));
            cov_inv_ef(k,:,3)=ef(k,:,1)'.*(cov(:,2).*cov(:,5)-cov(:,3).*cov(:,4)) ...
               -ef(k,:,2)'.*(cov(:,1).*cov(:,5)-cov(:,3).*cov(:,2)) ...
               +ef(k,:,3)'.*(cov(:,1).*cov(:,4)-cov(:,2).^2);
         end
      end
   end
   
   for k=1:numcontrasts
      if which_stats(k,1)
         T2=sum(ef(k,:,:).*cov_inv_ef(k,:,:),3)'./(det+(det<=0)).*(det>0) ...
            /covef(k,k);
         out.file_name=[deblank(output_file_base(k,:)) '_T2.' ext];
         out.dim=[d.dim(1:3) 0];
         out.data=reshape(T2,d.dim(1),d.dim(2));
         fmris_write_image(out,slice,1);
      end
   end
   
   if which_stats(1,4)
      a=zeros(d.dim(1)*d.dim(2),L,L);
      for k=1:L
         covef_inv_ef=covef_inv*ef(:,:,k);
         for j=1:L
            a(:,j,k)=sum(cov_inv_ef(:,:,j).*covef_inv_ef,1)';
         end
      end
      if L==1
         roy=a(:,1,1)./(det+(det<=0)).*(det>0)/p;
      end
      if L==2
         a0=a(:,1,1).*a(:,2,2)-a(:,1,2).*a(:,2,1);
         a1=(a(:,1,1)+a(:,2,2))/2;
         roy=(a1+abs(real(sqrt(a1.^2-a0))))./(det+(det<=0)).*(det>0)/p;
      end
      if L==3
         a0=-a(:,1,1).*(a(:,2,2).*a(:,3,3)-a(:,2,3).*a(:,3,2)) ...
            +a(:,1,2).*(a(:,2,1).*a(:,3,3)-a(:,2,3).*a(:,3,1)) ...
            -a(:,1,3).*(a(:,2,1).*a(:,3,2)-a(:,2,2).*a(:,3,1));
         a1=a(:,2,2).*a(:,3,3)-a(:,2,3).*a(:,3,2) ...
            +a(:,1,1).*a(:,3,3)-a(:,1,3).*a(:,3,1) ... 
            +a(:,1,1).*a(:,2,2)-a(:,1,2).*a(:,2,1);
         a2=-(a(:,1,1)+a(:,2,2)+a(:,3,3));
         q=a1/3-a2.^2/9;
         r=(a1.*a2-3*a0)/6-a2.^3/27;
         s1=(r+sqrt(q.^3+r.^2)).^(1/3);
         z=zeros(d.dim(1)*d.dim(2),L);
         z(:,1)=2*real(s1)-a2/3;
         z(:,2)=-real(s1)-a2/3+sqrt(3)*imag(s1);
         z(:,3)=-real(s1)-a2/3-sqrt(3)*imag(s1);
         roy=max(z,[],2)./(det+(det<=0)).*(det>0)/p;
      end
      out.file_name=[deblank(output_file_base(1,:)) '_Roy.' ext];
      out.dim=[d.dim(1:3) 0];
      out.data=reshape(roy,d.dim(1),d.dim(2));
      fmris_write_image(out,slice,1);
   end

   if which_stats(1,9)
      % Finds an estimate of the fwhm for each of the 8 cube corners surrounding
      % a voxel, then averages. Results for components are averaged.
      
      for i=1:n
         for j=1:L
            vv=cov(:,diags(j))*df;
            wresid_slice(:,i,j)=Y(i,:,j)'.*(vv>0)./sqrt(vv+(vv<=0));
         end
      end
      
      if slice==1
         u=reshape(wresid_slice,I,J,n,L);
         ux=diff(u,1,1);
         uy=diff(u,1,2);
         Axx=squeeze(sum(ux.^2,3));
         Ayy=squeeze(sum(uy.^2,3));
         if D==2
            for index=1:4
               i=(1:Im)+ip(index);
               j=(1:Jm)+jp(index);
               axx=Axx(:,j,:);
               ayy=Ayy(i,:,:);
               if df_resid>df_limit
                  axy=squeeze(sum(ux(:,j,:,:).*uy(i,:,:,:),3))*is(index)*js(index);
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
               f(i,j,:)=f(i,j,:)+(detlamf>0).*(detlamf+(detlamf<=0)).^alphaf;
               r(i,j,:)=r(i,j,:)+(detlamr>0).*(detlamr+(detlamr<=0)).^alphar;
            end
         end
      else 
         uz=reshape(wresid_slice,I,J,n,L)-u;
         Azz=squeeze(sum(uz.^2,3));
         % The 4 upper cube corners:
         for index=1:4
            i=(1:Im)+ip(index);
            j=(1:Jm)+jp(index);
            axx=Axx(:,j,:);
            ayy=Ayy(i,:,:);
            azz=Azz(i,j,:);
            if df_resid>df_limit
               axy=squeeze(sum(ux(:,j,:,:).*uy(i,:,:,:),3))*is(index)*js(index);
               axz=squeeze(sum(ux(:,j,:,:).*uz(i,j,:,:),3))*is(index);
               ayz=squeeze(sum(uy(i,:,:,:).*uz(i,j,:,:),3))*js(index);
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
            f(i,j,:)=f(i,j,:)+(detlamf>0).*(detlamf+(detlamf<=0)).^alphaf;
            r(i,j,:)=r(i,j,:)+(detlamr>0).*(detlamr+(detlamr<=0)).^alphar;
         end
         f=consf/((slice>2)+1)*f./nxy;
         r=consr/((slice>2)+1)*r./nxy;
         out.file_name=[deblank(output_file_base(1,:)) '_fwhm.' ext];
         out.dim=[d.dim(1:3) 2];
         ff=mean(f,3);
         out.data=ff.*(ff<100)+100*(ff>=100);
         fmris_write_image(out,slice-1,1);
         out.data=mean(r,3);
         fmris_write_image(out,slice-1,2);
         
         f=zeros(I,J,L);
         r=zeros(I,J,L);
         u=reshape(wresid_slice,I,J,n,L);
         ux=diff(u,1,1);
         uy=diff(u,1,2);
         Axx=squeeze(sum(ux.^2,3));
         Ayy=squeeze(sum(uy.^2,3));
         % The 4 lower cube corners:
         for index=1:4
            i=(1:Im)+ip(index);
            j=(1:Jm)+jp(index);
            axx=Axx(:,j,:);
            ayy=Ayy(i,:,:);
            azz=Azz(i,j,:);
            if df_resid>df_limit
               axy= squeeze(sum(ux(:,j,:,:).*uy(i,:,:,:),3))*is(index)*js(index);
               axz=-squeeze(sum(ux(:,j,:,:).*uz(i,j,:,:),3))*is(index);
               ayz=-squeeze(sum(uy(i,:,:,:).*uz(i,j,:,:),3))*js(index);
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
            f(i,j,:)=f(i,j,:)+(detlamf>0).*(detlamf+(detlamf<=0)).^alphaf;
            r(i,j,:)=r(i,j,:)+(detlamr>0).*(detlamr+(detlamr<=0)).^alphar;
         end
      end
      if slice==d.dim(3)
         f=consf*f./nxy;
         r=consr*r./nxy;
         ff=mean(f,3);
         out.data=ff.*(ff<100)+100*(ff>=100);
         fmris_write_image(out,slice,1);
         out.data=mean(r,3);
         fmris_write_image(out,slice,2);
      end
   end
end  

return




