% Looking at the fMRI data using pca_image

input_file='c:/keith/data/brian_ha_19971124_1_100326_mri_MC.mnc';
mask_file='c:/keith/data/brian_ha_19971124_1_100326_mri_MC.mnc';
mask_thresh=fmri_mask_thresh(mask_file);
pca_image(input_file, [], 4, mask_file, mask_thresh);
saveas(gcf,'c:/keith/fmristat/figs_test/figpca1.jpg');

exclude=[1 2 3];
pca_image(input_file, exclude, 4, mask_file);
saveas(gcf,'c:/keith/fmristat/figs_test/figpca2.jpg');

X_remove=[ones(120,1) (1:120)'];
X_interest=repmat(eye(12),10,1);
pca_image(input_file, exclude, 4, mask_file, [], [], X_remove, X_interest);
saveas(gcf,'c:/keith/fmristat/figs_test/figpca3.jpg');

% Plotting the hemodynamic response function (hrf) using fmridesign

hrf_parameters=[5.4 5.2 10.8 7.35 0.35]
time=(0:240)/10;
hrf0=fmridesign(time,0,[1 0],[],hrf_parameters);
clf;
plot(time,squeeze(hrf0.X(:,1,1)),'LineWidth',2)
xlabel('time (seconds)')
ylabel('hrf')
title('Glover hrf model')
saveas(gcf,'c:/keith/fmristat/figs_test/fighrf0.jpg');

% Making the design matrices using fmridesign

frametimes=(0:119)*3;
slicetimes=[0.14 0.98 0.26 1.10 0.38 1.22 0.50 1.34 0.62 1.46 0.74 1.58 0.86];
eventid=repmat([1; 2],10,1);
eventimes=(0:19)'*18+9;
duration=ones(20,1)*9;
height=ones(20,1);
events=[eventid eventimes duration height] 
X_cache=fmridesign(frametimes,slicetimes,events,[],hrf_parameters);

S=kron(ones(10,1),kron([0 0; 1 0; 0 0; 0 1],ones(3,1)));
X_cache=fmridesign(frametimes,slicetimes,  []  , S,hrf_parameters);

X_cache=fmridesign(frametimes,slicetimes,events);

plot(squeeze(X_cache.X(:,:,1,4)),'LineWidth',2)
legend('Hot','Warm')
xlabel('frame number')
ylabel('response')
title('Hot and Warm responses')
saveas(gcf,'c:/keith/fmristat/figs_test/figdesign.jpg');

% Analysing one run with fmrilm

contrast=[1  0;
          0  1;
          1 -1];
exclude=[1 2 3];
which_stats=[1 1 1 1 1 0 0 0 1];

input_file='c:/keith/data/brian_ha_19971124_1_100326_mri_MC.mnc';
output_file_base=['c:/keith/results_test/ha_100326_hot';
                  'c:/keith/results_test/ha_100326_wrm';
                  'c:/keith/results_test/ha_100326_hmw']

[df1 p]=fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats);

% Visualizing the results using view_slices, glass_brain and blob_brain

t_file='c:/keith/results/ha_100326_hmw_mag_t.mnc';
m=fmris_read_image(t_file,4,1);
imagesc(m.data',[-6 6]); colorbar; axis xy; colormap(spectral);
saveas(gcf,'c:/keith/fmristat/figs_test/fmristat/figslice.jpg');

mask_file=input_file;
clf;
view_slices(t_file,mask_file,[],3,1,[-6 6]);
saveas(gcf,'c:/keith/fmristat/figs_test/figviewslice.jpg');

clf;
view_slices(t_file,mask_file,[],0:11,1,[-6 6]);
saveas(gcf,'c:/keith/fmristat/figs_test/figviewslices.jpg');

glass_brain(t_file,3,mask_file);
saveas(gcf,'c:/keith/fmristat/figs_test/figlassbrain.jpg');

clf;
blob_brain(t_file,5,'c:/keith/results_test/ha_100326_hmw_mag_ef.mnc');
title('Hot - warm, T>5, coloured by effect (%BOLD)');
saveas(gcf,'c:/keith/fmristat/figs_test/figblobrain.jpg');

% F-tests

contrast=[0  0  0 1 0 0;
          0  0  0 0 1 0;
          0  0  0 0 0 1]
which_stats=[0 0 0 1];
output_file_base='c:/keith/results_test/ha_100326_drift';
fwhm_rho='c:/keith/results_test/ha_100326_hot_rho.mnc';
fmrilm(input_file,output_file_base,X_cache,contrast,exclude,which_stats,fwhm_rho);
clf;
view_slices('c:/keith/results_test/ha_100326_drift_mag_F.mnc',mask_file,[],0:11,1,[0 50]);
stat_threshold(1000000,26000,8,[3 110])
saveas(gcf,'c:/keith/fmristat/figs_test/figdrift.jpg');

% A linear effect of temperature
 
temperature=[45.5 35.0 49.0 38.5 42.0 49.0 35.0 42.0 38.5 45.5 ...
             38.5 49.0 35.0 45.5 42.0 45.5 38.5 42.0 35.0 49.0]';
       
Temperature=35:3.5:49;
subplot(2,2,1)
plot(Temperature,1:5,'LineWidth',2);
xlabel('Temperature'); ylabel('Response');
title('(a) Linear temperature effect');
subplot(2,2,2)
plot(Temperature,10-((1:5)-4).^2,'LineWidth',2);
xlabel('Temperature'); ylabel('Response');
title('(b) Quadratic temperature effect');
subplot(2,2,3)
plot(Temperature,10-((1:5)-4).^2+((1:5)-3).^3,'LineWidth',2);
xlabel('Temperature'); ylabel('Response');
title('(c) Cubic temperature effect');       
subplot(2,2,4)
plot(Temperature,[1 4 3 5 3.5],'LineWidth',2);
xlabel('Temperature'); ylabel('Response');
title(['(d) Quartic = arbitrary temperature effect']);       
saveas(gcf,'c:/keith/fmristat/figs_test/figpoly.jpg');

events=[zeros(20,1)+1 eventimes duration ones(20,1);
        zeros(20,1)+2 eventimes duration temperature] 

contrast=[0  1;
          1 49;
          1 35;
          0 14];

events=[zeros(20,1)+1 eventimes duration ones(20,1);
        zeros(20,1)+2 eventimes duration temperature; 
        zeros(20,1)+3 eventimes duration temperature.^2] 
contrast=[0 0 1];

contrast=[1  0  0  0  0;
          0  1  0  0  0;
          0  0  1  0  0;
          0  0  0  1  0;
          0  0  0  0  1];
contrast=[eye(5) ones(5,4)];
which_stats=[0 0 0 1];

contrast=[.8 -.2 -.2 -.2 -.2;
         -.2  .8 -.2 -.2 -.2;
         -.2 -.2  .8 -.2 -.2;
         -.2 -.2 -.2  .8 -.2;
         -.2 -.2 -.2 -.2  .8];
contrast=[eye(5)-ones(5)/5 ones(5,4)];
which_stats=[0 0 0 1];

contrast=[-7.0 -3.5  0  3.5 7.0];
which_stats=[1 1 1];

% Combining runs/sessions/subjects with multistat

contrast=[1 -1];
which_stats=[1 1 1];
      
input_file='c:/keith/data/brian_ha_19971124_1_093923_mri_MC.mnc';
output_file_base='c:/keith/results_test/ha_093923_hmw';
[df2 p]=fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats);

input_file='c:/keith/data/brian_ha_19971124_1_101410_mri_MC.mnc';
output_file_base='c:/keith/results_test/ha_101410_hmw';
[df3 p]=fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats);
 
input_file='c:/keith/data/brian_ha_19971124_1_102703_mri_MC.mnc';
output_file_base='c:/keith/results_test/ha_102703_hmw';
[df4 p]=fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats);

X=[1 1 1 1]';
contrast=[1];
which_stats=[1 1 1];
input_files_ef=['c:/keith/results_test/ha_093923_hmw_mag_ef.mnc';
                'c:/keith/results_test/ha_100326_hmw_mag_ef.mnc';
                'c:/keith/results_test/ha_101410_hmw_mag_ef.mnc';
                'c:/keith/results_test/ha_102703_hmw_mag_ef.mnc'];
input_files_sd=['c:/keith/results_test/ha_093923_hmw_mag_sd.mnc';
                'c:/keith/results_test/ha_100326_hmw_mag_sd.mnc';
                'c:/keith/results_test/ha_101410_hmw_mag_sd.mnc';
                'c:/keith/results_test/ha_102703_hmw_mag_sd.mnc'];
output_file_base='c:/keith/results_test/ha_multi_hmw_fixed'
input_files_df=[df1 df2 df3 df4]

clf;
view_slices('c:/keith/results_test/ha_100326_hot_fwhm.mnc',mask_file,0,3,1,[0 20]);
saveas(gcf,'c:/keith/fmristat/figs_test/figfwhm.jpg');
input_files_fwhm=8;

df=multistat(input_files_ef,input_files_sd,input_files_df, ...
   input_files_fwhm,X,contrast,output_file_base,which_stats,Inf)

output_file_base='c:/keith/results_test/ha_multi_hmw_random'
df=multistat(input_files_ef,input_files_sd,input_files_df, ...
   input_files_fwhm,X,contrast,output_file_base,which_stats,0)

% Fixed and random effects

which_stats=[1 1 1 1 0 0 0 0 1];
output_file_base='c:/keith/results_test/ha_multi_hmw'
[df, df_resid]=multistat(input_files_ef,input_files_sd,input_files_df, ...
   input_files_fwhm,X,contrast,output_file_base,which_stats)

clf;
view_slices('c:/keith/results_test/ha_multi_hmw_sdratio.mnc',mask_file,0,3,1);
saveas(gcf,'c:/keith/fmristat/figs_test/figsdratio.jpg');

clf;
view_slices('c:/keith/results_test/ha_multi_hmw_fwhm.mnc',mask_file,0,3,1,[0 20]);
saveas(gcf,'c:/keith/fmristat/figs_test/figfwhmulti.jpg');

blured_fwhm=gauss_blur('c:/keith/results_test/ha_multi_hmw_fwhm.mnc',10);
clf;
view_slices(blured_fwhm,mask_file,0,3,1,[0 20]);
saveas(gcf,'c:/keith/fmristat/figs_test/figfwhmultiblur.jpg');

% Thresholding the tstat image with stat_threshold and fdr_threshold

stat_threshold(1000000,26000,8,100)
mask_file='c:/keith/data/brian_ha_19971124_1_100326_mri_MC.mnc';
[search_volume, num_voxels]=mask_vol(mask_file)
stat_threshold(search_volume,num_voxels,8,100)
stat_threshold(search_volume,num_voxels,8,[100; 100])
stat_threshold(search_volume,num_voxels,8,[100 0; 3 100])
t_file='c:/keith/results_test/ha_multi_hmw_t.mnc';
fdr_threshold(t_file,[],mask_file,[],100)

% Finding the exact resels of a search region

mask_file='c:/keith/data/brian_ha_19971124_1_100326_mri_MC.mnc';
resels=mask_resels(8,[],mask_file)
stat_threshold(resels,num_voxels,1,100)
stat_threshold(resels,Inf,1,100)
stat_threshold(search_volume,Inf,8,100)

which_stats=[0 0 0 0 0 0 1];
contrast=[];
output_file_base='c:/keith/results_test/ha_100326_hot';
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats);
fwhm_info='c:/keith/results_test/ha_100326_hot_wresid.mnc';
resels=mask_resels(fwhm_info,[],mask_file)
stat_threshold(resels,num_voxels,1,[100; 100])

which_stats=[0 0 0 0 0 0 1];
contrast=[];
output_file_base='c:/keith/results_test/ha_multi_hmw';
[df, df_resid]=multistat(input_files_ef,input_files_sd,input_files_df, ...
   input_files_fwhm,X,contrast,output_file_base,which_stats);
fwhm_info='c:/keith/results_test/ha_multi_hmw_wresid.mnc';
resels=mask_resels(fwhm_info,[df_resid df],mask_file)
stat_threshold(resels,num_voxels,1,[100 0; 3 100])
    
% Curvature (not on web page):

output_file_base='c:/keith/results_test/ha_100326_hot_wresid';    
mesh_tet([output_file_base '.mnc'], output_file_base);   
tet_curv([output_file_base '_tet.mnc'], output_file_base);
gauss_blur([output_file_base '_curv.mnc'],20);
clf;
view_slices([output_file_base '_curv_blur20.mnc'], ...
   mask_file,0,0:11,1,[-0.03 0.03]);

output_file_base='c:/keith/results_test/ha_multi_hmw_wresid';    
mesh_tet([output_file_base '.mnc'], output_file_base);   
tet_curv([output_file_base '_tet.mnc'], output_file_base, [3 110]);
gauss_blur([output_file_base '_curv.mnc'],20);
clf;
view_slices([output_file_base '_curv_blur20.mnc'], ...
   mask_file,0,0:11,1,[-0.03 0.03]);

% Locating peaks and clusters with locmax

lm=locmax(t_file, 4.93, mask_file)
pval=stat_threshold(1000000,26000,8,100,lm(:,1))

% Producing an SPM style summary with stat_summary

t_file='c:/keith/results_test/ha_multi_hmw_t.mnc';
fwhm_file='c:/keith/results_test/ha_multi_hmw_fwhm.mnc';
stat_summary(t_file, fwhm_file, [100 0; 3 100], mask_file);
saveas(gcf,'c:/keith/fmristat/figs_test/figlassbrainmulti.jpg');

% Confidence regions for the spatial location of local maxima using conf_region

conf_region(t_file,4.93,mask_file)
conf_file='c:/keith/results_test/ha_multi_hmw_t_95cr.mnc';
clf;
blob_brain(conf_file,9,conf_file,9)
title('Approx 95% confidence regions for peak location, coloured by peak height')
saveas(gcf,'c:/keith/fmristat/figs_test/figconfregion.jpg');

% Conjunctions

which_stats=[0 0 0 0 1];
output_file_base='c:/keith/results_test/ha_multi_hmw'
multistat(input_files_ef,input_files_sd,input_files_df, ...
   input_files_fwhm,X,contrast,output_file_base,which_stats,20)
stat_summary('c:/keith/results_test/ha_multi_hmw_conj.mnc', fwhm_file, ...
    [100 0; 3 100], mask_file, [], 0.001, 1, 4);
saveas(gcf,'c:/keith/fmristat/figs_test/figconj.jpg');
stat_threshold(search_volume,num_voxels,8,100,0.05,0.001,0.05,4)

% Extracting values from a minc file using extract

voxel=[62  69   3]
ef=extract(voxel,'c:/keith/results_test/ha_multi_hmw_ef.mnc')  
sd=extract(voxel,'c:/keith/results_test/ha_multi_hmw_sd.mnc')
ef/sd

[df,p,spatial_av]=fmrilm(input_file,[],[],[],exclude);
ref_data=squeeze(extract(voxel,input_file))./spatial_av*100;
ef_hot=extract(voxel,'c:/keith/results_test/ha_100326_hot_mag_ef.mnc')
ef_wrm=extract(voxel,'c:/keith/results_test/ha_100326_wrm_mag_ef.mnc')
fitted=mean(ref_data)+ef_hot*X_cache.X(:,1,1,voxel(3)+1) ...
                     +ef_wrm*X_cache.X(:,2,1,voxel(3)+1);
clf;
plot(frametimes,[ref_data fitted],'LineWidth',2); 
legend('Reference data','Fitted values');
xlabel('time (seconds)');
ylabel('fMRI response, percent');
title(['Observed (reference) and fitted data, ignoring trends, at voxel ' num2str(voxel)]);
saveas(gcf,'c:/keith/fmristat/figs_test/figfit.jpg');

% Estimating the time course of the response

eventid=kron(ones(10,1),(1:12)');
eventimes=frametimes';
duration=ones(120,1)*3;
height=ones(120,1);
events=[eventid eventimes duration height]
X_bases=fmridesign(frametimes,slicetimes,events,[],zeros(1,5));
contrast=[eye(12)-ones(12)/12];
num2str(round(contrast*100)/100)
exclude=[1 2 3];
which_stats=[0 1 1 1];
input_file='c:/keith/data/brian_ha_19971124_1_100326_mri_MC.mnc';
output_file_base=['c:/keith/results_test/ha_100326_time01';
'c:/keith/results_test/ha_100326_time02';
'c:/keith/results_test/ha_100326_time03';
'c:/keith/results_test/ha_100326_time04';
'c:/keith/results_test/ha_100326_time05';
'c:/keith/results_test/ha_100326_time06';
'c:/keith/results_test/ha_100326_time07';
'c:/keith/results_test/ha_100326_time08';
'c:/keith/results_test/ha_100326_time09';
'c:/keith/results_test/ha_100326_time10';
'c:/keith/results_test/ha_100326_time11';
'c:/keith/results_test/ha_100326_time12'];
[df, p]=fmrilm(input_file,output_file_base,X_bases,contrast,exclude,which_stats,fwhm_rho)

stat_threshold(search_volume,num_voxels,8,[p df])
lm=locmax([output_file_base(1,:) '_mag_F.mnc'],5.2);
num2str(lm)

voxel=[62 69 3]
values=extract(voxel,output_file_base,'_mag_ef.mnc')
sd=extract(voxel,output_file_base,'_mag_sd.mnc')

b_hot=extract(voxel,'c:/keith/results_test/ha_100326_hot_mag_ef.mnc')
b_wrm=extract(voxel,'c:/keith/results_test/ha_100326_wrm_mag_ef.mnc')
time=(1:360)/10;
X_hrf=fmridesign(time,0,[1 9 9 1]);
hrf=squeeze(X_hrf.X(:,1,1));
plot((0:12)*3+slicetimes(voxel(3)+1),values([1:12 1]),'k', ...
   [0:11; 0:11]*3+slicetimes(voxel(3)+1), [values+sd; values-sd],'g', ...
   time,[zeros(1,90) ones(1,90) zeros(1,180)],'r', ...
   time,[zeros(1,270) ones(1,90)],'b', ...
   time,hrf*b_hot+hrf([181:360 1:180])*b_wrm,'g','LineWidth',2);
legend('Estimated response','Modeled response');
text(10,0.5,'Hot')
text(28,0.5,'Warm')
xlabel('time (seconds) from start of epoch');
ylabel('fMRI response, percent');
title(['Estimated and modeled response at voxel ' num2str(voxel)]);
saveas(gcf,'c:/keith/fmristat/figs_test/figmodelresp.jpg');

% Estimating the delay 

contrast=[1  0;
          0  1;
          1 -1]
which_stats=[1 1 1]
input_file='c:/keith/data/brian_ha_19971124_1_100326_mri_MC.mnc'; 
output_file_base=['c:/keith/results_test/ha_100326_hot';
                  'c:/keith/results_test/ha_100326_wrm';
                  'c:/keith/results_test/ha_100326_hmw']
n_trends=[];
fwhm_rho='c:/keith/results_test/ha_100326_hot_rho.mnc';
confounds=[];
contrast_is_delay=[1 1 1];
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, ...
    which_stats, fwhm_rho, n_trends, confounds, contrast_is_delay);

slice=2; 
subplot(2,2,1);
view_slices('c:/keith/results_test/ha_100326_hot_mag_t.mnc',mask_file,[],slice,1,[-6 6]);
subplot(2,2,2);
view_slices('c:/keith/results_test/ha_100326_hot_del_ef.mnc',mask_file,[],slice,1,[-3 3]);
subplot(2,2,3);
view_slices('c:/keith/results_test/ha_100326_hot_del_sd.mnc',mask_file,[],slice,1,[0 6]);
subplot(2,2,4);
view_slices('c:/keith/results_test/ha_100326_hot_del_t.mnc',mask_file,[],slice,1,[-6 6]);
saveas(gcf,'c:/keith/fmristat/figs_test/figdelay.jpg');

clf;
blob_brain('c:/keith/results_test/ha_100326_hot_mag_t.mnc',5, ...
   'c:/keith/results_test/ha_100326_hot_del_ef.mnc');
title('Delay (secs) of hot stimulus where T > 5')
saveas(gcf,'c:/keith/fmristat/figs_test/figdelay3D.jpg');

% Efficiency and choosing the best design 

efficiency(X_cache, contrast, exclude)

[sd_ef, Y]=efficiency(X_cache, contrast, exclude);
slice=4;
plot(squeeze(X_cache.X(:,:,1,slice)),'LineWidth',2)
hold on; plot(Y(:,:,slice)'/max(max(Y(:,:,slice))),':','LineWidth',2); hold off
xlim([0 40])
legend('Hot resp','Warm resp','Hot wt','Warm wt','Hot - Warm wt')
xlabel('frame number')
ylabel('response')
title('Hot and Warm responses and weights')
saveas(gcf,'c:/keith/fmristat/figs_test/figrespwt.jpg');

events=[1 90 90 1; 2 270 90 1] 
X_cache1=fmridesign(frametimes,slicetimes,events);
[sd_ef, Y]=efficiency(X_cache1, contrast, exclude);
sd_ef
slice=4;
plot(squeeze(X_cache1.X(:,:,1,slice)),'LineWidth',2)
hold on; plot(Y(:,:,4)'/max(max(Y(:,:,slice))),':','LineWidth',2); hold off
legend('Hot resp','Warm resp','Hot wt','Warm wt','Hot - Warm wt')
xlabel('frame number')
ylabel('response')
title('Hot and Warm responses and weights')
saveas(gcf,'c:/keith/fmristat/figs_test/figrespwt1.jpg');

rho=0.3;
n_temporal=3;
confounds=[];
contrast_is_delay=[1 1 1];
efficiency(X_cache, contrast, exclude, rho, n_temporal, confounds, contrast_is_delay)

%%%%%%%%%%%%%%%%%%%% Extra code for the optimum design figure

slicetimes1=mean(slicetimes);
SI=1:20;
ISI=0:20;
sd_ef_mag=zeros(3,length(SI),length(ISI));
sd_ef_del=zeros(3,length(SI),length(ISI));
for i=1:length(SI)
   si=SI(i)  
   for j=1:length(ISI)
      isi=ISI(j);
      nb=ceil(length(frametimes)*(frametimes(2)-frametimes(1))/(si+isi)/2);
      eventid=kron(ones(nb,1),[1; 2]);
      eventimes=((1:2*nb)-1)'*(si+isi)+isi;
      duration=ones(2*nb,1)*si;
      height=ones(2*nb,1);
      events=[eventid eventimes duration height];
      X_cache1=fmridesign(frametimes,slicetimes1,events);
      sd_ef_mag(:,i,j)=efficiency(X_cache1, contrast, exclude, rho);
      mag_t=1.0./sd_ef_mag(1:2,i,j)';
      sd_ef_del(:,i,j)=efficiency(X_cache1, contrast, exclude, rho, n_poly, ...
   confounds, contrast_is_delay, mag_t);
   end
end

clf;
%whitebg
range=[0 0.5];
subplot(2,2,1)
imagesc(SI,ISI,(squeeze(sd_ef_mag(1,:,:))'),range); axis xy; colorbar; 
xlabel('Stimulus Duration (sec)')
ylabel('InterStimulus Interval (sec)')
title('(a) Sd of magnitude of hot stimulus')
subplot(2,2,2)
imagesc(SI,ISI,(squeeze(sd_ef_mag(3,:,:))'),range); axis xy; colorbar; 
xlabel('Stimulus Duration (sec)')
ylabel('InterStimulus Interval (sec)')
title('(b) Sd of magnitude of hot-warm')
range=[0 1];
subplot(2,2,3)
m=sd_ef_del(1,:,:).*(sd_ef_mag(1,:,:)<=0.25); 
imagesc(SI,SI,(squeeze(m)'),range); axis xy; colorbar; 
xlabel('Stimulus Duration (sec)')
ylabel('InterStimulus Interval (sec)')
title('(c) Sd of delay of hot stimulus (sec)')
subplot(2,2,4)
m=sd_ef_del(3,:,:).*(sd_ef_mag(1,:,:)<=0.25);
imagesc(SI,ISI,(squeeze(m)'),range); axis xy; colorbar; 
xlabel('Stimulus Duration (sec)')
ylabel('InterStimulus Interval (sec)')
title('(d) Sd of delay of hot-warm (sec)')
colormap(spectral);
saveas(gcf,'c:/keith/fmristat/figs_test/figsd_ef.jpg');

nevents=round(360./(5:20));
tevents=360./nevents
contrast=[1];
rho=0.3;
n_poly=3;
confounds=[];
contrast_is_delay=[1];
sd_ef=zeros(length(nevents),6);
for i=1:length(nevents)
   n=nevents(i);
   events=[ones(n,1) ((1:n)'-0.5)/n*359 zeros(n,1) ones(n,1)];
   X_cache1=fmridesign(frametimes,slicetimes1,events);
   sd_ef(i,1)=efficiency(X_cache1, contrast, exclude, rho);
   mag_t=5/sd_ef(i,1);
   sd_ef(i,2)=efficiency(X_cache1, contrast, exclude, rho, n_poly, ...
      confounds, contrast_is_delay, mag_t);
   events=[ones(n,1) rand(n,1)*359 zeros(n,1) ones(n,1)];
   X_cache1=fmridesign(frametimes,slicetimes1,events);
   sd_ef(i,3)=efficiency(X_cache1, contrast, exclude, rho);
   mag_t=5/sd_ef(i,3);
   sd_ef(i,4)=efficiency(X_cache1, contrast, exclude, rho, n_poly, ...
      confounds, contrast_is_delay, mag_t);
   events=[1 180 0 n];
   X_cache1=fmridesign(frametimes,slicetimes1,events);
   sd_ef(i,5)=efficiency(X_cache1, contrast, exclude, rho);
   mag_t=5/sd_ef(i,5);
   sd_ef(i,6)=efficiency(X_cache1, contrast, exclude, rho, n_poly, ...
      confounds, contrast_is_delay, mag_t);
end

clf;
%whitebg
plot(tevents,sd_ef(:,[1 3 5])/2,'LineWidth',2)
hold on; plot(tevents,sd_ef(:,[2 4 6]),':','LineWidth',2); hold off;
ylim([0 1]/2);
legend('uniform . . . . . . . . .','random ..  .  ...  .. .','concentrated   :')
xlabel('Average time between events (secs)')
ylabel('Sd of effect (secs for delays)')
title('Efficiency of magnitudes (solid, x 0.5) and delays (dotted)')
saveas(gcf,'c:/keith/fmristat/figs_test/figsd_ef_event.jpg');

% Effective connectivity of all voxels with a reference voxel 

output_file_base='c:/keith/results_test/ha_100326_connect';
contrast=[0  0  0  0  0  0  0  1];
which_stats=[1 1 1];
voxel=[62 69 3];
ref_times=frametimes'+slicetimes(voxel(3)+1);
[df,p,spatial_av]=fmrilm(input_file,[],[],[],exclude);
ref_data=squeeze(extract(voxel,input_file))./spatial_av*100;
confounds=fmri_interp(ref_times,ref_data,frametimes,slicetimes);
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, ...
    which_stats, fwhm_rho, [3 1 1 spatial_av'], confounds);
clf;
view_slices('c:/keith/results_test/ha_100326_connect_mag_t.mnc',mask_file,[],0:11,1,[-6 6]);
saveas(gcf,'c:/keith/fmristat/figs_test/figconnect.jpg');

% Higher order autoregressive models 

clf;
view_slices('c:/keith/results_test/ha_100326_hot_rho.mnc',mask_file,0,3,1,[-0.15 0.35]); 
saveas(gcf,'c:/keith/fmristat/figs_test/figrho.jpg');

which_stats=[1 1 1 0 1 0 0 1];
contrast=[1 -1];
output_file_base='c:/keith/results_test/ha_100326_hmw_ar4';
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, ...
    which_stats, [], [], [], [], [], [], 4);
clf;
view_slices('c:/keith/results_test/ha_100326_hmw_ar4_A.mnc',mask_file,0,3,1:4,[-0.15 0.35]); 
saveas(gcf,'c:/keith/fmristat/figs_test/figar4.jpg');

% vbm using fmristat

fwhm_file='c:/keith/fmristat/vbm/TLE_R_vs_NC_lin_GM_10mm_fwhm.mnc';
mask_file='c:/keith/fmristat/vbm/TLE_R_vs_NC_lin_GM_10mm_mask.mnc';
clf;
view_slices(fwhm_file, mask_file, [0.05 0.95], 90)
saveas(gcf,'c:/keith/fmristat/figs/figvbmfwhm.jpg');

t_file='c:/keith/fmristat/vbm/TLE_R_vs_NC_lin_GM_10mm_tstat.mnc';
[c,p]=stat_summary(t_file, fwhm_file, [83; 83], mask_file, [0.05 0.95], 0.001, -1);
saveas(gcf,'c:/keith/fmristat/figs/figvbmt.jpg');



 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vol2exp commands:

vol2exp('c:/keith/results_test/ha_093923_hmw_mag_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_100326_hmw_mag_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_101410_hmw_mag_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_102703_hmw_mag_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_multi_hmw_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_multi_hmw_conj.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_multi_hot_del_ef.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_multi_hot_del_sd.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_multi_hot_del_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_multi_hot_mag_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_100326_hot_del_ef.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_100326_hot_del_sd.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_100326_hot_del_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/ha_100326_hot_mag_t.mnc',[34 1 95; 22 1 92; 1 1 13]);
vol2exp('c:/keith/results_test/random_mag_t.mnc',[34 1 95; 22 1 92; 1 1 13]);


vol2exp('c:/keith/results_test/ha_100326_pca.mnc',[34 1 95; 22 1 92; 1 1 13], ...
   'c:/keith/results_test/ha_100326_pca1',1);
vol2exp('c:/keith/results_test/ha_100326_pca.mnc',[34 1 95; 22 1 92; 1 1 13], ...
   'c:/keith/results_test/ha_100326_pca2',2);
vol2exp('c:/keith/results_test/ha_100326_pca.mnc',[34 1 95; 22 1 92; 1 1 13], ...
   'c:/keith/results_test/ha_100326_pca3',3);

vol2exp('c:/keith/results_test/ha_100326_hot_mag_t.mnc');
vol2exp('c:/keith/results_test/ha_100326_hot_del_ef.mnc');
!explorer -map c:/keith/fmristat_old/delay.map




