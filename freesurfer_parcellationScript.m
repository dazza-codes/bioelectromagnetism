clear

if isunix,
    nri_dir = '/data/acetylcholine2/nri/';
else
    nri_dir = 'e:\nri\';
end
subjectsPath = nri_dir;

if exist(subjectsPath) ~= 7,
    msg = sprintf('cannot find subjects directory:\n%s',subjectsPath);
    error(msg)
end

nri_demographics

% just for initial testing
subjects.id = subjects.id(1:15);
subjects.group = subjects.group(1:15);
subjects.age = subjects.age(1:15);

for subjectN = 1:length(subjects.id),

    subject = lower(['nri_',subjects.id{subjectN}]);
    group = subjects.group{subjectN};
    gender = subjects.gender{subjectN};
    age = subjects.age{subjectN};

    subjectPath = fullfile(subjectsPath,subject,'');
    labelPath = fullfile(subjectPath,'label','');
    avwPath = fullfile(subjectPath,'mri','analyze','');

    % load freesurfer surfaces
    subjectFS = freesurfer_read_subject(subjectPath);

    % histogram of thickness
    bins = 0:0.1:6;
    subplot(2,1,1); hist(subjectFS.lh.thickness,bins); ylabel('LH')
    title(upper(strrep(subject,'nri_','')))
    subplot(2,1,2); hist(subjectFS.rh.thickness,bins); ylabel('RH')
    xlabel('Cortical Thickness (mm)')
    imagefile = fullfile(subjectPath,[subject,'.histogram.thickness.png']);
    save_png(imagefile,gcf,300)
    close(gcf)

    % evaluate thickness by curvature
    lh.GyrusIndex = find(subjectFS.lh.curv < 0);
    lh.SulcusIndex = find(subjectFS.lh.curv > 0);
    rh.GyrusIndex = find(subjectFS.rh.curv < 0);
    rh.SulcusIndex = find(subjectFS.rh.curv > 0);

    % extract parcellation values
    aparc = unique(subjectFS.lh.aparc);
    pfcAparcIndex = aparc([1,4,7,9,10,14,16,18,22,23,27,34,37,40,42,48]);

    
    % maybe 32
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing section to identify location of aparc regions
    locate = 0;
    if locate,
        index = find( subjectFS.lh.aparc == aparc(48));
        findparc = subjectFS.lh.curv;
        findparc(index) = findparc(index) * 20;
        freesurfer_plot_thickness(subjectFS.lh.inflated, findparc);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    lh.pfcVertexIndex = [];
    rh.pfcVertexIndex = [];
    for i = 1:length(pfcAparcIndex),
        aparcValue = pfcAparcIndex(i);
        % find lh
        index = find( subjectFS.lh.aparc == aparcValue );
        lh.pfcVertexIndex = union(lh.pfcVertexIndex,index);
        % find rh
        index = find( subjectFS.rh.aparc == aparcValue );
        rh.pfcVertexIndex = union(rh.pfcVertexIndex,index);
    end
    
    lh.pfcThickness = subjectFS.lh.thickness(lh.pfcVertexIndex);
    rh.pfcThickness = subjectFS.lh.thickness(lh.pfcVertexIndex);
    subplot(2,1,1); hist(lh.pfcThickness,bins); ylabel('LH')
    title(upper(strrep(subject,'nri_','Prefrontal Thickness - ')))
    subplot(2,1,2); hist(rh.pfcThickness,bins); ylabel('RH')
    xlabel('Cortical Thickness (mm)')
    imagefile = fullfile(subjectPath,[subject,'.histogram.pfcThickness.png']);
    save_png(imagefile,gcf,300)
    close(gcf)


    %     % plot thickness overlayed on curvature
    %     lh.findparc = subjectFS.lh.curv;
    %     lh.findparc(lh.pfcVertexIndex) = subjectFS.lh.thickness(lh.pfcVertexIndex);
    %     freesurfer_plot_thickness(subjectFS.lh.inflated, lh.findparc);

    lh.pfcGyrusIndex  = intersect(lh.pfcVertexIndex,lh.GyrusIndex);
    lh.pfcSulcusIndex = intersect(lh.pfcVertexIndex,lh.SulcusIndex);
    rh.pfcGyrusIndex  = intersect(rh.pfcVertexIndex,rh.GyrusIndex);
    rh.pfcSulcusIndex = intersect(rh.pfcVertexIndex,rh.SulcusIndex);

    % plot curv, use 1, otherwise 0
    plotcurv = 1;

    % pfc thickness for curv < 0
    lh.findparc = subjectFS.lh.curv * plotcurv;
    lh.findparc(lh.pfcGyrusIndex) = subjectFS.lh.thickness(lh.pfcGyrusIndex);

    rh.findparc = subjectFS.rh.curv * plotcurv;
    rh.findparc(rh.pfcGyrusIndex) = subjectFS.rh.thickness(rh.pfcGyrusIndex);

    freesurfer_plot_thickness(subjectFS.lh.inflated, lh.findparc);
    view([-90,0])
    imagefile = fullfile(subjectPath,[subject,'.lh.inflated.pfcThicknessByGyrus.png']);
    save_png(imagefile,gcf,300)
    close(gcf)
    freesurfer_plot_thickness(subjectFS.rh.inflated, rh.findparc);
    view([90,0])
    imagefile = fullfile(subjectPath,[subject,'.rh.inflated.pfcThicknessByGyrus.png']);
    save_png(imagefile,gcf,300)
    close(gcf)

    lh.pfcThickness = subjectFS.lh.thickness(lh.pfcGyrusIndex);
    rh.pfcThickness = subjectFS.rh.thickness(rh.pfcGyrusIndex);
    subplot(2,1,1); hist(lh.pfcThickness,bins); ylabel('LH')
    title(upper(strrep(subject,'nri_','Prefrontal Thickness at Curv < 0 - ')))
    subplot(2,1,2); hist(rh.pfcThickness,bins); ylabel('RH')
    xlabel('Cortical Thickness (mm)')
    imagefile = fullfile(subjectPath,[subject,'.histogram.pfcThicknessByGyrus.png']);
    save_png(imagefile,gcf,300)
    close(gcf)

    % pfc thickness for curv > 0
    lh.findparc = subjectFS.lh.curv * plotcurv;
    lh.findparc(lh.pfcSulcusIndex) = subjectFS.lh.thickness(lh.pfcSulcusIndex);

    rh.findparc = subjectFS.rh.curv * plotcurv;
    rh.findparc(rh.pfcSulcusIndex) = subjectFS.rh.thickness(rh.pfcSulcusIndex);

    freesurfer_plot_thickness(subjectFS.lh.inflated, lh.findparc);
    view([-90,0])
    imagefile = fullfile(subjectPath,[subject,'.lh.inflated.pfcThicknessBySulcus.png']);
    save_png(imagefile,gcf,300)
    close(gcf)
    freesurfer_plot_thickness(subjectFS.rh.inflated, rh.findparc);
    view([90,0])
    imagefile = fullfile(subjectPath,[subject,'.rh.inflated.pfcThicknessBySulcus.png']);
    save_png(imagefile,gcf,300)
    close(gcf)

    lh.pfcThickness = subjectFS.lh.thickness(lh.pfcSulcusIndex);
    rh.pfcThickness = subjectFS.rh.thickness(rh.pfcSulcusIndex);
    subplot(2,1,1); hist(lh.pfcThickness,bins); ylabel('LH')
    title(upper(strrep(subject,'nri_','Prefrontal Thickness at Curv > 0 - ')))
    subplot(2,1,2); hist(rh.pfcThickness,bins); ylabel('RH')
    xlabel('Cortical Thickness (mm)')
    imagefile = fullfile(subjectPath,[subject,'.histogram.pfcThicknessBySulcus.png']);
    save_png(imagefile,gcf,300)
    close(gcf)

    % compile summary statistics
    pfcThickness.mean(subjectN,1) = mean( subjectFS.lh.thickness(lh.pfcVertexIndex) );
    pfcThicknessByGyrus.mean(subjectN,1) = mean( subjectFS.lh.thickness(lh.pfcGyrusIndex) );
    pfcThicknessBySulcus.mean(subjectN,1) = mean( subjectFS.lh.thickness(lh.pfcSulcusIndex) );
    
    pfcThickness.mean(subjectN,2) = mean( subjectFS.rh.thickness(rh.pfcVertexIndex) );
    pfcThicknessByGyrus.mean(subjectN,2) = mean( subjectFS.rh.thickness(rh.pfcGyrusIndex) );
    pfcThicknessBySulcus.mean(subjectN,2) = mean( subjectFS.rh.thickness(rh.pfcSulcusIndex) );

    pfcThickness.std(subjectN,1) = std( subjectFS.lh.thickness(lh.pfcVertexIndex) );
    pfcThicknessByGyrus.std(subjectN,1) = std( subjectFS.lh.thickness(lh.pfcGyrusIndex) );
    pfcThicknessBySulcus.std(subjectN,1) = std( subjectFS.lh.thickness(lh.pfcSulcusIndex) );

    pfcThickness.std(subjectN,2) = std( subjectFS.rh.thickness(rh.pfcVertexIndex) );
    pfcThicknessByGyrus.std(subjectN,2) = std( subjectFS.rh.thickness(rh.pfcGyrusIndex) );
    pfcThicknessBySulcus.std(subjectN,2) = std( subjectFS.rh.thickness(rh.pfcSulcusIndex) );

end



[group,groupSortIndex]=sort(subjects.group);

fprintf('\npfcThickness\n');
fprintf('ID\tGP\tLH: M\tSD\tRH: M\tSD\n');
for i = 1:length(group),
    j = groupSortIndex(i);
    fprintf('%s\t%s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n',subjects.id{j},subjects.group{j},...
        pfcThickness.mean(j,1), pfcThickness.std(j,1), ...
        pfcThickness.mean(j,2), pfcThickness.std(j,2));
end

fprintf('\npfcThicknessByGyrus\n');
fprintf('ID\tGP\tLH: M\tSD\tRH: M\tSD\n');
for i = 1:length(group),
    j = groupSortIndex(i);
    fprintf('%s\t%s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n',subjects.id{j},subjects.group{j},...
        pfcThicknessByGyrus.mean(j,1), pfcThicknessByGyrus.std(j,1), ...
        pfcThicknessByGyrus.mean(j,2), pfcThicknessByGyrus.std(j,2));
end

fprintf('\npfcThicknessBySulcus\n');
fprintf('ID\tGP\tLH: Mean\tLH: SD\tRH: Mean\tRH: SD\n');
for i = 1:length(group),
    j = groupSortIndex(i);
    fprintf('%s\t%s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n',subjects.id{j},subjects.group{j},...
        pfcThicknessBySulcus.mean(j,1), pfcThicknessBySulcus.std(j,1), ...
        pfcThicknessBySulcus.mean(j,2), pfcThicknessBySulcus.std(j,2));
end


% lh_group = strrep(subjects.group,'C','LH_C');
% lh_group = strrep(lh_group,'SA','LH_S');
% rh_group = strrep(lh_group,'LH','RH');
% gp{1} = [lh_group, rh_group]';
% Y = [ pfcThickness.mean(:,1); pfcThickness.mean(:,2) ];
% anovan(Y,gp)
