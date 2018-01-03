

% plots_waveforms - script to plot/save ERP waveforms

clear

comp = 'wm';
data = 'link14hz';

% --- create individual electrode plots? NO = 0, YES = 1
doelec = 1;


% ---- Set the defaults for the titles and labels
if strmatch('scd14hz',data),
    titles.ylabel = '\muA/m^3';
else
    titles.ylabel = '\muV';
end
titles.title  = '';
titles.xlabel = 'msec';
titles.xalign = 'center';
titles.xxpos  = [];
titles.xypos  = [];
titles.yalign = 'right';
titles.yypos  = [];
titles.yxpos  = -300;


switch comp,
case 'sa',
    exp = 'oac'; con = 'ouc'; dif = 'oac-ouc';
case 'wm',
    exp = 'tac'; con = 'oac'; dif = 'tac-oac';
case 'ea',
    exp = 'oat'; con = 'oac'; dif = 'oat-oac';
case 'dt',
    exp = 'tud'; con = 'tuc'; dif = 'tud-tuc';
end

plots_loaddata

if mkdir('elec_plots'),
    cd('elec_plots');
    if mkdir(comp),
        cd(comp);
    else
        msg = sprintf('Cannot create %s directory',comp);
        error(msg);
    end
else
    error('Cannot create elec_plots directory');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the exp/con conditions at all electrodes

expvolt = eval(CONT.expvolt);
convolt = eval(CONT.convolt);
CONTexpvolt = expvolt;
CONTconvolt = convolt;

expvolt = eval(PTSD.expvolt);
convolt = eval(PTSD.convolt);
PTSDexpvolt = expvolt;
PTSDconvolt = convolt;


% -- Determine common scale

xmin = min(timeArray);
xmax = max(timeArray);

CONT.erpmax = max([CONTexpvolt CONTconvolt]);
CONT.erpmin = min([CONTexpvolt CONTconvolt]);
CONT.ymax = ceil( max(CONT.erpmax));
CONT.ymin = floor(min(CONT.erpmin));

PTSD.erpmax = max([PTSDexpvolt PTSDconvolt]);
PTSD.erpmin = min([PTSDexpvolt PTSDconvolt]);
PTSD.ymax = ceil( max(PTSD.erpmax));
PTSD.ymin = floor(min(PTSD.erpmin));

ymin = min(CONT.ymin,PTSD.ymin);
ymax = max(CONT.ymax,PTSD.ymax);

scale = [xmin,xmax,ymin,ymax];


% -- Plot the data

interval = [100 1]; % tick marks 100 ms, 1 uV

F = eeg_plot(1,timeArray,CONTexpvolt,scale,interval,titles);
file = strcat(exp,'_page_cont.png');
saveas(F,file)
close(F)

F = eeg_plot(1,timeArray,PTSDexpvolt,scale,interval,titles);
file = strcat(exp,'_page_ptsd.png');
saveas(F,file)
close(F)

F = eeg_plot(1,timeArray,CONTconvolt,scale,interval,titles);
file = strcat(con,'_page_cont.png');
saveas(F,file)
close(F)

F = eeg_plot(1,timeArray,PTSDconvolt,scale,interval,titles);
file = strcat(con,'_page_ptsd.png');
saveas(F,file)
close(F)

interval = [200 1]; % tick marks 200 ms, 1 uV

F = eeg_plot(0,timeArray,CONTexpvolt,scale,interval,titles);
file = strcat(exp,'_column_cont.png');
saveas(F,file)
close(F)

F = eeg_plot(0,timeArray,PTSDexpvolt,scale,interval,titles);
file = strcat(exp,'_column_ptsd.png');
saveas(F,file)
close(F)

F = eeg_plot(0,timeArray,CONTconvolt,scale,interval,titles);
file = strcat(con,'_column_cont.png');
saveas(F,file)
close(F)

F = eeg_plot(0,timeArray,PTSDconvolt,scale,interval,titles);
file = strcat(con,'_column_ptsd.png');
saveas(F,file)
close(F)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the difference data at all electrodes

CONTdifvolt = eval(CONT.difvolt);

PTSDdifvolt = eval(PTSD.difvolt);

% -- Determine common scale

xmin = min(timeArray);
xmax = max(timeArray);

CONT.difmax = max(CONTdifvolt);
CONT.difmin = min(CONTdifvolt);
CONT.ymax = ceil( max(CONT.difmax));
CONT.ymin = floor(min(CONT.difmin));

PTSD.difmax = max(PTSDdifvolt);
PTSD.difmin = min(PTSDdifvolt);
PTSD.ymax = ceil( max(PTSD.difmax));
PTSD.ymin = floor(min(PTSD.difmin));

ymax = ceil( max([CONT.ymax PTSD.ymax]));
ymin = floor(min([CONT.ymin PTSD.ymin]));

scale = [xmin,xmax,ymin,ymax];


% -- Plot the data

interval = [100 1]; % tick marks 100 ms, 1 uV

F = eeg_plot(1,timeArray,CONTdifvolt,scale,interval,titles);
file = strcat(comp,'_dif_page_cont.png');
saveas(F,file)
close(F)

F = eeg_plot(1,timeArray,PTSDdifvolt,scale,interval,titles);
file = strcat(comp,'_dif_page_ptsd.png');
saveas(F,file)
close(F)

interval = [200 1]; % tick marks 200 ms, 1 uV

F = eeg_plot(0,timeArray,CONTdifvolt,scale,interval,titles);
file = strcat(comp,'_dif_column_cont.png');
saveas(F,file)
close(F)

F = eeg_plot(0,timeArray,PTSDdifvolt,scale,interval,titles);
file = strcat(comp,'_dif_column_ptsd.png');
saveas(F,file)
close(F)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return now unless individual electrode plots required

if ~doelec, return; end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the exp/con conditions at each electrode

for elec = 1:124,
    
    expvolt = eval(CONT.expvolt);
    convolt = eval(CONT.convolt);
    CONT.erp(:,1) = expvolt(:,elec);
    CONT.erp(:,2) = convolt(:,elec);
    
    expvolt = eval(PTSD.expvolt);
    convolt = eval(PTSD.convolt);
    PTSD.erp(:,1) = expvolt(:,elec);
    PTSD.erp(:,2) = convolt(:,elec);
    
    
    % -- Determine common scale
    
    xmin = min(timeArray);
    xmax = max(timeArray);
    
    CONT.erpmax = max(CONT.erp);
    CONT.erpmin = min(CONT.erp);
    CONT.ymax = ceil( max(CONT.erpmax));
    CONT.ymin = floor(min(CONT.erpmin));
    
    PTSD.erpmax = max(PTSD.erp);
    PTSD.erpmin = min(PTSD.erp);
    PTSD.ymax = ceil( max(PTSD.erpmax));
    PTSD.ymin = floor(min(PTSD.erpmin));
    
    ymin = min(CONT.ymin,PTSD.ymin);
    ymax = max(CONT.ymax,PTSD.ymax);
    
    scale = [xmin,xmax,ymin,ymax];
    
    
    % -- Plot the data
    
    if mkdir('page'),
        cd('page');
    else
        error('Cannot create page directory');
    end
    
    interval = [100 1]; % tick marks 100 ms, 1 uV
    
    F = eeg_plot(1,timeArray,CONT.erp,scale,interval,titles);
    file = strcat(comp,'_page_elec',sprintf('%03d',elec),'_cont.png');
    saveas(F,file)
    close(F)
    
    F = eeg_plot(1,timeArray,PTSD.erp,scale,interval,titles);
    file = strcat(comp,'_page_elec',sprintf('%03d',elec),'_ptsd.png');
    saveas(F,file)
    close(F)
    
    cd ..
    if mkdir('column'),
        cd('column');
    else
        error('Cannot create column directory');
    end
    
    interval = [200 1]; % tick marks 200 ms, 1 uV
    
    F = eeg_plot(0,timeArray,CONT.erp,scale,interval,titles);
    file = strcat(comp,'_column_elec',sprintf('%03d',elec),'_cont.png');
    saveas(F,file)
    close(F)
    
    F = eeg_plot(0,timeArray,PTSD.erp,scale,interval,titles);
    file = strcat(comp,'_column_elec',sprintf('%03d',elec),'_ptsd.png');
    saveas(F,file)
    close(F)
    
    cd ..
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the difference data at each electrode

for elec = 1:124,
    
    difvolt = eval(CONT.difvolt);
    erp(:,1) = difvolt(:,elec);
    
    difvolt = eval(PTSD.difvolt);
    erp(:,2) = difvolt(:,elec);
    
    % -- Determine common scale
    
    xmin = min(timeArray);
    xmax = max(timeArray);
    
    erpmax = max(erp);
    erpmin = min(erp);
    ymax = ceil( max(erpmax));
    ymin = floor(min(erpmin));
    
    scale = [xmin,xmax,ymin,ymax];
    
    
    % -- Plot the data
    
    if mkdir('page'),
        cd('page');
    else
        error('Cannot create page directory');
    end
    
    interval = [100 1]; % tick marks 100 ms, 1 uV
    
    F = eeg_plot(1,timeArray,erp,scale,interval,titles);
    file = strcat(comp,'_dif_page_elec',sprintf('%03d',elec),'.png');
    saveas(F,file)
    close(F)
    
    cd ..
    if mkdir('column'),
        cd('column');
    else
        error('Cannot create column directory');
    end
    
    interval = [200 1]; % tick marks 200 ms, 1 uV
    
    F = eeg_plot(0,timeArray,erp,scale,interval,titles);
    file = strcat(comp,'_dif_column_elec',sprintf('%03d',elec),'.png');
    saveas(F,file)
    close(F)
    
    cd ..
    
end



return
