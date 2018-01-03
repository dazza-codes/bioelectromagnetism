function [] = importBDF(refChan, hiChan, timeChunk) 

% This bit of code will load all of the bdf files in a directory, break
% them into segments that can load without EEGLAB crashing, reference the 
% data to the specified channel, look up channel information, then merge and 
% save the data.
%
% This function requires three inputs: reference channel #, highest channel, 
% and length of time chunks in seconds.

% Bradley Voytek
% Copyright (c) 2007
% University of California, Berkeley
% Helen Wills Neuroscience Institute
% btvoytek@berkeley.edu

bdfDir = [pwd '/'];
mkdir([bdfDir 'SetFiles/']);

% Load file names into array
bdfFiles = dir(bdfDir);

% Determines if files in a directory are *.bdf files and, if so, adds a
% "process = true" tag to them
bdfStr = ...
    'bdfFiles(bdfIt).name((length(bdfFiles(bdfIt).name)-3):(length(bdfFiles(bdfIt).name)))';
for bdfIt = 1:length(bdfFiles)
    if  bdfFiles(bdfIt).isdir == 0
        if strcmp(eval(bdfStr), '.bdf') == 1
            bdfFiles(bdfIt).process = 1;
        elseif strcmp(eval(bdfStr), '.bdf') == 0
            bdfFiles(bdfIt).process = 0;
        end
    elseif bdfFiles(bdfIt).isdir == 1
        bdfFiles(bdfIt).process = 0;
    end
end
clear bdfIt bdfStr;
% **********

% Create filename structure containing only *.bdf files & get bytes
procCount = 1;

for bdfIt = 1:length(bdfFiles)
    if bdfFiles(bdfIt).process == 1
        bdf(procCount).bdfName = bdfFiles(bdfIt).name;
        bdf(procCount).Bytes = bdfFiles(bdfIt).bytes;
        procCount = procCount + 1;
    end
end
clear bdfStr bdfIt bdfFiles procCount;
% **********

% Add leading '0' to subject number for file-naming purposes
for zeroIt = 1:length(bdf)
    if zeroIt < 10
        bdf(zeroIt).subNum = ['0' num2str(zeroIt)];
    elseif zeroIt >= 10
        bdf(zeroIt).subNum = [num2str(zeroIt)];
    end
end

clear zeroIt;
% **********

% Loop load across subjects in directory
for subIt = 1:length(bdf)
    
    tempDir = [bdfDir bdf(subIt).subNum '/'];
    mkdir(tempDir);

    % Get SR, # of channels, and time from header
    dat = sopen([bdfDir bdf(subIt).bdfName]);
    bdf(subIt).SampRate = dat.SampleRate;
    bdf(subIt).Time = dat.NRec * dat.Dur;
    bdf(subIt).TotChannels = length(dat.Label);
    bdf(subIt).Samples = bdf(subIt).Time * bdf(subIt).SampRate;

    clear dat;
    % **********
        
    % Determines number of segments in which to break the data, in 
    bdf(subIt).Segments = ceil(bdf(subIt).Time / timeChunk);

    % Create start and end-times for data segmentation
    bdf(subIt).startArray = [];
    bdf(subIt).endArray = [];

    for segArrayIt = 1:bdf(subIt).Segments
        bdf(subIt).startArray = ...
            [bdf(subIt).startArray ((timeChunk * (segArrayIt - 1)) - (1 / bdf(subIt).SampRate))];
        bdf(subIt).endArray = [bdf(subIt).endArray (timeChunk * segArrayIt)];
    end

    bdf(subIt).startArray(1) = 0;
    bdf(subIt).endArray(segArrayIt) = bdf(subIt).Time;

    clear segArrayIt;
    % **********

    for segIt = 1:bdf(subIt).Segments
        EEG = zeros([bdf(subIt).TotChannels (bdf(subIt).SampRate * timeChunk)]);

        % Load raw data segment & reference to refChan
        EEG = ...
            pop_biosig([bdfDir bdf(subIt).bdfName], ...
            'blockrange', [bdf(subIt).startArray(segIt) bdf(subIt).endArray(segIt)], ...
            'ref', refChan, 'rmeventchan', 'off');
        
        % Remove extra channels
        EEG = pop_select(EEG, 'channel', [1:hiChan]);
        
        % Save dataset
        EEG = ...
            pop_saveset(EEG, 'filename', [bdf(subIt).subNum '-' num2str(segIt)], ...
            'filepath', tempDir);
    end
    % **********
    
    % Clear info for next subject
    clear EEG segIt;
    % **********
    
    % Load first dataset
    EEG = pop_loadset('filename', [bdf(subIt).subNum '-1.set'], 'filepath', tempDir);
    
    for segIt = 2:bdf(subIt).Segments
        % Load next datasets
        mergeEEG = ...
            pop_loadset('filename', [bdf(subIt).subNum '-' num2str(segIt) '.set'], ...
            'filepath', tempDir);
        
        % *******THIS WHOLE SECTION IS TAKEN MOSTLY FROM pop_mergeset.m******
        % Merge datasets
            % Concatenate data
            % ----------------
            if EEG.trials > 1 | mergeEEG.trials > 1
                EEG.data(:,:,end+1:end+size(mergeEEG.data,3)) = mergeEEG.data;
            else
                EEG.data(:,end+1:end+size(mergeEEG.data,2)) = mergeEEG.data;
            end;

            EEG.setname = 'Merged datasets';
            EEGtrials = EEG.trials;
            mergeEEGtrials = mergeEEG.trials;
            EEGpnts   = EEG.pnts;
            mergeEEGpnts   = mergeEEG.pnts;

            if EEG.trials > 1 | mergeEEG.trials > 1 % epoched data
                EEG.trials  =  EEG.trials + mergeEEG.trials;

            else % continuous data
                EEG.pnts = EEG.pnts + mergeEEG.pnts;
            end;

            if isfield(EEG, 'reject')
                EEG = rmfield(EEG, 'reject' );
            end;
            EEG.specicaact = [];
            EEG.specdata = [];
            EEG.icaact = [];
            EEG.icawinv = [];
            EEG.icasphere = [];
            EEG.icaweights = [];
            if isfield(EEG, 'stats')
                EEG = rmfield(EEG, 'stats' );
            end;
        
            % concatenate events
            % ------------------
            if isempty(mergeEEG.event)
                % boundary event
                % -------------
                disp('Inserting boundary event...');
                EEG.event(end+1).type    = 'boundary'; % add boundary event between datasets
                EEG.event(end).latency = EEGpnts+0.5; % make boundary halfway between last,first pts

                % check urevents
                % --------------
                if ~isfield(EEG, 'urevent'), 
                    EEG.urevent = []; 
                    fprintf('Warning: first dataset has no urevent structure.\n');
                end;

                % add boundary urevent
                % --------------------
                disp('Inserting boundary urevent...');
                EEG.urevent(end+1).type    = 'boundary';  

                if length(EEG.urevent) > 1 % if previous EEG urevents
                    EEG.urevent(end  ).latency = max(EEGpnts, EEG.urevent(end-1).latency)+0.5;
                else
                    EEG.urevent(end  ).latency = EEGpnts+0.5;
                end;

            else
                % concatenate urevents
                % --------------------
                if isfield(mergeEEG, 'urevent') 
                  if ~isempty(mergeEEG.urevent) 

                    % insert boundary event
                    % ---------------------
                    disp('Inserting boundary event...');
                    EEG.urevent(end+1).type    = 'boundary';
                    EEG.urevent(end  ).latency = max(EEGpnts, EEG.urevent(end-1).latency)+0.5;

                    % update urevent indices for second dataset
                    % -----------------------------------------
                    disp('Concatenating urevents...');
                    orilen    = length(EEG.urevent);
                    for e=1:length(mergeEEG.event)
                        mergeEEG.event(e).urevent = mergeEEG.event(e).urevent + orilen;
                    end;

                    allfields = fieldnames(mergeEEG.urevent);
                    for i=1:length( allfields )
                        for e=1:length(mergeEEG.urevent)
                            tmpval = getfield(mergeEEG.urevent, { e }, allfields{i});
                            EEG.urevent = setfield(EEG.urevent, {orilen + e}, allfields{i}, tmpval);
                        end
                    end
                  else
                    fprintf('Warning: second dataset has empty urevent structure.\n');
                  end 
                end

                % concatenate events
                % ------------------
                disp('Concatenating events...');
                orilen = length(EEG.event);
                %allfields = fieldnames(mergeEEG.event);

                % ensure similar event structures
                % -------------------------------
                if ~isempty(mergeEEG.event)
                    fields1 = lower(fieldnames(EEG.event));
                    fields2 = lower(fieldnames(mergeEEG.event));
                    if length(fields1) > length(fields2)
                        for index = 1:length(fields1)
                            if isempty(strmatch(fields1{index}, fields2))
                                mergeEEG.event = setfield( mergeEEG.event, { orilen + 1}, fields1{index}, []);
                            end
                        end
                    elseif length(fields1) < length(fields2)
                        for index = 1:length(fields2)
                            if isempty(strmatch(fields2{index}, fields1))
                                EEG.event = setfield( EEG.event, { 1 }, fields2{index}, []);
                            end
                        end
                    end
                end

                for e=1:length(mergeEEG.event)
                    EEG.event(orilen + e) = mergeEEG.event(e);
                    if isfield(EEG.event,'latency') & isfield(mergeEEG.event,'latency')
                       EEG.event(orilen + e).latency = mergeEEG.event(e).latency + EEGpnts * EEGtrials;
                    end
                    if isfield(EEG.event,'epoch') & isfield(mergeEEG.event,'epoch')
                       EEG.event(orilen + e).epoch = mergeEEG.event(e).epoch + EEGtrials;
                    end
                 end

                EEG.epoch = []; % epoch info regenerated below by 'eventconsistency' in eeg_checkset()

                % add discontinuity event if continuous
                % -------------------------------------
                if EEGtrials  == 1 & mergeEEGtrials == 1
                    disp('Adding boundary event...');
                    EEG.event = eeg_insertbound(EEG.event, EEG.pnts, EEGpnts+1, 0); % +1 since 0.5 is subtracted
                end
            end

            if isfield(EEG, 'epoch') && isfield(mergeEEG, 'epoch') && ~isempty(EEG.epoch) && ~isempty(mergeEEG.epoch)
                try 
                    EEG.epoch(end+1:end+mergeEEG.trials) = mergeEEG.epoch(:);
                catch
                    disp('pop_mergetset: epoch info removed (information not consistent across datasets)');
                end
            else
                EEG.epoch =[];
            end

            % check consistency of merged dataset, regenerate epoch structure
            % ---------------------------------------------------------------
            if ~isempty(mergeEEG.event)
                EEG.pnts = size(EEG.data,2);
                disp('Reconstituting epoch information...');
                EEG = eeg_checkset(EEG, 'eventconsistency');
            end
            
            clear EEGpnts EEGtrials allfields e fields1 fields2 i mergeEEGpnts mergeEEGtrials orilen tmpval;
        %*************************************
        
        clear mergeEEG;
    end

    % Edit set metadata
    EEG = pop_editset(EEG, 'setname', ['S' bdf(subIt).subNum], 'subject', bdf(subIt).subNum);
    % **********
    
    % Look up channel information
    EEG = pop_chanedit(EEG, 'lookup',...
        '/Applications/MATLAB74/toolbox/eeglab_current/eeglab2007May01_beta/plugins/dipfit2.1/standard_BESA/standard-10-5-cap385.elp');
    % **********
    
    % Save raw, re-referenced dataset
    EEG = pop_saveset(EEG, 'filename', bdf(subIt).subNum, 'filepath', [bdfDir 'SetFiles/']);
    % **********
    
    % Delete temp files
    delete([tempDir '*.set']);
    rmdir(tempDir);
    
    clear EEG segIt;
end