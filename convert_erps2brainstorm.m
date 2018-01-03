
% convert_erps2brainstorm - script to convert ascii to brainstorm ERP data


%clear all


data = 'd:\matlab\brainstorm_v1\studies\';

cd(data);

datapath = pwd;

% Load the averaged data
sample_rate = 2.5;
epoch_start = -200;
epoch_end = 1500;
points = 681;

% --- Setup data structures for timing array

timeArray = meshgrid(epoch_start:sample_rate:epoch_end,1)';
timeNonZero = find(timeArray);
timeZero = find(timeArray == 0);

% Load data
for g = {'c','p'},
    
    if strcmp(g,'c'),
        
        for s = 1:10,
            
            sdir = sprintf('%s%02d',char(g),s);
            cd(sdir)
            
            for c = {'o','oac','oat','ouc','out','t','tac','tad','tat','tuc','tud','tut'},
                
                file = sprintf('%s%02d%s.dat',char(g),s,char(c));
                volt = sprintf('%s%02d%s',    char(g),s,char(c));
                
                channelfile = sprintf('%s%02d_channel.mat',char(g),s)
                if exist(channelfile) == 2,
                    load(channelfile,'Channel')
                    channelfile = sprintf('%s%02d_%s_volts_channel.mat',char(g),s,char(c))
                    save(channelfile,'Channel')
                    clear Channel
                end
                
                if ~exist(volt,'var')
                    
                    fprintf('loading ascii %s\n',file);
                    load(file);
                    
                    % Interpolate the Zero value
                    V = eval(volt)';
                    InterpZero = interp1( timeArray(timeNonZero), V, 0, 'cubic' );
                    V = [V(1:timeZero-1,:); InterpZero; V(timeZero:end,:)];
                    eval(strcat(volt,' = V'';'));
                    
                    % Convert and save into brainstorm file
                    
                    data.Time = timeArray' ./ 1000; % in seconds
                    data.F    = eval(volt) ./ 10^6; % in Volts
                    
                    %data.F(125,:) = zeros(1,size(data.F,2)); % create zero ref data
                    %data.F(125,:) = mean(data.F(1:124,:));   % create avg  ref data??
                    
                    fileprefix = sprintf('%s%02d_%s_volts',char(g),s,char(c));
                    eeg_write_brainstorm(fileprefix,data);
                    
                else
                    fprintf('Already loaded %s\n',file);
                end
                
                
            end
            cd ..
        end
        
    else
        
        for s = [2 4:9],
            
            sdir = sprintf('%s%02d',char(g),s);
            cd(sdir)
            
            for c = {'o','oac','oat','ouc','out','t','tac','tad','tat','tuc','tud','tut'},
                
                file = sprintf('%s%02d%s.dat',char(g),s,char(c));
                volt = sprintf('%s%02d%s',    char(g),s,char(c));
                
                channelfile = sprintf('%s%02d_channel.mat',char(g),s)
                if exist(channelfile) == 2,
                    load(channelfile,'Channel')
                    channelfile = sprintf('%s%02d_%s_volts_channel.mat',char(g),s,char(c))
                    save(channelfile,'Channel')
                    clear Channel
                end
                
                if ~exist(volt,'var')
                    fprintf('loading ascii %s\n',file);
                    load(file);
                    
                    % Interpolate the Zero value
                    V = eval(volt)';
                    InterpZero = interp1( timeArray(timeNonZero), V, 0, 'cubic' );
                    V = [V(1:timeZero-1,:); InterpZero; V(timeZero:end,:)];
                    eval(strcat(volt,' = V'';'));
                    
                    % Convert and save into brainstorm file
                    
                    data.Time = timeArray' ./ 1000; % in seconds
                    data.F    = eval(volt) ./ 10^6; % in Volts
                    
                    %data.F(125,:) = zeros(1,size(data.F,2)); % create zero ref data
                    %data.F(125,:) = mean(data.F(1:124,:));   % create avg  ref data??
                    
                    fileprefix = sprintf('%s%02d_%s_volts',char(g),s,char(c));
                    eeg_write_brainstorm(fileprefix,data);
                    
                else
                    fprintf('Already loaded %s\n',file);
                end
                
            end
            cd ..
        end
    end
end
clear timeZero timeNonZero g c file volt V InterpZero
