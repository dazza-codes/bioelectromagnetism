
% plots_loaddata - script to load and plot ERPs

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
        CONT.expfile = sprintf('%s%s%s%s.txt',char(g),exp,'_',data);
        CONT.expvolt = sprintf('%s%s%s%s',    char(g),exp,'_',data);
        CONT.confile = sprintf('%s%s%s%s.txt',char(g),con,'_',data);
        CONT.convolt = sprintf('%s%s%s%s',    char(g),con,'_',data);
        
        CONT.diffile = sprintf('%s%s%s%s%s%s.txt',char(g),exp,'-',con,'_',data);
        CONT.difvolt = sprintf('%s%s%s%s%s%s',    char(g),exp,'_',con,'_',data);
    else
        PTSD.expfile = sprintf('%s%s%s%s.txt',char(g),exp,'_',data);
        PTSD.expvolt = sprintf('%s%s%s%s',    char(g),exp,'_',data);
        PTSD.confile = sprintf('%s%s%s%s.txt',char(g),con,'_',data);
        PTSD.convolt = sprintf('%s%s%s%s',    char(g),con,'_',data);
        
        PTSD.diffile = sprintf('%s%s%s%s%s%s.txt',char(g),exp,'-',con,'_',data);
        PTSD.difvolt = sprintf('%s%s%s%s%s%s',    char(g),exp,'_',con,'_',data);
    end
    
    for c = {exp,con,dif},
        
        if findstr('-',char(c)),
            file = sprintf('%s%s%s%s%s%s.txt',char(g),exp,'-',con,'_',data);
            volt = sprintf('%s%s%s%s%s%s',    char(g),exp,'_',con,'_',data);
        else
            file = sprintf('%s%s%s%s.txt',char(g),char(c),'_',data);
            volt = sprintf('%s%s%s%s',    char(g),char(c),'_',data);
        end
        
        if ~exist(volt,'var')
            fprintf('loading %s\n',file);
            load(file);
            % Interpolate the Zero value
            V = eval(volt);
            InterpZero = interp1( timeArray(timeNonZero), V, 0, 'cubic' );
            V = [V(1:timeZero-1,:); InterpZero; V(timeZero:end,:)];
            eval(strcat(volt,' = V;'));
		else
            fprintf('Already loaded %s\n',file);
		end
    end
end

clear timeZero timeNonZero g c file volt V InterpZero
