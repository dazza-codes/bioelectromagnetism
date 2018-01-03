function progress_bar(action,progress,title)

% progress_bar - Display a 'Progress Bar'
% 
% Usage: progress_bar('init',progress,title)
%         Initialises the progress bar.  Input 'progress' is
%         optional, assumed zero, otherwise a decimal
%         value between 0:1.  Input 'title' is optional,
%         but can be used to indicate the process of the
%         progress bar.
%
% Usage: progress_bar('Set',progress)
%         Updates the progress bar.  Input 'progress' is
%         a decimal value between 0:1.
%
% Usage: progress_bar('Clear')
%         Clears the progress bar.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% @(#)spm_progress_bar.m	2.1 John Ashburner 99/05/17
% Modified 03/2002, Darren.Weber_at_radiology.ucsf.edu
%                   - removed spm specific references
%                   - modified inputs/input handling to my liking
%                   - progress bar is now a patch object and the
%                     erasemode is xor with no figure backingstore
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isequal(nargin,0),
	action = 'init';
end

action = lower(action);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
if strcmp(action,'init'),
    
    tim = clock;
    timestr = sprintf('  Began %02.0f:%02.0f:%02.0f',tim(4),tim(5),tim(6));
	
    if exist('title','var'),
        if ~isempty(title),
            name = strcat(sprintf('%s  -',title),timestr);
        else
            name = strcat('Progress  -',timestr);
        end
    else
        name = strcat('Progress  -',timestr);
    end
    
    fg = figure('MenuBar','none',...
                'NumberTitle','off',...
                'Name',name,...
                'Tag','ProgressBar',...
                'BackingStore','off',...
                'Pointer','watch',...
                'Position',[1 1 300 75]);
    movegui(fg,'center');
    ax = axes('Position',[0.1 0.25 0.8 0.4],...
              'YTick',[],...
              'Xlim', [0 1],'Ylim',[0 1],...
              'Box','on',...
              'FontSize',8,...
              'Tag','ProgressBarAxis',...
              'Parent',fg);
    
    xlab = get(ax,'xticklabel');
    xlab = str2num(xlab) * 100;
    xlab = num2str(xlab);
    set(ax,'xticklabel',xlab)
    
    if exist('progress','var'),
        if ~isempty(progress), setpb(fg,progress);
        else
            progress = 0;
            setpb(fg,progress);
        end
    else
        progress = 0;
        setpb(fg,progress);
    end
    
    drawnow;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set
elseif strcmp(action,'set'),
    
    fg = findobj('Tag','ProgressBar');
    
    if ~isempty(fg), setpb(fg,progress); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
elseif strcmp(action,'clear'),
    fg = findobj('Tag','ProgressBar');
    if ~isempty(fg), close(fg);	end;
end;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setpb(fig,progress)
    
    pbaxis = findobj(fig,'Tag','ProgressBarAxis');
    
    if ~isempty(pbaxis),
        
        vert = [0 0; progress 0; progress 1; 0 1];
        face = [1 2 3 4];
        
        pbpatch = findobj(fig,'Tag','ProgressBarPatch');
        
        if ~isempty(pbpatch),
            set(pbpatch,'Vertices',vert);
        else
            pbpatch = patch('Faces',face,'Vertices',vert,'FaceColor','r',...
                'Tag','ProgressBarPatch',...
                'EraseMode','none',...
                'Parent',pbaxis);
        end
        
        title = get(pbaxis,'Title');
        set(title,'string',sprintf('%5.1f%% Complete',100*progress),'EraseMode','xor');
        
        drawnow;
        
        figure(fig);
    end;
return
