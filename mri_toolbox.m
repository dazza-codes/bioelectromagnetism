function [ varargout ] = mri_toolbox(command)

% mri_toolbox - Graphical user interface (GUI) to various MRI tools
% 
% The main gui is the primary store for general parameters 
% and provides access to other tools.
%

% $Revision: 1.2 $ $Date: 2007/02/17 03:04:30 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2003, Darren.Weber_at_radiology.ucsf.edu
%                    adapted from eeg_toolbox
% 
% Depends:  various tools in the mri_toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('command','var'),
    command = 'init';
elseif isempty(command),
    command = 'init';
end

command = lower(command);

% first establish GUI or get userdata
switch command,
case 'init',
    MRITOOLBOX = init;    
otherwise,
    MRITOOLBOX = get(gcbf,'Userdata');
end


% now run commands, if required
switch command,
case 'init',
    % taken care of above
    
case 'openavw',
    gui_avw_open(MRITOOLBOX.mri,'init',MRITOOLBOX.gui);
    
case 'opencor',
    gui_cor_open(MRITOOLBOX.mri,'init',MRITOOLBOX.gui);
    
case 'openge',
    gui_ge_open(MRITOOLBOX.mri,'init',MRITOOLBOX.gui);
    
case 'defaultreturn',
    % keep this command so the return action is tidy
    %mri = MRITOOLBOX.mri; % done below just before return
    
case 'defaultreset',
    MRITOOLBOX.mri = mri_toolbox_defaults('create');
    
case 'defaultsave',
    mri_toolbox_defaults('write',MRITOOLBOX.mri);
    
case 'saveas',
    mri_toolbox_defaults('write_other',MRITOOLBOX.mri);
    
case 'recent',
    
    % -- get recent files list
    
    recentfiles = mri_toolbox_recent;
    
    % -- remove current menu items for recent files
    
    if isfield(MRITOOLBOX.menu,'recentfiles'),
        handleIndex = find(ishandle(MRITOOLBOX.menu.recentfiles));
        delete(MRITOOLBOX.menu.recentfiles(handleIndex));
    end
    MRITOOLBOX.menu.recentfiles = [];
    
    % -- recreate recent files menu items
    
    if and(size(recentfiles,2) == 1, isempty(recentfiles{1})),
        if ishandle(MRITOOLBOX.menu.recent),
            set(MRITOOLBOX.menu.recent,'Label','No Recent Files');
        end
    else
        if ishandle(MRITOOLBOX.menu.recent),
            set(MRITOOLBOX.menu.recent,'Label','Recent Files');
        end
        
        % -- add recent files to menu and setup their callbacks
        
        for i=1:length(recentfiles),
            if ~isempty(recentfiles{i}),
                MRITOOLBOX.menu.recentfiles(i) = uimenu(MRITOOLBOX.menu.recent,...
                    'Label',recentfiles{i},...
                    'Callback',strcat('[recentfiles,mri] = mri_toolbox_recent(''',...
                    recentfiles{i},''',''load''); ',...
                    'MRITOOLBOX = get(gcbf,''Userdata''); ',...
                    'gui_mri_open(mri,''init'',MRITOOLBOX.gui); ',...
                    'clear MRITOOLBOX recentfiles;'));
            end
        end
        
        % -- add recent files clear command
        
        MRITOOLBOX.menu.recentfiles(i+1) = uimenu(MRITOOLBOX.menu.recent,...
            'Label','Clear All',...
            'Callback',strcat('mri_toolbox_recent('''',''clear''); ',...
            'mri_toolbox(''recent''); '));
    end
    
case 'exit',
    close gcbf;
    
otherwise,
    fprintf('...invalid command to mri_toolbox\n\n');
    
end


switch command,
case 'exit',
otherwise,
    set(MRITOOLBOX.gui,'UserData',MRITOOLBOX);
end

if nargout > 0,
    if isfield(MRITOOLBOX,'mri'),
        if ~isempty(MRITOOLBOX.mri),
            varargout{1} = MRITOOLBOX.mri;
        end
    end
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paint the GUI
function [H] = init()

% Parameters are supplied in the file defaultfile.
H.mri = mri_toolbox_defaults('read');

GUIwidth  = 250;
GUIheight = 50;

H.gui = figure('Name','MRI Toolbox',...
               'Tag','MRI_TOOLBOX',...
               'NumberTitle','off',...
               'MenuBar','none');
set(H.gui,'Position',[1 1 GUIwidth GUIheight]);  % Activate GUI Figure
movegui(H.gui, 'center');

% -- file menu

H.menu.file_menu = uimenu(H.gui,'Label','File',...
    'Callback','mri_toolbox(''recent'');');

if exist('avw_read.m') == 2,
    H.menu.open_mri = uimenu(H.menu.file_menu,'Label','Open Analyze 7.5',...
        'Callback','mri = mri_toolbox(''openAVW'');','Accelerator','a');
end

if exist('cor_img_read.m') == 2,
    H.menu.open_mri = uimenu(H.menu.file_menu,'Label','Open FreeSurfer COR',...
        'Callback','mri = mri_toolbox(''openCOR'');','Accelerator','f');
end

if exist('ge_series_read.m') == 2,
    H.menu.open_mri = uimenu(H.menu.file_menu,'Label','Open GE Signa 5.x/LX',...
        'Callback','mri = mri_toolbox(''openGE'');','Accelerator','g');
end

H.menu.recent = uimenu(H.menu.file_menu,'Label','Recent');
H.menu.quit   = uimenu(H.menu.file_menu,'Label','Exit',...
    'Callback','mri_toolbox(''exit'');','Accelerator','x');

% -- Parameters menu

H.menu.p_menu = uimenu(H.gui,'Label','Parameters');
H.menu.show   = uimenu(H.menu.p_menu,'Label','Return to Workspace',...
    'Callback','mri = mri_toolbox(''defaultreturn'')');
H.menu.reset  = uimenu(H.menu.p_menu,'Label','Reset to Defaults',...
    'Callback','mri = mri_toolbox(''defaultreset'');');
H.menu.save   = uimenu(H.menu.p_menu,'Label','Save Defaults',...
    'Callback','mri = mri_toolbox(''defaultsave'');');
H.menu.saveas = uimenu(H.menu.p_menu,'Label','Save As Data Workspace',...
    'Callback','mri = mri_toolbox(''saveas'');');

% -- help menu

H.menu.Help = uimenu(H.gui,'Label','Help');
H.menu.help = uimenu(H.menu.Help,'Label','Help','Callback','doc mri_toolbox;');
H.menu.gpl  = uimenu(H.menu.Help,'Label','GPL','Callback','eeg_gpl;');

return
