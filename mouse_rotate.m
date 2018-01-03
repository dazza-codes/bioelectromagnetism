function mouse_rotate(command,H)

% mouse_rotate - Left click to 3Drotate a matlab figure graphic
%
% mouse_rotate(command,H)
% 
% The mouse 'up' command switches between rotating and not
% rotating the figure.  The mouse 'rotate' command, when enabled,
% rotates the figure in the direction of mouse movement 
% within the figure.
% 
% H is the figure userdata, unless provided otherwise.  It
% contains various fields, see the .m file for details.
% 
% The figure context menu holds access to this function
% even after 'exit' (ie, right click in figure window
% away from any axis, data, colorbar, etc objects).
% 

% $Revision: 1.2 $ $Date: 2007/04/30 18:00:49 $

%Licence:            GNU GPL, no implied or express warranties
%Created:  10/1995 - Eric Soroos
%Modified: 02/2002 - Darren.Weber_at_radiology.ucsf.edu
%                  - integrated all commands into one function .m file
%                  - changed rotation to mouse direction rather than
%                    anti-mouse direction and during mouse up rather
%                    than mouse down.
%                  - changed target to the 'view' property of gca
%                  - added elevation,azimuth,view2/3,exit uicontrols
%                    and uicontext menu.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('command','var'), command = 'init'; end


command = lower(command);

switch command,
  
  case 'init',
    
    if exist('H','var'),
      if ~isempty(H),
        H = INIT(H);
      end
    else
      H = get(gcbf,'userdata');
      H = INIT(H);
    end
    
  case 'up',
    
    if ~exist('H','var'), H = get(gcbf,'userdata'); end
    if isempty(H),        H = get(gcbf,'userdata'); end
    if isempty(H),        H = get(gcf, 'userdata'); end
    
    %This switches between rotation and no rotation.
    %On a button press, it checks whether the motion function is
    %set - if yes, it clears it, otherwise it adds it.
    if isempty(get(H.gui,'WindowButtonMotionFcn')),
      set(H.gui,'WindowButtonMotionFcn','mouse_rotate(''rotate'');');
      %set(H.gui,'WindowButtonDownFcn','mouse_rotate(''up''); ');
      %set(H.info,'Value',1);
    else
      H.pos.current = get(H.axis,'View');
      %fprintf('\nCurrent Figure View: Az = %6.4f\tEl = %6.4f\n\n',H.pos.current);
      set(H.gui,'WindowButtonMotionFcn',[]);
      %set(H.gui,'WindowButtonDownFcn',[]);
      %set(H.info,'Value',0);
    end
    
  case 'rotate',
    
    if ~exist('H','var'), H = get(gcbf,'userdata'); end
    if isempty(H),        H = get(gcbf,'userdata'); end
    if isempty(H),        H = get(gcf, 'userdata'); end
    
    % rotates the figure in the direction of mouse movement
    % i.e. in the center, 0 altitude, 0 azimuth. Going up and down
    % changes the altitude, left and right changes the azimuth.
    H.pos.fig = get(H.gui,'position');
    H.pos.mouse = get(0,'pointerlocation');
    H.pos.relative = -1 * (((H.pos.mouse(1:2) - H.pos.fig(1:2))./ H.pos.fig(3:4))-.5).*[360,190];
    set(H.axis,'View',H.pos.relative);
    if ishandle(H.Az), set(H.Az,'String',sprintf('Az = %4.0f',H.pos.relative(1))); end
    if ishandle(H.El), set(H.El,'String',sprintf('El = %4.0f',H.pos.relative(2))); end
    
    % Adjust the light angle
    if isfield(H,'light'),
      if ishandle(H.light), delete(H.light); end
      H.light = camlight('headlight','infinite');
      set(H.light,'color',[1 1 1]/length(H.light)/1.2);
    end
    %     else
    %         H.light = camlight('headlight','infinite');
    %         set(H.light,'color',[1 1 1]/length(H.light)/1.2);
    %     end
    
    set(H.gui,'userdata',H);
    
  case 'view',
    
    if ~exist('H','var'), H = get(gcbf,'userdata'); end
    if isempty(H),        H = get(gcbf,'userdata'); end
    if isempty(H),        H = get(gcf, 'userdata'); end
    
    v = get(H.View,'Value');
    set(H.axis,'View',H.data.view(v,:));
    set(H.Az,'String',sprintf('Az = %4.0f',H.data.view(v,1)));
    set(H.El,'String',sprintf('El = %4.0f',H.data.view(v,2)));
    
    % Find and remove all light children of axes
    AC = get(H.axis,'Children');
    for c = 1:length(AC),
      ctype = get(AC(c),'type');
      if isequal(ctype,'light'),
        delete(AC(c));
      end
    end
    
    % Create a new light
    H.light = camlight('headlight','infinite');
    set(H.light,'color',[1 1 1]/length(H.light)/1.2);
    
    set(H.gui,'userdata',H);
    
  otherwise,
    
end

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = INIT(H),

%H = get(H,'userdata');

% only one per figure
if isfield(H,'El'),
  if ~isempty(H.El),
    return;
  end
end

H.gui = gcf;
H.axis = gca;

set(H.gui,'DoubleBuffer','on','WindowButtonDownFcn','mouse_rotate(''up''); ');

% enable right click access to mouse_rotate
menu=uicontextmenu;
a=uimenu(menu,'Label','Rotate','Callback','mouse_rotate; ');
if isempty(get(H.gui,'uicontextmenu')),
  set(H.gui,'uicontextmenu',menu);
end
if isempty(get(H.axis,'uicontextmenu')),
  set(H.axis,'uicontextmenu',menu);
  axsibs = get(H.axis,'Children');
  for i=1:length(axsibs),
    type = get(axsibs(i),'Type');
    if isequal(type,'patch'),
      if isempty(get(axsibs(i),'uicontextmenu')),
        set(axsibs(i),'uicontextmenu',menu);
      end
    end
  end
end

% Match background figure colour
bgcolor = get(H.gui,'Color');
% Try to adapt the foreground colour a little
black = find(bgcolor <= .6);
fgcolor = [0 0 0]; %black text
if length(black)>2, fgcolor = [1 1 1]; end

Font.FontName   = 'Helvetica';
Font.FontUnits  = 'Pixels';
Font.FontSize   = 12;
Font.FontWeight = 'bold';
Font.FontAngle  = 'normal';

H.Az  = uicontrol('Parent',H.gui,'Style','text',Font,...
  'Units','Normalized','Position',[.05 .0 .1 .05],...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','','TooltipString','Azimuth','HorizontalAlignment','left');
H.El  = uicontrol('Parent',H.gui,'Style','text',Font,...
  'Units','Normalized','Position',[.15 .0 .1 .05],...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','','TooltipString','Elevation','HorizontalAlignment','left');

H.Info = uicontrol('Parent',H.gui,'Style','radiobutton',Font,...
  'Units','Normalized','Position',[.65 .0 .1 .05],...
  'HorizontalAlignment','left',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','Rotate','Value',1,...
  'TooltipString','When on, click in figure to switch on/off 3D rotation.',...
  'Callback',strcat('H = get(gcbf,''userdata''); ',...
  'rotate = get(H.Info,''Value''); ',...
  'if rotate, ',...
  'set(H.gui,''WindowButtonDownFcn'',''mouse_rotate(''''up'''');''); ',...
  'else, ',...
  'set(H.gui,''WindowButtonDownFcn'',[]); ',...
  'set(H.gui,''WindowButtonMotionFcn'',[]); ',...
  'end; ',...
  'clear H rotate;'));

H.data.viewstr = {'2D' '3D' 'Top' 'Bottom' 'Front' 'Back' 'Left' 'Right'};
H.data.view = [ 0 90; -37.5 30; 0 90; 0 -90; 180 0; 0 0; -90 0; 90 0 ];

H.View = uicontrol('Parent',H.gui,'Style','Popup',...
  'Units','Normalized','Position',[.78 .0 .1 .05],...
  'Tag','VIEW','BackGroundColor','w',...
  'TooltipString','Select view.',...
  'String',H.data.viewstr,...
  'CallBack','mouse_rotate(''view'');');

H.RClose = uicontrol('Parent',H.gui,'Style','pushbutton',...
  'Units','Normalized','Position',[.9 .0 .1 .05],...
  'String','Close','Value',0,...
  'TooltipString','Use right click context menu after ''Close'' to get 3D rotation back.',...
  'Callback',strcat('H = get(gcbf,''userdata''); ',...
  'delete(H.RClose); delete(H.Info); ',...
  'delete(H.El); delete(H.Az); ',...
  'delete(H.View); ',...
  'fields = fieldnames(H); ',...
  'trashfields(1) = strmatch(''El'',   fields,''exact''); ',...
  'trashfields(2) = strmatch(''Az'',   fields,''exact''); ',...
  'trashfields(3) = strmatch(''Info'', fields,''exact''); ',...
  'trashfields(4) = strmatch(''View'', fields,''exact''); ',...
  'trashfields(5) = strmatch(''RClose'',fields,''exact''); ',...
  'trashfields(6) = strmatch(''pos'',fields,''exact''); ',...
  'H = setfield(H,char(fields(trashfields(1))),[]); ',...
  'H = setfield(H,char(fields(trashfields(2))),[]); ',...
  'H = setfield(H,char(fields(trashfields(3))),[]); ',...
  'H = setfield(H,char(fields(trashfields(4))),[]); ',...
  'H = setfield(H,char(fields(trashfields(5))),[]); ',...
  'H = setfield(H,char(fields(trashfields(6))),[]); ',...
  'H.data = setfield(H.data,''view'',[]); ',...
  'H.data = setfield(H.data,''viewstr'',[]); ',...
  'set(H.gui,''WindowButtonMotionFcn'',[]); ',...
  'set(H.gui,''WindowButtonDownFcn'',[]); ',...
  'set(H.gui,''WindowButtonUpFcn'',[]); ',...
  'set(H.gui,''userdata'',H); clear H trashfields fields; '));

H.pos.current = get(H.axis,'View');
set(H.Az,'String',sprintf('Az = %4.0f',H.pos.current(1)));
set(H.El,'String',sprintf('El = %4.0f',H.pos.current(2)));

set(H.gui,'userdata',H);
return
