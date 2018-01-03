function avw_view3(avw,FV),

% AVW_VIEW3 - Create slices of Analyze file and overlay surface vertices
% 
% avw_view3(avw,FV)
% 
% avw - a struct, created by avw_img_read
% FV  - a surface struct with fields FV.vertices, FV.faces
% 
% This provides a transparent 3D view of ortho slices, to
% provide a platform for overlay of mesh/tesselation vertices,
% given in FV.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/2002, Darren.Weber@flinders.edu.au
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('avw','var'),
    msg = sprintf('AVW_VIEW3: No input avw - see help avw_view\n');
    error(msg);
end
if ~exist('FV','var'),
    FV.vertices = [];
end


% GUI General Parameters
GUIwidth  = 150;
GUIheight = 50;

name = 'AVW View 3D';
if isfield(avw,'fileprefix'),
    if ~isempty(avw.fileprefix),
        format = strcat('%+',sprintf('%d',length(avw.fileprefix)+1),'s');
        name = strcat('AVW View 3D - ',sprintf(format,avw.fileprefix));
    end
end

GUI = figure('Name',name,'Tag','AVWVIEW3','units','characters',...
             'NumberTitle','off',...
             'MenuBar','figure','Position',[1 1 GUIwidth GUIheight]);
movegui(GUI,'center');

AVWVIEW3.gui = GUI;


Font.FontName   = 'Helvetica';
Font.FontUnits  = 'Pixels';
Font.FontSize   = 12;
Font.FontWeight = 'normal';
Font.FontAngle  = 'normal';


shading flat

xdim = size(avw.img,1);
ydim = size(avw.img,2);
zdim = size(avw.img,3);

SagSlice = 1;
CorSlice = 1;
AxiSlice = 1;
if xdim > 1, SagSlice = floor(xdim/2); end
if ydim > 1, CorSlice = floor(ydim/2); end
if zdim > 1, AxiSlice = floor(zdim/2); end


% data aspect should be proportional to image dimensions
% Should check header for these values


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The slice command swaps x/y axes - really annoying!
G.H = slice(avw.img,SagSlice,CorSlice,AxiSlice,'nearest');
set(G.H,'edgecolor','none','facealpha',.7);
axis tight, daspect([1,1,1]);
colormap('gray');
xlabel('Coronal')
ylabel('Sagittal')
zlabel('Axial')
view(3)

%axis off

% Overlay plot of surface vertices
if ~isempty(FV.vertices),
    hold on;
    % The above slice command swaps x/y axes, so we have to
    % mess around here to get the vertices in the right place
    
    % rotate the vertices -90 degrees around Z
    FV.vertices = Rz(FV.vertices,-90,'degrees');
    % Now move vertices from axis origin to the true MRI volume center
    FV.vertices(:,1) = FV.vertices(:,1) + SagSlice;
    FV.vertices(:,2) = FV.vertices(:,2) + CorSlice;
    FV.vertices(:,3) = FV.vertices(:,3) + AxiSlice;
    % Now shift it according to the x/y flip
    FV.vertices(:,1) = FV.vertices(:,1) - (SagSlice - CorSlice);
    FV.vertices(:,2) = FV.vertices(:,2) - (CorSlice - SagSlice);
    
    %plot3(FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3),'r.');
    patch('vertices',FV.vertices,'faces',FV.faces,...
          'FaceColor',[0.0 0.0 0.2],'Edgecolor','none',...
          'FaceAlpha',0.5);
end

rotate3d;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intensity Value at Mouse Click

%G.Timval = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
%    'Position',[.05 .90 .20 .05], 'HorizontalAlignment', 'center',...
%    'BusyAction','queue',...
%    'String','Image Intensity');
%G.imval = uicontrol('Parent',GUI,'Style','text','Units','Normalized', Font, ...
%    'Position',[.25 .90 .10 .05], 'HorizontalAlignment', 'right',...
%    'BusyAction','queue',...
%    'String','?');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Font.FontWeight = 'bold';

% OK: Return the avw!
G.Bhdr = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
    'Position',[.8 .95 .08 .04],...
    'String','HDR','BusyAction','queue',...
    'TooltipString','View the hdr parameters.',...
    'BackgroundColor',[0.0 0.0 0.5],...
    'ForegroundColor',[1 1 1], 'HorizontalAlignment', 'center',...
    'Callback',strcat('AVWVIEW3 = get(gcbf,''Userdata''); ',...
    'avw_view_hdr(AVWVIEW3.avw);',...
    'clear AVWVIEW3;'));

% Cancel
G.Bquit = uicontrol('Parent',GUI,'Style','pushbutton','Units','Normalized', Font, ...
    'Position',[.9 .95 .08 .04],...
    'String','CLOSE','BusyAction','queue',...
    'BackgroundColor',[0.75 0.0 0.0],...
    'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
    'Callback','close gcbf;');


% Update the gui_struct handles for this gui
AVWVIEW3.avw = avw;
AVWVIEW3.FV = FV;
AVWVIEW3.handles = G;
set(AVWVIEW3.gui,'Userdata',AVWVIEW3);
%set(AVWVIEW3.gui,'HandleVisibility','callback');

return

