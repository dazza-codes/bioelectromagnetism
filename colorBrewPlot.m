function colorBrewPlot(command) %cname,ncolor,ctype)

% colorBrewPlot
%
% This product includes color specifications and designs developed by
% Cynthia Brewer (http://colorbrewer.org/).
%
% Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania
% State University. All rights reserved.  Redistribution and use in source
% and binary forms, with or without modification, are permitted provided
% that the following conditions are met: 1. Redistributions as source code
% must retain the above copyright notice, this list of conditions and the
% following disclaimer.  2. The end-user documentation included with the
% redistribution, if any, must include the following acknowledgment: This
% product includes color specifications and designs developed by Cynthia
% Brewer (http://colorbrewer.org/).  Alternately, this acknowledgment may
% appear in the software itself, if and wherever such third-party
% acknowledgments normally appear.  4. The name "ColorBrewer" must not be
% used to endorse or promote products derived from this software without
% prior written permission. For written permission, please contact Cynthia
% Brewer at cbrewer@psu.edu.  5. Products derived from this software may not
% be called "ColorBrewer", nor may "ColorBrewer" appear in their name,
% without prior written permission of Cynthia Brewer.  THIS SOFTWARE IS
% PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT
% NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CYNTHIA BREWER,
% MARK HARROWER, OR THE PENNSYLVANIA STATE UNIVERSITY BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
% OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
% SUCH DAMAGE.
%


if ~exist('command','var'),
  command = 'init';
elseif isempty(command),
  command = 'init';
end

command = lower(command);

switch command,
 case 'init',
  CB = init;
  command = 'plot';
 otherwise,
  CB = get(gcbf,'Userdata');
end


switch command,
  
 case 'init',
  
 case 'plot',
  
  typeIndex = strmatch(CB.ctype,{CB.cs(:).colorType},'exact');
  ct = CB.cs(typeIndex);
  
  nameIndex = strmatch(CB.cname,{ct(:).colorName},'exact');
  cn = ct(nameIndex);
  
  titleStr = ['Color Scales (type: ', CB.ctype, ', name: ', CB.cname, ')'];
  set(CB.gui,'Name',titleStr)
  
  % create one figure for all the numColor definitions
  
  figure(CB.gui)
  set(CB.gui,'NextPlot','replace')
  
  
  numColors = unique([cn(:).numColors]);
  n = length(numColors);
  for c = 1:n,
    cm = cn(c).colorRGB;
    colormap(cm)
    h = subplot(1,n,c);
    cla(h)
    set(h,'NextPlot','replace')
    axis off
    cb = colorbar('peer',h);
    yticks = 1:size(cm,1)+1;
    set(cb,'YTick',yticks);
    %set(cb,'YTickLabel',num2str(yticks)');
  end
  
  
 case 'exit',
  close gcbf;
  
  
 otherwise,
  fprintf('...invalid command to colorBrewPlot\n\n');
  
end


switch command,
 case 'exit',
 otherwise,
  set(CB.gui,'UserData',CB);
end

if nargout > 0,
  if isfield(CB,'map'),
    if ~isempty(CB.map),
      varargout{1} = CB.map;
    end
  end
end

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paint the GUI
function [H] = init()

% Parameters are supplied in the file defaultfile.
H.cs = colorBrewScales;

H.ctype = H.cs(1).colorType;            % color type
H.cname = H.cs(1).colorName;            % color name
H.ncolor = [];                          % number of colors


GUIwidth  = 500;
GUIheight = 500;

name = 'Color Scales';

H.gui = figure('Name',name,...
               'Tag','COLORBREW',...
               'NumberTitle','off',...
               'NextPlot','replace',...
               'MenuBar','none');

set(H.gui,'Position',[1 1 GUIwidth GUIheight]);  % Activate GUI Figure
movegui center;


license = ['This product includes color specifications\n',...
           'and designs developed by Cynthia Brewer\n',...
           '(http://colorbrewer.org/).\n\n'];
fprintf(license);


% -- color scale select menu

H.menu.select = uimenu(H.gui,'Label','Select');

colorTypes = unique({H.cs(:).colorType});
for i = 1:length(colorTypes),
  ctype = colorTypes{i};
  typeIndex = strmatch(ctype,{H.cs(:).colorType},'exact');
  ct = H.cs(typeIndex);
  
  switch ctype,
   case 'div',
    label = 'diverging';
   case 'seq',
    label = 'sequential';
   case 'qual',
    label = 'qualitative';
  end
  H.menu.(ctype) = uimenu(H.menu.select,'Label',label);
  
  % loop over the names
  colorNames = unique({ct(:).colorName});
  for j = 1:length(colorNames),
    
    cname = colorNames{j};
    
    callback = [...
        'CB = get(gcbf,''UserData''); ',...
        'CB.ctype = ''', ctype, ''';',...
        'CB.cname = ''', cname, ''';',...
        'set(gcbf,''UserData'',CB);',...
        'colorBrewPlot(''plot'');'];
    
    H.menu.(cname) = uimenu(H.menu.(ctype),'Label',cname,...
                            'Callback', callback);
    
  end
end




%H.menu.plot = uimenu(H.gui,'Label','Plot',...
%                     'Callback','colorBrewPlot(''plot'');');

H.menu.exit = uimenu(H.gui,'Label','Exit',...
                     'Callback','colorBrewPlot(''exit'');');

% -- help menu

H.menu.Help = uimenu(H.gui,'Label','Help');
H.menu.help = uimenu(H.menu.Help,'Label','Help',...
                     'Callback','doc colorBrewPlot;');
H.menu.gpl  = uimenu(H.menu.Help,'Label','License',...
                     'Callback','colorBrewLicense;');

return
