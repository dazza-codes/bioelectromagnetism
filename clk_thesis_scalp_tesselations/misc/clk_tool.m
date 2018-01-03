function clk_tool(varargin)
% CLK_TOOL	Miscellaneous MATLAB graphics tools
%		CLK_TOOL(tool [, ...])
%		    'all'   - all available tools
%		    'none'  - remove all UI controls
%
%		    'axes'  - menu for axes labels
%		    'cmap'  - menu for color maps
%		    'misc'  - menu for misc actions
%		    'shade' - menu for shading style
%		    'shift' - menu for shifting axis
%		    'zoom'  - menu for zoom modes
%		    'spin'  - sliders for azimuth and elevation

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% functions defined at end of file: 
%   ax_kids, cb_kids, uicntl_kids, uimenu_kids, is2dplot, is3dplot

%%% THINGS TO DO
% ? default to black lines / text
% ? optional argument to select figure / subplot / axes
% ? text options: (left|center|right, top|middle|bottom, normal|italic|oblique)
% ? create other (filled?) shapes - square, circle, pgon, points
% ? menu(s) to set color/style of points/lines?
% ? extend shade to allow edges w/o faces, points only, etc.
% ? use set(gcf,'KeyPressFcn') for accelerator keys

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin == 0) help clk_tool; return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% step through each argument and perform actions

for k = 1:nargin
    tool = varargin{k};    
    if ~ischar(tool) error(['Invalid argument.']); end
    
    switch tool
	case 'all',
	    if is2dplot
		clk_tool  axes cmap misc shade shift zoom;
	    else
		clk_tool  axes cmap misc shade spin  zoom;
	    end

    	case 'list',
    	    disp(get(uicntl_kids, 'String'));
    	    disp(get(uimenu_kids, 'Label'));
    	    
    	case 'new',
    	    clk_tool none all;

	case 'none',
	    delete([uicntl_kids, uimenu_kids]);

    	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% menu for axis labels in figure

    	case 'axes',
    	    clk_tool noaxes;
    	    mh = clk_tool_uimenu('&Axes', 'clk_tool cbaxes', ...
    	    	    	    	{'u1,u2','f1,f2','x,y,z','X,Y,Z', ...
    	    	    	    	'(none)','OFF','ON'});
    	    set(findobj(mh,'Label','(none)'),'Checked','on');
    	    	
    	case 'cbaxes',
    	    % clear and set check mark
    	    set(findobj(get(gcbo,'Parent'),'Checked','on'),'Checked','off');
    	    set(gcbo,'Checked','on');
    	    switch get(gcbo,'Label')
    	    	case  'u1,u2', xval='u1'; yval='u2'; zval='';
    	    	case  'f1,f2', xval='f1'; yval='f2'; zval='';
    	    	case  'x,y,z', xval='x';  yval='y';  zval='z';
    	    	case  'X,Y,Z', xval='X';  yval='Y';  zval='Z';
    	    	case '(none)', xval='';   yval='';   zval='';
    	    end % switch
	    for j = ax_kids
	    	subplot(j);
		switch get(gcbo,'Label')
		    case 'OFF', axis off;
		    case 'ON',  axis on;
		    otherwise, xlabel(xval); ylabel(yval); zlabel(zval);
		end % switch
	    end % for
    	    
    	case 'noaxes',
    	    delete(findobj(uimenu_kids, 'Label', '&Axes'));


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% menu for colormap in figure

    	case 'cmap',
    	    clk_tool nocmap;
    	    mh = clk_tool_uimenu('&CMap', 'clk_tool cbcmap', ...
    	    	    	    	{'hsv', 'hot', 'gray', 'bone', 'copper', 'pink', 'white', ...
    	    	     	    	 'flag', 'lines', 'colorcube', 'jet', 'prism', 'cool', ...
    	    	     	    	 'autumn', 'spring', 'winter', 'summer', ...
    	    	     	    	 'gray(4)', 'gray(16)', 'gray(64)','gray(256)'});
    	    set(findobj(mh,'Label','jet'),'Checked','on');
    	
    	case 'cbcmap',
    	    % set check mark for menu item
    	    set(findobj(get(gcbo,'Parent'),'Checked','on'),'Checked','off');
    	    set(gcbo,'Checked','on');
    	    colormap(get(gcbo,'Label'));
    	    for j = cb_kids, colorbar(j); end
    	
    	case 'nocmap',
    	    delete(findobj(uimenu_kids, 'Label', '&CMap'));


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% menu for miscellaneous actions

    	case 'misc',
    	    clk_tool nomisc;
    	    mh = clk_tool_uimenu('&Misc', 'clk_tool cbmisc', ...
    	    	    	    	{'apply...', 'axis...', 'caxis...', ...
    	    	    	    	 'brighten', 'darken', 'EPSF C', 'EPSF BW', ...
    	    	     	    	 'ginput', 'gline', 'gtext', 'grid on', 'grid off', ...
    	    	     	    	 {'toggle', 'image', 'line', 'patch', 'surface', 'text'}, ...
    	    	     	    	 'size*2', 'size/2', 'no tools'});
    	    set(findobj(mh, 'Label','brighten'), 'Accelerator', 'b');
    	    set(findobj(mh, 'Label','darken'),   'Accelerator', 'd');
    	    
	case 'cbmisc',
	    switch get(gcbo, 'Label')
	    	case 'apply...',
	    	    beep;
	    	    val = input('Command to apply to all axes in figure: ','s');
	    	    for j = ax_kids, subplot(j); eval(val); end
	    	
	    	case 'axis...',
	    	    beep(2); val = input('Enter [xmin xmax ymin ymax zmin zmax]: ');
	    	    for j = ax_kids, subplot(j); axis(val); end
	    	case 'caxis...',
	    	    beep(2); val = input('Enter [cmin cmax]: ');
	    	    for j = ax_kids, subplot(j); caxis(val); end
    	    	    for j = cb_kids, colorbar(j); end
	    	    
		case 'brighten',    brighten( 0.25);
		case 'darken',	    brighten(-0.25);
		
		case 'EPSF C',
		    disp('Printing to color EPS file ''clk_out.eps''.');
		    print -noui -depsc clk_out.eps
		case 'EPSF BW',
		    disp('Printing to B/W EPS file ''clk_out.eps''.');
		    print -noui -deps clk_out.eps

		case 'ginput',
		    disp('Click on desired locations; hit return in figure to stop.');
		    ginput
		case 'gline',
		    disp('Click on vertices with button 1; other button to stop.');
		    ox = NaN; oy = NaN;
		    while (1)
			[x, y, b] = ginput(1);
			line([ox, x], [oy, y]); ox = x; oy = y;
			if b ~= 1 break; end
		    end
		case 'gtext',
		    str = input('Enter text string to place: ', 's');
		    disp('Click on desired location.');
		    [x, y] = ginput(1); text(x, y, str);

		case 'grid on',
		    set(ax_kids, 'Xgrid','on',  'Ygrid','on',  'Zgrid','on');
		case 'grid off',
		    set(ax_kids, 'Xgrid','off', 'Ygrid','off', 'Zgrid','off');

    	    	case {'image','line','patch','surface','text'},
    	    	    val = get(gcbo, 'Label');
    	    	    kkids = findobj(ax_kids, 'Type', val);
    	    	    if ~isempty(kkids)
    	    	    	vkkids = strcmp(get(kkids, 'Visible'), 'on');
    	    	    	set(kkids( vkkids), 'Visible','off');
    	    	    	set(kkids(~vkkids), 'Visible','on');
    	    	    end
    	    	        	    	    
		case 'size*2',
		    set(gcbf, 'Position', [1 1 2.0 2.0] .* get(gcbf, 'Position'));
		case 'size/2',
		    set(gcbf, 'Position', [1 1 0.5 0.5] .* get(gcbf, 'Position'));

		case 'no tools',
		    clk_tool none
	    end

	case 'nomisc',
    	    delete(findobj(uimenu_kids, 'Label','&Misc'));


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% menu for axis shadings in figure

	case 'shade',
	    clk_tool noshade;
	    mh = clk_tool_uimenu('&Shade', 'clk_tool cbshade', ...
    	    	    	    	{{'shading', ' faceted', ' flat', ' interp'}, ...
    	    	     	    	 {'lighting', ' none', ' flat', ' gouraud', ' phong'}, ...
    	    	     	    	 {'material', ' default', ' shiny', ' dull', ' metal'}});
    	    set(findobj(findobj(mh, 'Label','shading' ), 'Label',' faceted'), 'Checked','on');
    	    set(findobj(findobj(mh, 'Label','lighting'), 'Label',' flat'   ), 'Checked','on');
    	    set(findobj(findobj(mh, 'Label','material'), 'Label',' default'), 'Checked','on');

	case 'cbshade',
    	    % set check mark for menu item
    	    set(findobj(get(gcbo,'Parent'), 'Checked','on'), 'Checked','off');
    	    set(gcbo,'Checked','on');
    	    val = [get(get(gcbo,'Parent'), 'Label'), get(gcbo, 'Label')];
	    for j = ax_kids, subplot(j); eval(val); end

	case 'noshade',
    	    delete(findobj(uimenu_kids, 'Label','&Shade'));
   	    	    	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% menu for axis shift control in figure

    	case 'shift',
	    clk_tool noshift;
	    mh = clk_tool_uimenu('&Shift', 'clk_tool cbshift', ...
	    	    	    	{'left', 'down', 'up', 'right'});
	    % use VI-style controllers 
	    set(findobj(mh, 'Label',  'left'), 'Accelerator', 'h');
	    set(findobj(mh, 'Label',  'down'), 'Accelerator', 'j');
	    set(findobj(mh, 'Label',    'up'), 'Accelerator', 'k');
	    set(findobj(mh, 'Label', 'right'), 'Accelerator', 'l');
	
	case 'cbshift',
	    a1 = axis;
	    switch get(gcbo, 'Label')
	    	case  'left', a1(1:2) = a1(1:2) - 0.1 * diff(a1(1:2));
		case  'down', a1(3:4) = a1(3:4) - 0.1 * diff(a1(3:4));	
		case    'up', a1(3:4) = a1(3:4) + 0.1 * diff(a1(3:4));			
		case 'right', a1(1:2) = a1(1:2) + 0.1 * diff(a1(1:2));
    	    end % switch on label
	    axis(a1);
	    
	case 'noshift',
	    delete(findobj(uimenu_kids, 'Label', '&Shift'));
	    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% menu for axis zoom control in figure

	case 'zoom',
	    clk_tool nozoom;
	    if is2dplot
	    	zoom on;
	    	mh = clk_tool_uimenu('&Zoom', 'clk_tool cbzoom2', ...
	    	    	    	    {'on', 'xon', 'yon', 'reset', 'out', 'off'});
    	    	set(findobj(mh,'Label','on'),'Checked','on');
	    	
	    else
		mh = clk_tool_uimenu('&Zoom', 'clk_tool cbzoom3', ...
		    	       {'0.5', '0.707', '1.414', '2.0'});
    		set(findobj(mh, 'Label','0.707'), 'Accelerator', 'z');
    		set(findobj(mh, 'Label','1.414'), 'Accelerator', 'a');
	    end

	case 'cbzoom2',
    	    % set check mark for menu item
    	    set(findobj(get(gcbo,'Parent'), 'Checked','on'), 'Checked','off');
    	    set(gcbo,'Checked','on');
    	    feval('zoom',gcbf,get(gcbo,'Label'));

    	case 'cbzoom3',
	    val = get(gcbo,'Label');
	    for j = ax_kids, subplot(j); feval('zoom',str2num(val)); end

	case 'nozoom',
    	    delete(findobj(uimenu_kids, 'Label','&Zoom'));


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Sliders for axis azimuth and elevation in figure

	case 'spin',
	    clk_tool nospin;
	    % check for 3D plot
	    if is3dplot
		[azi, ele] = view;
		% azimuth from -180 to 179
		uicontrol(  'Style',	'Slider', ...
		    	    'String',	'Azimuth', ...
			    'Units',	'normalized', ...
			    'Position',	[.10, .00, .80, .05], ...
			    'Min',	-180, ...
			    'Max',	 179, ...
			    'SliderStep',[0.01 0.05], ...
			    'Value',	-mod(azi+180,360)+180, ...
			    'Callback',	'clk_tool cbspin_azi' );
		% elevation from -90 to 89
		uicontrol(  'Style',	'Slider', ...
		    	    'String',	'Elevation', ...
			    'Units',	'normalized', ...
			    'Position',	[.00, .10, .05, .80], ...
			    'Min',	-90, ...
			    'Max',	 89, ...
			    'SliderStep',[0.01 0.05], ...
			    'Value',	-mod(ele+90,180)+90, ...
			    'Callback',	'clk_tool cbspin_ele' );
	    else
		error('Spin only works for 3D plots.');
	    end

	case 'cbspin_azi',
	    [azi, ele] = view;
	    for j = ax_kids
		subplot(j); view(-get(gco, 'Value'), mod(ele+180,360)-180); 
	    end

	case 'cbspin_ele',
	    [azi, ele] = view;
	    for j = ax_kids
	    	subplot(j); view(mod(azi+90,180)-90, -get(gco, 'Value'));
	    end

	case 'nospin',
	    delete(findobj(uicntl_kids, 'CallBack','clk_tool cbspin_azi'));
	    delete(findobj(uicntl_kids, 'CallBack','clk_tool cbspin_ele'));


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	otherwise, 
	    error(['Invalid clk_tool control: ', tool]);
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting functions start here
%     AX_KIDS uses COLORBAR method to identify colorbars so we can ignore them
%     (colorbars have handles in UserData)

function mh = clk_tool_uimenu(label,callback,entries)
    mh = uimenu('Label',label);
    for i = entries, i=i{:};
    	% deal with submenus
    	if iscell(i) & (length(i) > 1)
    	    mh1 = uimenu(mh, 'Label',i{1});
    	    for j = i(2:end), j=j{:};
    	    	uimenu(mh1, 'Callback',callback, 'Label',j);
    	    end
    	else
    	    uimenu(mh, 'Callback',callback, 'Label',i);
    	end
    end

function kids = ax_kids
    ax_kids = findobj(gcf, 'Type','axes')';
    mask = ones(size(ax_kids));
    for j = 1:length(ax_kids)
    	ax_d = get(ax_kids(j), 'UserData');
    	if prod(size(ax_d)) == 1 & ishandle(ax_d) mask(j) = 0; end
    end
    kids = ax_kids( logical(mask));
   
function kids = cb_kids
    ax_kids = findobj(gcf, 'Type','axes')';
    mask = ones(size(ax_kids));
    for j = 1:length(ax_kids)
    	ax_d = get(ax_kids(j), 'UserData');
    	if prod(size(ax_d)) == 1 & ishandle(ax_d) mask(j) = 0; end
    end
    kids = ax_kids(~logical(mask));
   
function kids = uicntl_kids
    kids = findobj(get(gcf, 'Children'), 'flat', 'Type','uicontrol')';

function kids = uimenu_kids
    kids = findobj(get(gcf, 'Children'), 'flat', 'Type','uimenu')';

function tf = is2dplot
    tf = all(get(gca, 'view') == [0 90]);
    
function tf = is3dplot
    tf = any(get(gca, 'view') ~= [0 90]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

