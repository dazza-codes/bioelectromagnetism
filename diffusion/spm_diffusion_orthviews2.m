function varargout = spm_diffusion_orthviews2(action,varargin)
% spm_diffusion_orthviews2: part of the SPM diffusion tensor toolbox
% Russ Poldrack, 11/19/00
%
% - This routine allows display of parallel orthogonal sections
%   with both color-coded diffusion direction maps, and a 
%   colocalized anatomical image.  
%   It is an adapted version of spm_orthviews, which is used for 
%   orthogonal section display with anatomical/functional data.  
%
% Arguments:
% spm_diffusion_orthviews2 uses the action string method.  See the 
% help for spm_orthviews for more information on action strings not
% described here.  I can't guarantee that all will work, but most of 
% them should - please submit a bug report if you find an action string
% that doesn't work.
%
% The primary method for visualization should use the "Image" action
% string, as follows:
%
% spm_diffusion_orthviews2('Image',anat_img,col1_img,col2_img,col3_img)
%
% anat_img: anatomical image for co-localization
% col[1-3]_img: three images representing the red, green, and blue 
% components of the color map (floating point images with values between 
% 0 and 1).  These images are created by spm_diffusion_calc.
%
% Effects:
% Displays an interactive colormap alongside the specified anatomical image.
%
% CVS repository and development information can be found at:
% http://spm-toolbox.sourceforge.net 
%
% The SPM diffusion toolbox is distributed freely under the GNU Public License.
%


global st col1;

if isempty(st), reset_st; end;

spm('Pointer','watch');

if nargin == 0, action = ''; end;
action = lower(action);


switch lower(action),

 case 'image',

        if (length(varargin)<4),
	  fprintf('Must specify four images!\n');
	  return;
	end;
	clf;
	H = specify_color_image(varargin{1},varargin{2},varargin{3}, ...
				varargin{4});

	if ~isempty(H)
	        st.colorhandle=H;
		if length(varargin)>=5, 
		  st.vols{H}.area = varargin{5}; 
		else
		  st.vols{H}.area=[0.0 0.55 1 0.45];
		end;
		if isempty(st.bb), st.bb = maxbb; end;
		bbox;
		redraw(H);
	end;
	varargout{1} = H;
	H = specify_color_image(varargin{1},varargin{2},varargin{3},varargin{4});
	if ~isempty(H)
		if length(varargin)>=5, 
		  st.vols{H}.area = varargin{5}; 
		else
		  st.vols{H}.area=[0.0 0.05 1 0.45];
		end;
		if isempty(st.bb), st.bb = maxbb; end;
		bbox;
		redraw(H);
	end;
	varargout{1} = H;
%	spm_orthviews('image',varargin{1},[0.0 0.05 1 0.45]);

case 'bb',
	if length(varargin)> 0 & all(size(varargin{1})==[2 3]), st.bb = varargin{1}; end;
	bbox;
	redraw_all;

case 'redraw',
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;

case 'reposition',

	if length(varargin)<1, tmp = findcent;
	else, tmp = varargin{1}; end;
	if length(tmp)==3, st.centre = tmp(:); end;
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;

case 'setcoords',
        fprintf('callback to setcoords:\n');
	st.centre = varargin{1};
	st.centre = st.centre(:);
	redraw_all;
	eval(st.callback);

case 'space',
	if length(varargin)<1,
		st.Space = eye(4);
		st.bb = maxbb;
		redraw_all;
	else,
		space(varargin{1});
		redraw_all;
	end;
	bbox;

case 'maxbb',
	st.bb = maxbb;
	bbox;
	redraw_all;

case 'resolution',
	resolution(varargin{1});
	bbox;
	redraw_all;

case 'window',
	if length(varargin)<2,
		win = 'auto';
	elseif length(varargin{2})==2,
		win = varargin{2};
	end;
	for i=valid_handles(varargin{1}),
		st.vols{i}.window = win;
	end;
	redraw(varargin{1});

case 'delete',
	my_delete(varargin{1});

case 'move',
	move(varargin{1},varargin{2});
	% redraw_all;

case 'reset',
	my_reset;

case 'pos',
	if isempty(varargin),
		H = st.centre(:);
	else,
		H = pos(varargin{1});
	end;
	varargout{1} = H;

case 'interp',
	st.hld = varargin{1};
	redraw_all;

case 'xhairs',
	xhairs(varargin{1});

case 'register',
	register(varargin{1});

case 'addimage',
	addimage(varargin{1}, varargin{2});
	redraw(varargin{1});

case 'addcolouredimage',
	addcolouredimage(varargin{1}, varargin{2},varargin{3});

otherwise,
	warning('Unknown action string')
end;

spm('Pointer');
return;


%_______________________________________________________________________
%_______________________________________________________________________
function addimage(handle, fname)
global st
for i=valid_handles(handle),
	vol = spm_vol(fname);
	mat = vol.mat;
	st.vols{i}.blobs=cell(1,1);
	if st.mode == 0,
		axpos = get(st.vols{i}.ax{2}.ax,'Position');
	else,
		axpos = get(st.vols{i}.ax{1}.ax,'Position');
	end;
	ax = axes('Parent',st.fig,'Position',[(axpos(1)+axpos(3)+0.05) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
		'Box','on');
	mx = max([eps maxval(vol)]);
	image([0 mx/32],[0 mx],[1:64]' + 64,'Parent',ax);
	set(ax,'YDir','normal','XTickLabel',[]);
	st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'cbar',ax,'max',mx);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredimage(handle, fname,colour)
global st
for i=valid_handles(handle),

	vol = spm_vol(fname);
	mat = vol.mat;
	if ~isfield(st.vols{i},'blobs'),
		st.vols{i}.blobs=cell(1,1);
		bset = 1;
	else,
		bset = length(st.vols{i}.blobs)+1;
	end;
	axpos = get(st.vols{i}.ax{2}.ax,'Position');
	mx = max([eps maxval(vol)]);
	st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'colour',colour);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function register(hreg)
global st
tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles(1:24);
if ~isempty(h),
	tmp = st.vols{h(1)}.ax{1}.ax;
	st.registry = struct('hReg',hreg,'hMe', tmp);
	spm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'spm_diffusion_orthviews2');
else,
	warning('Nothing to register with');
end;
st.centre = spm_XYZreg('GetCoords',st.registry.hReg);
st.centre = st.centre(:);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function xhairs(arg1),
global st
st.xhairs = 0;
opt = 'on';
if ~strcmp(arg1,'on'),
	opt = 'off';
else,
	st.xhairs = 1;
end;
for i=valid_handles(1:24),
	for j=1:3,
		set(st.vols{i}.ax{j}.lx,'Visible',opt); 
		set(st.vols{i}.ax{j}.ly,'Visible',opt);  
	end; 
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = pos(arg1)
global st
H = [];
for arg1=valid_handles(arg1),
	is = inv(st.vols{arg1}.premul*st.vols{arg1}.mat);
	H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_reset
global st
if ~isempty(st) & isfield(st,'registry') & ishandle(st.registry.hMe),
	delete(st.registry.hMe); st = rmfield(st,'registry');
end;
my_delete(1:24);
reset_st;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_delete(arg1)
global st

for i=valid_handles(arg1),
	kids = get(st.fig,'Children');
	for j=1:3,
		if any(kids == st.vols{i}.ax{j}.ax),
			set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
			delete(st.vols{i}.ax{j}.ax);
		end;
	end;
	st.vols{i} = [];
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function resolution(arg1)
global st
res      = arg1/mean(svd(st.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st.Space = st.Space*Mat;
st.bb    = st.bb/res;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function move(handle,pos)
global st
for handle = valid_handles(handle),
	st.vols{handle}.area = pos;
end;
bbox;
% redraw(valid_handles(handle));
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bb = maxbb
global st
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles(1:24),
	bb = [[1 1 1];st.vols{i}.dim(1:3)];
	c = [	bb(1,1) bb(1,2) bb(1,3) 1
		bb(1,1) bb(1,2) bb(2,3) 1
		bb(1,1) bb(2,2) bb(1,3) 1
		bb(1,1) bb(2,2) bb(2,3) 1
		bb(2,1) bb(1,2) bb(1,3) 1
		bb(2,1) bb(1,2) bb(2,3) 1
		bb(2,1) bb(2,2) bb(1,3) 1
		bb(2,1) bb(2,2) bb(2,3) 1]';
	tc = st.Space\(st.vols{i}.premul*st.vols{i}.mat)*c;
	tc = tc(1:3,:)';
	mx = max([tc ; mx]);
	mn = min([tc ; mn]);
end;
bb = [mn ; mx];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function space(arg1)
global st
if ~isempty(st.vols{arg1})
	num = arg1;
	Mat = st.vols{num}.premul(1:3,1:3)*st.vols{num}.mat(1:3,1:3);
	Mat = diag([sqrt(sum(Mat.^2)) 1]);
	Space = (st.vols{num}.premul*st.vols{num}.mat)/Mat;
	bb = [1 1 1;st.vols{num}.dim(1:3)];
	bb = [bb [1;1]];
	bb=bb*Mat';
	bb=bb(:,1:3);
	st.Space  = Space;
	st.bb = bb;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = specify_color_image(arg1,arg2,arg3,arg4)
global st
H=[];
ok = 1;
eval('V = spm_vol(arg1);','ok=0;');
if ok == 0,
	fprintf('Can not use image "%s"\n', arg1);
	return;
end;

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end;

DeleteFcn = ['spm_diffusion_orthviews2(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3,
	ax = axes('Visible','off','DrawMode','fast','Parent',st.fig,'DeleteFcn',DeleteFcn,...
		'YDir','normal');
	d  = image(0,'Tag','Transverse','Parent',ax,...
		'DeleteFcn',DeleteFcn);
	set(ax,'Ydir','normal');
	lx = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	ly = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	if ~st.xhairs,
		set(lx,'Visible','off');
		set(ly,'Visible','off');
	end;
	V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';
st.vols{ii} = V;


% repeat this for each of the three colormaps - each in a separate img
% colormaps are contained in st.cmap1-3

ok = 1;
eval('V = spm_vol(arg2);','ok=0;');
if ok == 0,
	fprintf('Can not use image "%s"\n', arg2);
	return;
end;

V.ax = cell(3,1);
for i=1:3,

  V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';

st.cmap1{ii}=V;

ok = 1;
eval('V = spm_vol(arg3);','ok=0;');
if ok == 0,
	fprintf('Can not use image "%s"\n', arg3);
	return;
end;

V.ax = cell(3,1);
for i=1:3,

	V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';

st.cmap2{ii}=V;


ok = 1;
eval('V = spm_vol(arg4);','ok=0;');
if ok == 0,
	fprintf('Can not use image "%s"\n', arg4);
	return;
end;

V.ax = cell(3,1);
for i=1:3,

	V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';

st.cmap3{ii}=V;


H = ii;
return;

%_______________________________________________________________________
%_______________________________________________________________________

function H = specify_image(arg1)
global st
H=[];
ok = 1;
eval('V = spm_vol(arg1);','ok=0;');
if ok == 0,
	fprintf('Can not use image "%s"\n', arg1);
	return;
end;

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end;

DeleteFcn = ['spm_diffusion_orthviews2(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3,
	ax = axes('Visible','off','DrawMode','fast','Parent',st.fig,'DeleteFcn',DeleteFcn,...
		'YDir','normal');
	d  = image(0,'Tag','Transverse','Parent',ax,...
		'DeleteFcn',DeleteFcn);
	set(ax,'Ydir','normal');
	lx = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	ly = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	if ~st.xhairs,
		set(lx,'Visible','off');
		set(ly,'Visible','off');
	end;
	V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';
st.vols{ii} = V;


H = ii;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bbox
global st
Dims = diff(st.bb)'+1;

TD = Dims([1 2])';
CD = Dims([1 3])';
if st.mode == 0, SD = Dims([3 2])'; else, SD = Dims([2 3])'; end;

un    = get(st.fig,'Units');set(st.fig,'Units','Pixels');sz=get(st.fig,'Position');set(st.fig,'Units',un);
sz    = sz(3:4);
sz(2) = sz(2)-40;

for i=valid_handles(1:24),
	area = st.vols{i}.area(:);
	area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
	if st.mode == 0,
		sx   = area(3)/(Dims(1)+Dims(3))/1.02;
	else,
		sx   = area(3)/(Dims(1)+Dims(2))/1.02;
	end;
	sy   = area(4)/(Dims(2)+Dims(3))/1.02;
	s    = min([sx sy]);

	offy = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
	sky = s*(Dims(2)+Dims(3))*0.02;
	if st.mode == 0,
		offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(3))*0.02;
	else,
		offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(2))*0.02;
	end;

	DeleteFcn = ['spm_diffusion_orthviews2(''Delete'',' num2str(i) ');'];

	% Transverse
	set(st.vols{i}.ax{1}.ax,'Units','pixels', ...
		'Position',[offx offy s*Dims(1) s*Dims(2)],...
		'Units','normalized','Xlim',[0 TD(1)]+0.5,'Ylim',[0 TD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Coronal
	set(st.vols{i}.ax{2}.ax,'Units','Pixels',...
		'Position',[offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)],...
		'Units','normalized','Xlim',[0 CD(1)]+0.5,'Ylim',[0 CD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Saggital
	if st.mode == 0,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
			'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	else,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
			'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function redraw_all
global st
redraw(1:24);
return;
%_______________________________________________________________________
function mx = maxval(vol)
if isstruct(vol),
	mx = -Inf;
	for i=1:vol.dim(3),
		tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
		mx  = max(mx,max(tmp(:)));
	end;
else,
	mx = max(vol(:));
end;

%_______________________________________________________________________
function redraw(arg1)
global st
bb   = st.bb;
Dims = diff(bb)';
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);

for i = valid_handles(arg1),
	M = st.vols{i}.premul*st.vols{i}.mat;
	TM0 = [	1 0 0 -bb(1,1)
		0 1 0 -bb(1,2)
		0 0 1 -cent(3)
		0 0 0 1];
	TM = inv(TM0*(st.Space\M));
	TD = Dims([1 2]);

	CM0 = [	1 0 0 -bb(1,1)
		0 0 1 -bb(1,3)
		0 1 0 -cent(2)
		0 0 0 1];
	CM = inv(CM0*(st.Space\M));
	CD = Dims([1 3]);

	if st.mode ==0,
		SM0 = [	0 0 1 -bb(1,3)
			0 1 0 -bb(1,2)
			1 0 0 -cent(1)
			0 0 0 1];
		SM = inv(SM0*(st.Space\M)); SD = Dims([3 2]);
	else,
		SM0 = [	0  1 0 -bb(1,2)
			0  0 1 -bb(1,3)
			1  0 0 -cent(1)
			0  0 0 1];
		SM0 = [	0 -1 0 +bb(2,2)
			0  0 1 -bb(1,3)
			1  0 0 -cent(1)
			0  0 0 1];
		SM = inv(SM0*(st.Space\M));
		SD = Dims([2 3]);
	end;

	ok=1;
	eval('imgt  = (spm_slice_vol(st.vols{i},TM,TD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);end;
	eval('imgt_c1  = (spm_slice_vol(st.cmap1{i},TM,TD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap1{i}.fname);end;
	eval('imgt_c2  = (spm_slice_vol(st.cmap2{i},TM,TD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap2{i}.fname);end;
	eval('imgt_c3  = (spm_slice_vol(st.cmap3{i},TM,TD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap3{i}.fname);end;
	eval('imgc  = (spm_slice_vol(st.vols{i},CM,CD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);end;
	eval('imgc_c1  = (spm_slice_vol(st.cmap1{i},CM,CD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap1{i}.fname);end;
	eval('imgc_c2  = (spm_slice_vol(st.cmap2{i},CM,CD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap2{i}.fname);end;
	eval('imgc_c3  = (spm_slice_vol(st.cmap3{i},CM,CD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap3{i}.fname);end;
	eval('imgs  = (spm_slice_vol(st.vols{i},SM,SD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);end;
	eval('imgs_c1  = (spm_slice_vol(st.cmap1{i},SM,SD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap1{i}.fname);end;
	eval('imgs_c2  = (spm_slice_vol(st.cmap2{i},SM,SD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap2{i}.fname);end;
	eval('imgs_c3  = (spm_slice_vol(st.cmap3{i},SM,SD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.cmap3{i}.fname);end;

        if (ok ~= 0)

          if (i == st.colorhandle),
	  % create color maps from individual color images
	       col_t=zeros([size(imgt_c1) 3]);
	       col_c=zeros([size(imgc_c1) 3]);
	       col_s=zeros([size(imgs_c1) 3]);
	       imgt(find(imgt>1))=1;
	       imgc(find(imgc>1))=1;
	       imgs(find(imgs>1))=1;
	       fa_thresh=0.25;
	       imgt(find(imgt<fa_thresh))=0;
	       imgc(find(imgc<fa_thresh))=0;
	       imgs(find(imgs<fa_thresh))=0;

               col_t(:,:,1)=imgt_c1.*imgt;
	       col_t(:,:,2)=imgt_c2.*imgt;
	       col_t(:,:,3)=imgt_c3.*imgt;
	       col_c(:,:,1)=imgc_c1.*imgc;
	       col_c(:,:,2)=imgc_c2.*imgc;
	       col_c(:,:,3)=imgc_c3.*imgc;
	       col_s(:,:,1)=imgs_c1.*imgs;
	       col_s(:,:,2)=imgs_c2.*imgs;
	       col_s(:,:,3)=imgs_c3.*imgs;

	       col_t(find(col_t < 0))=0.0;
	       col_c(find(col_c < 0))=0.0;
	       col_s(find(col_s < 0))=0.0;
	     
	     imgt=col_t;
	     imgc=col_c;
	     imgs=col_s;
	   else,
	            if strcmp(st.vols{i}.window,'auto'),
			mx = -Inf; mn = Inf;
			if ~isempty(imgt),
				mx = max([mx max(max(imgt))]);
				mn = min([mn min(min(imgt))]);
			end;
			if ~isempty(imgc),
				mx = max([mx max(max(imgc))]);
				mn = min([mn min(min(imgc))]);
			end;
			if ~isempty(imgs),
				mx = max([mx max(max(imgs))]);
				mn = min([mn min(min(imgs))]);
			end;
			if mx==mn, mx=mn+eps; end;
		else,
			mx = st.vols{i}.window(2);
			mn = st.vols{i}.window(1);
			r=min([mn mx]);imgt = max(imgt,r); r=max([mn mx]);imgt = min(imgt,r);
			r=min([mn mx]);imgc = max(imgc,r); r=max([mn mx]);imgc = min(imgc,r);
			r=min([mn mx]);imgs = max(imgs,r); r=max([mn mx]);imgs = min(imgs,r);
		end;
	    	scal = 64/(mx-mn);
		dcoff = -mn*scal;
		imgt = imgt*scal+dcoff;
		imgc = imgc*scal+dcoff;
		imgs = imgs*scal+dcoff;

	    
	   end;

% put data from separate color images into single array



		callback = 'spm_diffusion_orthviews2(''Reposition'');';
		set(st.vols{i}.ax{1}.d,'ButtonDownFcn',callback, 'Cdata',imgt);
%		set(st.vols{i}.ax{1}.d,'ButtonDownFcn',callback, 'Cdata',col_t);
		set(st.vols{i}.ax{1}.lx,'ButtonDownFcn',callback,...
			'Xdata',[0 TD(1)],'Ydata',[1 1]*(cent(2)-bb(1,2)));
		set(st.vols{i}.ax{1}.ly,'ButtonDownFcn',callback,...
			'Ydata',[0 TD(2)],'Xdata',[1 1]*(cent(1)-bb(1,1)));

		set(st.vols{i}.ax{2}.d,'ButtonDownFcn',callback, 'Cdata',imgc);
%		set(st.vols{i}.ax{2}.d,'ButtonDownFcn',callback, 'Cdata',col_c);
		set(st.vols{i}.ax{2}.lx,'ButtonDownFcn',callback,...
			'Xdata',[0 CD(1)],'Ydata',[1 1]*(cent(3)-bb(1,3)));
		set(st.vols{i}.ax{2}.ly,'ButtonDownFcn',callback,...
			'Ydata',[0 CD(2)],'Xdata',[1 1]*(cent(1)-bb(1,1)));

		set(st.vols{i}.ax{3}.d,'ButtonDownFcn',callback,'Cdata',imgs);
%		set(st.vols{i}.ax{3}.d,'ButtonDownFcn',callback,'Cdata',col_s);
	        if st.mode ==0,
			set(st.vols{i}.ax{3}.lx,'ButtonDownFcn',callback,...
				'Xdata',[0 SD(1)],'Ydata',[1 1]*(cent(2)-bb(1,2)));
			set(st.vols{i}.ax{3}.ly,'ButtonDownFcn',callback,...
				'Ydata',[0 SD(2)],'Xdata',[1 1]*(cent(3)-bb(1,3)));
		else,
			set(st.vols{i}.ax{3}.lx,'ButtonDownFcn',callback,...
				'Xdata',[0 SD(1)],'Ydata',[1 1]*(cent(3)-bb(1,3)));
			set(st.vols{i}.ax{3}.ly,'ButtonDownFcn',callback,...
				'Ydata',[0 SD(2)],'Xdata',[1 1]*(bb(2,2)-cent(2)));
		end;
	end; % if ok~=0
end; % for handles
drawnow;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function centre = findcent
global st
obj    = get(st.fig,'CurrentObject');

centre = [];
cent   = [];
cp     = [];

for i=valid_handles(1:24),
    fprintf('doing handle %d\n',i);
	for j=1:3,

		if ~isempty(obj),
			if any([st.vols{i}.ax{j}.d  ...
				st.vols{i}.ax{j}.lx ...
				st.vols{i}.ax{j}.ly]== obj)
				cp = get(get(obj,'Parent'),'CurrentPoint');
			elseif (st.vols{i}.ax{j}.ax == obj),
				cp = get(obj,'CurrentPoint');
			end;
		end;


		if ~isempty(cp),
			cp   = cp(1,1:2);
			is   = inv(st.Space);
			cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
			switch j,
				case 1,
				cent([1 2])=[cp(1)+st.bb(1,1) cp(2)+st.bb(1,2)];
				case 2,
				cent([1 3])=[cp(1)+st.bb(1,1) cp(2)+st.bb(1,3)];
				case 3,
				if st.mode ==0,
					cent([3 2])=[cp(1)+st.bb(1,3) cp(2)+st.bb(1,2)];
				else,
					cent([2 3])=[st.bb(2,2)-cp(1) cp(2)+st.bb(1,3)];
				end;
			end;
			break;
		end;
	end;

	if ~isempty(cent), break; end;
end;

if ~isempty(cent), centre = st.Space(1:3,1:3)*cent(:) + st.Space(1:3,4); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function handles = valid_handles(handles)
global st;
handles = handles(:)';
handles = handles(find(handles<=24 & handles>=1 & ~rem(handles,1)));
for h=handles,
	if isempty(st.vols{h}), handles(find(handles==h))=[]; 
%	  fprintf('handle %d is empty\n',h);
	else
%	  fprintf('handle %d is good\n',h);
	  
	end;
end;

return;
%_______________________________________________________________________
%_______________________________________________________________________
function reset_st
global st
fig     = spm_figure('FindWin','Graphics');
bb      = [ [-78 78]' [-112 76]' [-50 85]' ];
st      = struct('n', 0, 'vols',[], 'cmap1',[],'cmap2',[],'cmap3',[],'bb',bb,'Space',eye(4),'centre',[0 0 0],'callback',';','xhairs',1,'hld',1,'fig',fig,'mode',1,'colorhandle',[]);
st.vols = cell(24,1);
return;







