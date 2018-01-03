function [hReg,SPM,VOL,xX,xCon,xSDM]=spm_tbx_roi
% SPM Diffusion tensor analysis toolbox
% Russ Poldrack, 8/15/00
% - provides preprocessing and visualization tools for diffusion
%   tensor MRI data
% 
% CVS repository and development information can be found at:
% http://spm-toolbox.sourceforge.net 
%
% Options:
% - Combine interleaved series (spm_diffusion_combine.m) 
%	Combines two series of Siemens .ima images acquired separately, 
%	in which interslice interval = slice thickness and slices are 
%	shifted by slice thickness between the two scans
% - Convert diffusion images (spm_diffusion_convert.m)
%	Converts Siemens .ima images for 6-axis diffusion tensor 
%	sequence into Analyze-format (.img) files 
%	NOTE: Ordering of images is sequence specific; this module
%	is specifically written to convert data acquired using the
%	DTI sequence written by Jim Moore at Siemens and currently
%	in use at the MGH-NMR Center
% - Calculate tensor (spm_diffusion_calc.m)
%	Calculates the diffusion tensor and related data based upon
%	a set of .img files representing diffusion along each of 6 
%	directions along with a low-b (non-diffusion-weighted) image.
%	NOTE: As with the .ima conversion routines, this routine
%	assumes the specific gradient combinations and orderings used
%	in the sequuence mentioned above.  These are described more 
%	fully in spm_diffusion_calc.m
% - Display colormap orthviews (spm_diffusion_orthviews2.m)
%	Display parallel orthogonal projections of color-coded 
%	tensor direction maps and anatomical image
% - Display vector map (spm_diffusion_quiv.m)
%	Display overlay of vector map (representing a specified
%	eigenvector) on an anatomical image
%
% The SPM diffusion toolbox is distributed under the GNU Public License.
% Please see spm_LICENSE.man in the SPM distribution for full details.
%

global BCH; %- used as a flag to know if we are in batch mode or not.
SCCSid  = '2.31';


%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Diffusion analysis',0);
spm_clf(Fgraph);
spm_clf(Finter);

spm_help('!ContextHelp',mfilename)

%hReg = uicontrol(Fgraph,'Style','Text','String','hReg',...
%  	'Position',[100 200 100 025],...
%  	'FontName','Times','FontSize',14,'FontWeight','Bold',...
%  	'HorizontalAlignment','Center');

% get design matrix and/or data
%=======================================================================
MType = {'Combine interleaved series (.ima)',...
	 'Convert diffusion images',...
           'Calculate tensor',...
	   'View data using rip',...
	   'Display colormap orthviews',...
           'Display vector map' ...
};
MT    = spm_input('What would you like to do?',1,'m',MType,...
                  'batch',{},'types');

%-Initialise output arguments in case return early
xX   = [];
Sess = [];

switch MT
%-----------------------------------------------------------------------

	case 1
	% combine ima's interleaved across two series
	%---------------------------------------------------------------
	spm_clf(Finter);
        imgs1 = spm_get([1,10000],'.ima','select .imas for first series');
        imgs2 = spm_get([1,10000],'.ima','select .imas for second series');
        spm_diffusion_combine(imgs1,imgs2);
	return

	case 2
	% convert .ima's
	%---------------------------------------------------------------
	spm_clf(Finter);
	%nslices=spm_input('How many slices','+1');
	%nacq=spm_input('How many acquisitions','+1');
	imgs   = spm_get([1,10000],'.ima','select .ima images');
        [ima_dir fname ext]=fileparts(imgs(1,:));
        first_img=sprintf('%s%s',fname,ext);
        spm_diffusion_convert(first_img,1);
	return


	case 3
	% Calculate diffusion
	%---------------------------------------------------------------
	spm_clf(Finter);
	bvalue=spm_input('B value','+1');
	imgs   = spm_get([1,10000],'avg.img','select diffusion images');
        spm_diffusion_calc(imgs,bvalue);

        return

 	case 4
	% Display tensor data using rip (outside of SPM)
	%---------------------------------------------------------------
	spm_clf(Finter);
	load_tensor_data
	return

  	case 5
	% Display colormap and anatomy in register
	%---------------------------------------------------------------
	spm_clf(Finter);
	a_imgs   = spm_get(1,'.img','select anatomical image');
        c_imgs   = spm_get(3,'*col*.img','select diffusion colormaps 1-3');
        spm_diffusion_orthviews2('Image',a_imgs(1,:),c_imgs(1,:),c_imgs(2,:),c_imgs(3,:));

      	return

  	case 6
	% Display vector map
	%---------------------------------------------------------------
	spm_clf(Finter);
	a_img   = spm_get(1,'.img','select anatomical image');
        fa_img=spm_get(1,'.img','select an FA image');
        e_imgs   = spm_get(3,'.img','select eigenvector images X, Y, & Z');
        spm_diffusion_quiv(a_img(1,:),fa_img(1,:),e_imgs(1,:),e_imgs(2,:),e_imgs(3,:));


end












