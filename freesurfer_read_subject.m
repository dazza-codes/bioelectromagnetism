function [subject] = freesurfer_read_subject(subjectPath,filetype)

% freesurfer_read_subject - Call functions to read subject surface data
%
% [subject] = freesurfer_read_subject(subjectPath,filetype)
%
% This a general wrapper for reading lh/rh surface data.  It calls other
% input functions, such as freesurfer_read_surf.  The main purpose of this
% function is to integrate the surface data into a coherent, consistent
% data structure.
%
% subjectPath is a full path to the subject folder, which contains all the
% freesurfer subdirectories, such as mri, surf, label etc.
%
% filetype is a cell array of strings to indicate what files to read, such
% as 'orig', 'smoothwm', 'inflated', 'white', 'pial', 'curv', 'thickness',
% etc.  The default 'all' will read surfaces and their curvature/thickness'
% that is:
%
% filetype = {'bem','orig','smoothwm','white','pial','inflated',...
%             'curv','thickness','aparc','xfm'}
%
% subject is a data structure with fields for the '.file' paths and the
% data, including .brain and cortical surfaces, the latter are separated
% into .lh and .rh fields and subfields for each file type read.  If the
% talairach.xfm file exists in the <subject>/mri/transforms, the matrix is
% returned into subject.TalairachXFM; to get the Talairach coordinates, use
% freesurfer_surf2tal
%
% See also freesurfer_read_surf, freesurfer_read_curv
%

% $Revision: 1.5 $ $Date: 2007/05/09 00:06:58 $

% Copyright (C) 2000  Darren L. Weber
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% History:  04/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.5 $ $Date: 2007/05/09 00:06:58 $';
fprintf('FREESURFER_READ_SUBJECT [v %s]\n',ver(11:15));

if ~exist('subjectPath','var'),
	error('no subjectPath specified');
end
if isempty(subjectPath),
	error('empty subjectPath specified');
end
if exist(subjectPath) ~= 7,
	error('subjectPath specified is not a directory');
end

surfPath = fullfile(subjectPath,'surf','');
labelPath = fullfile(subjectPath,'label','');
mriPath = fullfile(subjectPath,'mri','');
xfmPath = fullfile(mriPath,'transforms','');
bemPath = fullfile(subjectPath,'bem','');

%subject.avwPath = fullfile(mriPath,'analyze','');


if ~exist('filetype','var'),
	filetype = 'all';
end
if isempty(filetype),
	filetype = 'all';
end

filetype = lower(filetype);
if strmatch('all',filetype),
	filetype = {'bem','orig','smoothwm','white','pial','inflated','curv','thickness','aparc','xfm'};
end

for f = filetype,
	f = f{:};
	switch f,

		case {'bem'},

			% load all the BEM surface files
			d = dir(bemPath);
			for i = 1:length(d),
				if d(i).isdir,
					continue
				end
				if findstr('brain_surface',d(i).name),
					file.brain = fullfile(bemPath,d(i).name);
					[brain.vertices, brain.faces] = freesurfer_read_surf(file.brain);
					subject.brain = brain;
					continue
				end
				if findstr('inner_skull_surface',d(i).name),
					file.inner_skull = fullfile(bemPath,d(i).name);
					[inner_skull.vertices, inner_skull.faces] = freesurfer_read_surf(file.inner_skull);
					subject.inner_skull = inner_skull;
					continue
				end
				if findstr('outer_skull_surface',d(i).name),
					file.outer_skull = fullfile(bemPath,d(i).name);
					[outer_skull.vertices, outer_skull.faces] = freesurfer_read_surf(file.outer_skull);
					subject.outer_skull = outer_skull;
					continue
				end
				if findstr('outer_skin_surface',d(i).name),
					file.outer_skin = fullfile(bemPath,d(i).name);
					[outer_skin.vertices, outer_skin.faces] = freesurfer_read_surf(file.outer_skin);
					subject.outer_skin = outer_skin;
					continue
				end
			end

		case {'brain'},

			% This section loads an older format .tri file.  Recent
			% versions of mri_watershed will output surface files

			% load brain surface
			file.brain = fullfile(bemPath,'brain.tri');
			if exist(file.brain) == 2,
				[brain.vertices, brain.faces] = freesurfer_read_tri(file.brain);
				subject.brain = brain;
			else
				msg = sprintf('file does not exist: %s',file.brain);
				warning(msg);
			end

		case {'orig'},

			% load orig surfaces
			file.lh.orig = fullfile(surfPath,'lh.orig');
			if exist(file.lh.orig) == 2,
				[lh.orig.vertices, lh.orig.faces] = freesurfer_read_surf(file.lh.orig);
			else
				msg = sprintf('file does not exist: %s',file.lh.orig);
				warning(msg);
			end

			file.rh.orig = fullfile(surfPath,'rh.orig');
			if exist(file.rh.orig) == 2,
				[rh.orig.vertices, rh.orig.faces] = freesurfer_read_surf(file.rh.orig);
			else
				msg = sprintf('file does not exist: %s',file.rh.orig);
				warning(msg);
			end

		case {'smoothwm'},

			% load smoothwm surfaces
			file.lh.smoothwm = fullfile(surfPath,'lh.smoothwm');
			if exist(file.lh.smoothwm) == 2,
				[lh.smoothwm.vertices, lh.smoothwm.faces] = freesurfer_read_surf(file.lh.smoothwm);
			else
				msg = sprintf('file does not exist: %s',file.lh.smoothwm);
				warning(msg);
			end

			file.rh.smoothwm = fullfile(surfPath,'rh.smoothwm');
			if exist(file.rh.smoothwm) == 2,
				[rh.smoothwm.vertices, rh.smoothwm.faces] = freesurfer_read_surf(file.rh.smoothwm);
			else
				msg = sprintf('file does not exist: %s',file.rh.smoothwm);
				warning(msg);
			end

		case {'white'},

			% load white surfaces
			file.lh.white = fullfile(surfPath,'lh.white');
			if exist(file.lh.white) == 2,
				[lh.white.vertices, lh.white.faces] = freesurfer_read_surf(file.lh.white);
			else
				msg = sprintf('file does not exist: %s',file.lh.white);
				warning(msg);
			end

			file.rh.white = fullfile(surfPath,'rh.white');
			if exist(file.rh.white) == 2,
				[rh.white.vertices, rh.white.faces] = freesurfer_read_surf(file.rh.white);
			else
				msg = sprintf('file does not exist: %s',file.rh.white);
				warning(msg);
			end

		case {'pial'},

			% load pial surfaces
			file.lh.pial = fullfile(surfPath,'lh.pial');
			if exist(file.lh.pial) == 2,
				[lh.pial.vertices, lh.pial.faces] = freesurfer_read_surf(file.lh.pial);
			else
				msg = sprintf('file does not exist: %s',file.lh.pial);
				warning(msg);
			end

			file.rh.pial = fullfile(surfPath,'rh.pial');
			if exist(file.rh.pial) == 2,
				[rh.pial.vertices, rh.pial.faces] = freesurfer_read_surf(file.rh.pial);
			else
				msg = sprintf('file does not exist: %s',file.rh.pial);
				warning(msg);
			end

		case {'inflated'},

			% load inflated surfaces
			file.lh.inflated = fullfile(surfPath,'lh.inflated');
			if exist(file.lh.inflated) == 2,
				[lh.inflated.vertices, lh.inflated.faces] = freesurfer_read_surf(file.lh.inflated);
			else
				msg = sprintf('file does not exist: %s',file.lh.inflated);
				warning(msg);
			end

			file.rh.inflated = fullfile(surfPath,'rh.inflated');
			if exist(file.rh.inflated) == 2,
				[rh.inflated.vertices, rh.inflated.faces] = freesurfer_read_surf(file.rh.inflated);
			else
				msg = sprintf('file does not exist: %s',file.rh.inflated);
				warning(msg);
			end

		case {'curv'},

			% load curvature
			file.lh.curv = fullfile(surfPath,'lh.curv');
			if exist(file.lh.curv) == 2,
				lh.curv = freesurfer_read_curv(file.lh.curv);
			else
				msg = sprintf('file does not exist: %s',file.lh.curv);
				warning(msg);
			end

			file.rh.curv = fullfile(surfPath,'rh.curv');
			if exist(file.rh.curv) == 2,
				rh.curv = freesurfer_read_curv(file.rh.curv);
			else
				msg = sprintf('file does not exist: %s',file.rh.curv);
				warning(msg);
			end

		case {'thickness'},

			% load thickness
			file.lh.thickness = fullfile(surfPath,'lh.thickness');
			if exist(file.lh.thickness) == 2,
				lh.thickness = freesurfer_read_thickness(file.lh.thickness);
			else
				msg = sprintf('file does not exist: %s',file.lh.thickness);
				warning(msg);
			end

			file.rh.thickness = fullfile(surfPath,'rh.thickness');
			if exist(file.rh.thickness) == 2,
				rh.thickness = freesurfer_read_thickness(file.rh.thickness);
			else
				msg = sprintf('file does not exist: %s',file.rh.thickness);
				warning(msg);
			end

		case {'aparc'},

			% load parcellation
			file.lh.aparc = fullfile(labelPath,'lh.aparc.annot');
			if exist(file.lh.aparc) == 2,
				lh.aparc = freesurfer_read_annotation(file.lh.aparc);
			else
				msg = sprintf('file does not exist: %s',file.lh.aparc);
				warning(msg);
			end

			file.rh.aparc = fullfile(labelPath,'rh.aparc.annot');
			if exist(file.rh.aparc) == 2,
				rh.aparc = freesurfer_read_annotation(file.rh.aparc);
			else
				msg = sprintf('file does not exist: %s',file.rh.aparc);
				warning(msg);
			end

		case {'xfm'},

			% load Talairach transform
			file.xfm = fullfile(xfmPath,'talairach.xfm');
			if exist(file.xfm) == 2,
				subject.TalairachXFM = freesurfer_read_talxfm(file.xfm);
			else
				msg = sprintf('file does not exist: %s',file.xfm);
				warning(msg);
			end
	end
end

subject.file = file;
if exist('lh','var'),
	subject.lh = lh;
end
if exist('rh','var'),
	subject.rh = rh;
end

return
