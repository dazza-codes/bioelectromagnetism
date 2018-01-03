function [err, Img, Img_rand, Map] = Randomize_Images (List, Percv, Opt)
%
%   [err, Img, Img_rand, Map] = Randomize_Images (List, Perc, [Opt])
%
%Purpose:
%   Read and introduce randomization to an image
%   
%   
%Input Parameters:
%   List: a set of strings delimiting the images to treat
%        (see example)
%   Perc: a vector of N percentage scrambling (0..100)
%       images are created for each value in Perc
%   
%   Opt: is an optional options structure with the 
%       following optional fields
%       .Check: (0/[1])propmts the user to verify list of images
%       .NewExt: a string used as the new extension to the output images
%              The default is an automatic exetension with the percentage
%              scrambling included. This is only used with Opt.Write = 1
%       .Write: write out the scrambled eggs.
%       .ImType: specify the image type if matlab does not guess it correctly.
%              see help imread for available types
%       .Show: (0/[1]) show the images side by side 
%       .Pause: (0/[1]) pauses between successive images
% 
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   Img : Image matrix as read by imread function
%   Img_rand: the scrambled version of the Img
%   Map : The color map returne by imread 
%   
%      
%Key Terms:
%   
%More Info :
%   imread, imwrite
%
%   example:
%       List = {'1*.TIFF', '2*.TIFF'};
%       Perc = [10 35 55]; 
%       Opt.Write = 0;
%       Opt.Pause = 1;
%    [err, Img, Img_rand, Map] = Randomize_Images (List, Perc, Opt);
%
%   
%
%     Author : Ziad Saad
%     Date : Thu Nov 1 14:07:50 EST 2001
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'Randomize_Images';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
Img=[]; Img_rand=[]; Map = [];
if (nargin == 2),
	Opt.Check = [];
	Opt.NewExt = [];
	Opt.Show = [];
	Opt.ImType = [];
	Opt.Write = [];
	Opt.Pause = [];
end

if (~isfield(Opt,'Pause') | isempty(Opt.Pause)), Opt.Pause = 1; end
if (~isfield(Opt,'Check') | isempty(Opt.Check)), Opt.Check = 1; end
if (~isfield(Opt,'Show') | isempty(Opt.Show)), Opt.Show = 1; end
if (~isfield(Opt,'ImType') | isempty(Opt.ImType)), Opt.ImType = ''; end
if (~isfield(Opt,'Write') | isempty(Opt.Write)), Opt.Write = 0; end

%find the list of images to load
[err, ErrMessage, List] = zglobb (List);

N_files = length(List);
if (~N_files),
	err = 1; errordlg('No files found'); return;
end

fprintf (1,'\n%g files found:\n', N_files);
if (Opt.Check),
	for (i=1:1:N_files),
		fprintf(1,'%g-  %s  -%g\n', i, List(i).name, i);
	end
	chc = input('Do Ya like it ([y]/n): ', 's');
	if (isempty(chc)), chc = 'y'; end
	if (zdeblank(chc) ~= 'y'),
		fprintf(1,'So you don''t, fine.\n');
		return;
	end
end

%main thing
for (i=1:1:N_files),
	%load the image
	if (~isempty(Opt.ImType)),
		[Img, Map] = imread(List(i).name, Opt.ImType);
	else
		[Img, Map] = imread(List(i).name);
	end
	%determine dimensionality
	[Nr, Nc, Nd] = size(Img);
	
	for (ip=1:1:length(Percv)),
		Perc = Percv(ip);
		if (~isfield(Opt, 'NewExt') | isempty(Opt.NewExt)),
			NewExt = sprintf('rnd_%g_', Perc);
		else
			NewExt = Opt.NewExt;
		end
		%find the number of voxels to scramble
		Ns = round (Nr.*Nc.*Perc./100); if (Ns > Nr.*Nc), Ns = Nr.*Nc; elseif (Ns < 1), Ns = 1; end

		%create a random vector of the size of Ns
		rvc = randn(Ns, 1);
		[rvcs, irvcs] = sort(rvc);

		%create a random vector of the size of Nr.*Nc
		rv = randn(Nr.*Nc, 1);
		%sort it
		[rvs,irvs] = sort(rv);
		%Now take Ns pixels to randomize
		iToRand = irvs(1:Ns);
		if (Ns < Nr.*Nc), iNotToRand = irvs(Ns+1:Ns); else iNotToRand = []; end

		%create a blank slate
		Img_rand = zeros(size(Img));
		%copy what is not to be randomized to Img
		for (id=1:1:Nd),
			mt = Img(:,:,id);
			PixToRand = mt(iToRand);
			mt(iToRand) = PixToRand(irvcs);
			Img_rand(:,:,id) = mt;
		end

		%write out results
		[S,Ext] = RemoveExtension(List(i).name, '');
		NewName = sprintf('%s%s%s',NewExt, S,  Ext);
		if (Opt.Write),
			if (exist(NewName) == 2),
				err = 1; errordlg(sprintf('File %s exists, will not overwrite.', NewName));
				return;
			end
			if (~isempty(Opt.ImType)),
				imwrite(Img_rand, Map,  NewName, Opt.ImType);
			else
				imwrite(Img_rand, Map, NewName);
			end
		end
		
		%show results
		if (Opt.Show),
			figure(1); clf
			subplot (121); imshow(Img,Map); drawnow;
			title (List(i).name, 'Interpreter', 'none', 'FontSize',18);
			subplot (122); imshow(Img_rand, Map); drawnow;
			title (NewName, 'Interpreter', 'none', 'FontSize',18');
			
			if (Opt.Pause), 
				fprintf(1,'\nPausing ....(hit Enter to continue)\n'); pause; 
			end
		end
		
	end %loop across percentage
end %loop across files

err = 0;
return;

