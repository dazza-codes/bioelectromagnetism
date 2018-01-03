function AnalyzeMovie

% AnalyzeMovie - Save analyze images to a kine-loop AVI
% 
% Options are the size of the movie and the views.
% 


% NOTE ON UNIX THERE IS NO OPTION TO COMPRESS, SO THE
% FILES GET REALLY BIG!!!
%
% I'll work on compression for a windows version - maybe.
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:23 $

% Licence:  GNU GPL, no express or implied warranties
% 10/2002, Robert C. Welsh <rcwelsh@umich.edu>
%          University of Michigan, Radiology
% 02/2003, Darren.Weber@flinders.edu.au
%          - attempt to adapt for mri_toolbox, rather than
%            depend on spm
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SCCSid  = '0.5';

global BCH; %- used as a flag to know if we are in batch mode or not.

%-GUI setup
%-----------------------------------------------------------------------

SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Analyze Movie Maker',0);
fprintf('AnalyzeMovie Toolbox 0.5\n');

spm('FigName','Analyze Image Movie Maker',Finter,CmdLine);
% get the name of the rois file.

movieSize = spm_input_ui('Movie size?','+1','m',...
			  ['Small (200x200)|', ...
		           'Medium (300x300)|', ...
		           'Large (400x400)|'],[1 2 3],1);

imageView = spm_input_ui('Image View?','+1','m',...
			  ['Axial|', ...
		           'Coronal|', ...
		           'Sagitall|All'],[1 2 3 4],1);

pixSize = [200 300 400];

movieType= ['a' 'c' 's' 'o'];

% smoothing parameters
nMovies = spm_input('Number of movies','+1','i','1',1,[0,Inf]);

if nMovies < 1
  spm('alert','Exiting as you requested.','AnalyzeMovie',[],0);
  return
end


% change this...
spmImgFiles = {};


for iMovie = 1:nMovies,
    
    % change this...
    spmImgFiles{iMovie}  = spm_get([0,Inf],'*.img',sprintf(['Pick' ...
            ' Image files for movie %d'],iMovie),'./',0);
    if (length(spmImgFiles{iMovie})< 2)
        spm('alert','Exiting as you requested.','AnalyzeMovie',[],0);
        return
    end
    
end

% Now extract the files.


figMovie=figure;
close(figMovie);
figMovie=figure(figMovie);
set(figMovie,'visible','off');
set(figMovie,'DoubleBuffer','on');
set(figMovie,'color',[0 0 0]);
curPos = get(figMovie,'Position');
set(figMovie,'Position',[curPos(1) curPos(2) pixSize(movieSize) pixSize(movieSize)]);
for iMovie = 1:nMovies,
    [movieDir movieName] = fileparts(spmImgFiles{iMovie}(1,:));
    newMovieName = fullfile(movieDir,['kine-loop-' ...
            movieType(imageView:imageView) '.avi']);
    theMovie = avifile(newMovieName);
    
    
    % change this...
    spm_progress_bar('Init',size(spmImgFiles{iMovie},1),sprintf(['Movie #' ...
            ' %d of %d'],iMovie,nMovies),'Extracting data');
    
    
    
    % change spmImgFiles...
    for iFile = 1:size(spmImgFiles{iMovie},1),
        
        spm_progress_bar('Set',iFile);
        
        
        % change this...
        volHdr = spm_vol(spmImgFiles{iMovie}(iFile,:));
        
        
        slices = movieView(volHdr);
        %aVol = spm_read_vols(spm_vol(spmImgFiles{iMovie}(iFile,:)));
        %[nX nY nZ] = size(aVol);
        %midX = floor(nX/2);
        %midY = floor(nY/2);
        %midZ = floor(nZ/2);
        switch imageView
        case 1
            aSlice = slices.aS;
            %aSlice = squeeze(aVol(:,:,midZ));
        case 2
            aSlice = rot90(slices.cS,1);
            %aSlice = squeeze(aVol(:,midY,:));
        case 3
            aSlice = rot90(slices.sS,1);
            %aSlice = squeeze(aVol(midX,:,:));
        case 4 
            aSlice = slices.ortho;
        end
        aSlice = aSlice/max(max(max(aSlice)));
        aSlice = aSlice*256;
        image(aSlice);
        axis image;
        colormap(gray(256));
        curFrame = getframe(gca);
        theMovie = addframe(theMovie,curFrame);
    end
    theMovie=close(theMovie);
end

% Change this...
spm_progress_bar('Clear');
spm_clf(Finter);
spm('FigName','Finished',Finter,CmdLine);
spm('Pointer','Arrow');



close(figMovie);

fprintf('\nFinished  making movie : %s\n',newMovieName);

return






%-----------------------------------------------------------

function results = movieView(vHdr)

% A macro to build the views (based on spm_orthoviews)
%
% Robert C. Welsh <rcwelsh@umich.edu>
% Radiology, University of Michigan
%
% 2002.10.29, Version 0.1

bb = [-78 -112 -50; 78 76 85];
Dims = diff(bb)'+1;

TM0 = [ 1 0 0 -bb(1,1)+1;...
	0 1 0 -bb(1,2)+1;...
	0 0 1 0;...
	0 0 0 1];
CM0 = [ 1 0 0 -bb(1,1)+1;...
	0 0 1 -bb(1,3)+1;...
	0 1 1 0;...
	0 0 0 1];
SM0 = [ 0 -1 0 -bb(1,2)+1;...
	0 0 1 -bb(1,3)+1;...
	1 0 0 0;...
	0 0 0 1];

TD = [Dims(1) Dims(2)];
CD = [Dims(1) Dims(3)];
SD = [Dims(2) Dims(3)];

notSure = [1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 1];

TM = inv(TM0*(notSure\vHdr.mat));
CM = inv(CM0*(notSure\vHdr.mat));
SM = inv(SM0*(notSure\vHdr.mat));


%
% Change this...
axialSlice = spm_slice_vol(vHdr,TM,TD,0);
coronalSlice = spm_slice_vol(vHdr,CM,CD,0);
sagittalSlice = spm_slice_vol(vHdr,SM,SD,0);



results.aS = axialSlice;
results.cS = coronalSlice;
results.sS = sagittalSlice;

bigPix = zeros(Dims(2)+Dims(3),Dims(1)+Dims(2));

bigPix(1:Dims(2),1:Dims(1)) = rot90(axialSlice,1);
bigPix(Dims(2)+1:Dims(2)+Dims(3),1:Dims(1)) = rot90(coronalSlice);
bigPix(Dims(2)+1:Dims(2)+Dims(3),Dims(1)+1:Dims(1)+Dims(2)) = ...
    rot90(sagittalSlice);

results.ortho=bigPix;

return
