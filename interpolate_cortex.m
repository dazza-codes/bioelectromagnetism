function [InterpMatrix] = interpolate_cortex(sparseFV,denseFV,vertConn)
% INTERPOLATE_CORTEX Interpolate current density from a sparse tesselation to a dense tesselation
% function InterpMatrix = interpolate_cortex(sparseFV,denseFV,vertConn);
% sparseFV is the sparse Faces/Vertices structure
% denseFV is the dense Faces/Vertices structure
% vertConn is the vertices connectivity of the dense Faces/Vertices structure
%
% InterpMatrix returned is a sparse matrix (nDenseVertices x nSparseVertices) used for interpolation,
% where nDenseVertices in the number of vertices in the dense structure and nSparseVertices is the
% number of vertices in the sparse structure
% Interpolated current density is given by DenseCurrentDensity = InterpMatrix * SparseCurrentDensity
%
% Remarks: Every vertex in sparseFV should have a corresponding one in denseFV (i.e. same coordinates)
% sparseFV=reducepatch(denseFV) matlab function produces such tesselations
% Surfaces should be connected, i.e. there should exist a path through triangles connecting every pair
% of vertices on the surface
% interpolate_cortex performs interpolation in terms of surface neighbours and not in terms of volume
% neighbours
%
% See also REDUCEPATCH VERTICES_CONNECTIVITY PATCH_SWELL

% Dimitrios Pantazis, Ph.D.
% 11-Apr-02
% <copyright>
% <copyright>


%choose whether to display bars
if(~exist('VERBOSE','var')),
   VERBOSE = 1; % default non-silent running of waitbars
end

%Initialize variables
nDenseVertices=size(denseFV.vertices,1);
nSparseVertices=size(sparseFV.vertices,1);
existInSparse=zeros(nDenseVertices,1);
InterpMatrix=sparse(nDenseVertices,nSparseVertices);


%find vertices that have been preserved in sparse patch
%iDense: index of common vertices in denseFV
%iSparse: index of common vertices in sparseFV
[ignore,iDense,iSparse]=intersect(denseFV.vertices,sparseFV.vertices,'rows');
existInSparse(iDense) =iSparse; %helpful vector 

%For a strange reason it takes awful time if I write InterpMatrix(iDense,iSparse)=1;
for i=1:length(iDense)
    InterpMatrix(iDense(i),iSparse(i))=1;
end

%interpolate the rest of the vertices
if(VERBOSE)
    hwait = waitbar(0,sprintf('Interpolating the %.0f remaining vertices...',nDenseVertices-nSparseVertices));
    drawnow %flush the display
    step=round(nDenseVertices/20);
end
for i=1:nDenseVertices
    if(VERBOSE)
       if(~rem(i,step)) % twenty updates
          waitbar(i/nDenseVertices,hwait);
          drawnow %flush the display         
       end
    end
    if(~existInSparse(i))
        startPatch=i;
        nLayer=0;
        neighboursKnown=[];
        neighboursKnownScale=[];
        while(length(neighboursKnown)<3 & nLayer <10 )
            swelledPatch=patch_swell(startPatch,vertConn); %swell patch to find know vertices
            if(isempty(swelledPatch))
                %disp('surface not closed');
                break;
            end         
            temp=existInSparse(swelledPatch);
            neighboursNew=temp(temp>0)';
            nLayer=nLayer+1;
            neighboursKnown=[neighboursKnown neighboursNew];
            neighboursKnownScale=[neighboursKnownScale nLayer*ones(1,length(neighboursNew))];
            startPatch=[swelledPatch startPatch];
        end
        %update InterpMatrix
        InterpMatrix(i,neighboursKnown)=(1./neighboursKnownScale)/sum(1./neighboursKnownScale);
    end
end
if(VERBOSE)
    close(hwait);
end