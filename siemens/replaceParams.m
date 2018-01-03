function params = replaceParams( figHandle, params);

FoV            = get( findobj( figHandle, 'Tag', 'paramsFoVx'), 'String');
if ~isempty(FoV)
   params.FoV(1) = str2num( FoV);   
end
FoV            = get( findobj( figHandle, 'Tag', 'paramsFoVy'), 'String');
if ~isempty(FoV)
   params.FoV(2) = str2num( FoV);   
end


MatrixnX       = get( findobj( figHandle, 'Tag', 'paramsMatrixnX'), 'String');
if ~isempty( MatrixnX)
   params.matrix(1) = str2num( MatrixnX);
end

MatrixnY       = get( findobj( figHandle, 'Tag', 'paramsMatrixnY'), 'String');
if ~isempty( MatrixnY)
   params.matrix(2) = str2num( MatrixnY);
end

NSlices        = get( findobj( figHandle, 'Tag', 'paramsNSlices'), 'String');
if ~isempty( NSlices)
   params.nSlices = str2num( NSlices);
end

SliceThickness = get( findobj( figHandle, 'Tag', 'paramsSliceThickness'), 'String');
if ~isempty( SliceThickness)
   params.sliceThickness = str2num( SliceThickness);
end

DistFactor     = get( findobj( figHandle, 'Tag', 'paramsDistFactor'), 'String');
if ~isempty( DistFactor)
 params.distFactor = str2num( DistFactor);   
end
