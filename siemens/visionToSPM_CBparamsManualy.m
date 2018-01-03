function visionToSPM_CBparamsManualy( figHandle)

if get( findobj( figHandle, 'Tag', 'paramsManualy'), 'Value') == 0
   set( findobj( figHandle, 'Tag', 'frameParams'), 'Visible', 'off');
   
   set( findobj( figHandle, 'Tag', 'paramsDistFactor'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsDistFactorDesc'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsMatrixnY'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsMatrixnX'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsMatrixnXDesc'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsNSlices'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsNSlicesDesc'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsFoVx'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsFoVy'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsFoVDesc'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsDistFactor'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsDistFactorDesc'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsSliceThickness'), 'Visible', 'off');
   set( findobj( figHandle, 'Tag', 'paramsSliceThicknessDesc'), 'Visible', 'off');
else
   set( findobj( figHandle, 'Tag', 'frameParams'), 'Visible', 'on');
   
   set( findobj( figHandle, 'Tag', 'paramsDistFactor'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsDistFactorDesc'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsMatrixnY'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsMatrixnX'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsMatrixnXDesc'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsNSlices'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsNSlicesDesc'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsFoVx'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsFoVy'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsFoVDesc'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsDistFactor'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsDistFactorDesc'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsSliceThickness'), 'Visible', 'on');
   set( findobj( figHandle, 'Tag', 'paramsSliceThicknessDesc'), 'Visible', 'on');
   
end   
