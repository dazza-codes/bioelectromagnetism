function varargout = mrFlatMesh(varargin)
% MRFLATMESH Application M-file for mrFlatMesh.fig
%    FIG = MRFLATMESH launch mrFlatMesh GUI.
%    MRFLATMESH('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 16-Feb-2001 12:31:20

if (nargin == 0)  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
 
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if (nargout > 0)
		varargout{1} = fig;
	end % endif nargout

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end % end try/catch

end % end function


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

function retVar=browse1_Callback(h,eventData,handles)
    %thisFig=gcf;
    set (handles.status,'String','Waiting for gray matter segmentation file.');
    [fileName,inputGrayPath]=uigetfile('*.gray','Get gray file');
    disp (inputGrayPath);
    set(handles.inputGrayPath,'String',[inputGrayPath,fileName]);    
    set (handles.status,'String','Gray file set');
    return;
    

function retVar=browse2_Callback(h,eventData,handles)
    %thisFig=gcf;
    set (handles.status,'String','Waiting for mrGray mesh file.');
    [fileName,inputMeshPath]=uigetfile('*.mrm','Get mesh file');
    set(handles.inputMeshPath,'String',[inputMeshPath,fileName]);
    set (handles.status,'String','Mesh file set');
    return;

function retVar=browse3_Callback(h,eventData,handles)
    %thisFig=gcf;
    set (handles.status,'String','Setting output path.');
    [fileName,outputMeshPath]=uiputfile('*.mat','Select an output file');
    set(handles.inputSavePath,'String',[outputMeshPath,fileName]);
    set (handles.status,'String','Output path set');
return;

function retVar=go_Callback(h,eventData,handles)
    
    set (handles.status,'String','Running unfold.');
    % Get filenames, coords and unfold size
    meshFileName=get(handles.inputMeshPath,'String');
    grayFileName=get(handles.inputGrayPath,'String');
    flatFileName=get(handles.inputSavePath,'String');
	
    xPos=get(handles.editStartX,'String');
    yPos=get(handles.editStartY,'String');
    zPos=get(handles.editStartZ,'String');
    startPos=[eval(xPos),eval(yPos),eval(zPos)];
	
	sagMM=get(handles.editScaleSag,'String');
	axMM=get(handles.editScaleAx,'String');
	corMM=get(handles.editScaleCor,'String');
	scaleFactor=[eval(sagMM),eval(axMM),eval(corMM)];
    unfoldSize=eval(get(handles.editUnfoldSize,'String'));
 
	showFigures=get(handles.showFiguresCheck,'Value');
	saveExtra=get(handles.saveExtraCheck,'Value');
	truePerimDist=get(handles.perimDistCheck,'Value');
	
	set(handles.status,'UserData','unfoldMeshAlpha v1.0 2001');
	% Go, go ,go !!!!
	statusHandle=handles.status;
	busyHandle=handles.busyBar;
	unfoldMeshFromGUI(meshFileName,grayFileName,flatFileName,startPos,scaleFactor,unfoldSize,statusHandle,busyHandle,showFigures,saveExtra,truePerimDist);
	
return;

function retVal=cancel_Callback(h,eventData,handles)

	close(handles.figure1);
return;


function retVal=help_Callback(h,eventData,handles)
    web1 = 'web([''file:///'' which(''unfoldMeshHelp.html'')],''-browser'')';
    web2 = 'web([''file:///'' which(''unfoldMeshHelp.html'')])';
    
	eval(web1,web2);

return;









