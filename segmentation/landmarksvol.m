%landmarksvol

function landmarksvol(argu)
% selects points (landmarks) in volume vol

% mcc -B sglcpp landmarksvol.m -d d:\mymatlab\compres

% Rolando.Grave@hcuge.ch 23-6-2000

% estructura usada para los datos y guardada en el 'UserData' de gcf
% Handles:
% data.hfig ==> current figure
% data.hiomenu ==> Input Output menu principal
% data.hreadfld
% data.hreadslp
% data.hreadbin
% data.hreadana
% data.hxcte
% data.hsliderx
% data.himagex
% data.hycte
% data.hslidery
% data.hzcte
% data.hsliderz
% data.hfilein
% data.hfileout
% data.hsavepoint
% data.hedit
% data.hroiselect
% data.hsaveimage
% data.h
% data.hupdatecenter
% data.hplotsphere
% data.hselectslp 
% ********datos********
% data.nx
% data.x
% data.ny
% data.y
% data.nz
% data.z
% data.vol
% data.max
% data.min
% data.linearindex (used only for slp input data)
% data.filein
% data.fileout
% data.volchanged
% data.center (center of the sphere) only for analize files
% data.radio (radio of the sphere)   only for analize files
% data.header ( header of the analyze file or to save in analyze format)

if ~nargin argu='start'; end
argu=lower(argu);

switch argu,
case 'start',
   % image preferences
   iptsetpref('ImshowAxesVisible','on')
   iptsetpref('ImshowBorder','loose')
   % volume data initialization
   data.vol=repmat(uint16(1),[100,100,100]);
   data=nuevovolumen(data,argu);
   data.linearindex=[];
   data.center=[0,0,0];
   data.radio=-1;
   data.header=[];
   
   % handle of the figure
   data.hfig=figure('Name',[blanks(10),'Selecting points in a 3D volume (e.g. MRI)'],...
      'Units','Normalized','NumberTitle','off',...
      'DefaultUicontrolUnits','Normalized');
   
   % IO-menu
   data.hiomenu=uimenu(data.hfig,'Label','&IO menu ');
   data.hreadmri=uimenu(data.hiomenu,'Label','Input from &FLD file ','Callback','landmarksvol(''readfld'')');
   data.hreadslp=uimenu(data.hiomenu,'Label','Input from S&PI / SPR file ','Callback','landmarksvol(''readslp'')');
   data.hreadbin=uimenu(data.hiomenu,'Label','Input from &Binary file ','Callback','landmarksvol(''readbin'')');
   data.hreadana=uimenu(data.hiomenu,'Label','Input from &Analize file ','Callback','landmarksvol(''readana'')');
   data.hfileout=uimenu(data.hiomenu,'Label','&Output file ','Callback','landmarksvol(''fileout'')');
   
   % handle for the xcte axis, image and slider
   data.hxcte=axes; title(['X = ',num2str(data.x)])
   set(data.hxcte,'Tag','xcte',...
      'Position',[ 0.0232142857142857 0.597619047619048 0.3 0.342857142857143 ],...
      'Xtick',[],'Ytick',[],'Ztick',[],'Box','on','Units','Normalized')
   
   data.hsliderx=uicontrol(data.hfig,'Style','slider','Tag','sliderx',...
      'Position',[ 0.0232142857142857 0.545238095238095 0.3 0.0428571428571429 ],...
      'TooltipString','<-- change  X  value -->',...
      'Callback','landmarksvol(''cambiax'')');

   % handle for the ycte axis and the slider
   data.hycte=axes; title(['Y = ',num2str(data.y)])
   set(data.hycte,'Tag','ycte',...
      'Position',[ 0.355357142857143 0.6 0.3 0.342857142857143 ],...
      'Xtick',[],'Ytick',[],'Ztick',[],'Box','on','Units','Normalized')
   
   data.hslidery=uicontrol(data.hfig,'Style','slider','Tag','slidery',...
      'Position',[ 0.355357142857143 0.542857142857143 0.3 0.0428571428571429 ],...
      'TooltipString','<-- change  Y  value -->',...
      'Callback','landmarksvol(''cambiay'')');

   
   % handle for the zcte axis and slider
   data.hzcte=axes; title(['Z = ',num2str(data.z)])
   set(data.hzcte,'Tag','zcte',...
      'Position',[ 0.683928571428571 0.597619047619048 0.3 0.342857142857143 ],...
      'Xtick',[],'Ytick',[],'Ztick',[],'Box','on','Units','Normalized')
   
   data.hsliderz=uicontrol(data.hfig,'Style','slider','Tag','sliderz',...
      'Position',[ 0.683928571428571 0.538095238095238 0.3 0.0428571428571429 ],...
      'TooltipString','<-- change  Z  value -->',...
      'Callback','landmarksvol(''cambiaz'')');
   
   data=actualiza_sliders(data);
   
   % handle for the current value and edition
   kk3=uicontrol(data.hfig,'Style','Text',...
      'Position',[ 0.0267857142857143 0.454761904761905 0.253571428571429 0.0476190476190476 ],...
      'TooltipString',' Change current pixel value ',...
      'String','Edit Current Value ==>','Fontsize',9','FontWeight','bold');
   
   data.hedit=uicontrol(data.hfig,'Style','edit',...
      'Position',[ 0.280357142857143 0.447619047619048 0.105357142857143 0.0619047619047619 ],...
      'TooltipString',' Change current pixel value ',...
      'String','0','Callback','landmarksvol(''editval'')');
   data.volchanged=0;
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   
   
   % input and output files
   kk1=uicontrol(data.hfig,'Style','Text',...
      'Position',[ 0.674107142857143 0.428571428571429 0.305059523809524 0.0555555555555556 ],...
      'String','Input file (IO menu)','Fontsize',9','FontWeight','bold',...
      'TooltipString',' Select Input file in IO menu ');
   
   data.filein='No file selected';
   data.hfilein=uicontrol(data.hfig,'Style','Text',...
      'Position',[ 0.675595238095238 0.373015873015873 0.305059523809524 0.0555555555555556 ],...
      'String',data.filein,'Fontsize',10','FontWeight','bold',...
      'ForegroundColor','red','TooltipString',' Select Input file in IO menu ');
   
   kk2=uicontrol(data.hfig,'Style','Text',...
      'Position',[ 0.674107142857143 0.285714285714286 0.305059523809524 0.0555555555555556 ],...
      'String','Output file (IO menu)','Fontsize',9','FontWeight','bold',...
      'TooltipString',' Select Output file in IO menu ');
   
   data.fileout='No file selected';
   data.hfileout=uicontrol(data.hfig,'Style','Text',...
      'Position',[ 0.674107142857143 0.23015873015873 0.305059523809524 0.0555555555555556 ],...
      'String',data.fileout,'Fontsize',10','FontWeight','bold',...
      'TooltipString',' Select Output file in IO menu ',...
      'ForegroundColor','red');
   
   data.hsavepoint=uicontrol(data.hfig,'Style','push','Tag','savepoint',...
      'Position',[ 0.68452380952381 0.111111111111111 0.287202380952381 0.0674603174603175 ],...
      'String','Save current point (x,y,z)',...
      'Enable','off',...
      'TooltipString',' Save current X,Y,Z to Output file ',...
      'Callback','landmarksvol(''savepoint'')');
   
   data.hroiselect=uicontrol(data.hfig,'Style','push','Tag','roiselect',...
      'Position',[ 0.686011904761905 0.0214285714285714 0.2875 0.0690476190476191 ],...
      'String',' Select & Save  ROI ',...
      'Enable','on',...
      'TooltipString',' Select ROI with mouse and save in User file',...
      'Callback','landmarksvol(''roiselect'')');
   
   data.hsaveimage=uicontrol(data.hfig,'Style','push','Tag','saveimage',...
      'Position',[ 0.0392857142857143 0.030952380952381 0.2875 0.0690476190476191 ],...
      'String','Save Image',...
      'Enable','off',...
      'TooltipString',' Save image in Analyze format ',...
      'Callback','landmarksvol(''saveimage'')');
   
   data.hchangeimage=uicontrol(data.hfig,'Style','push','Tag','changeimage',...
      'Position',[ 0.0392857142857143 0.10952380952381 0.2875 0.0690476190476191 ],...
      'String','Change Image values',...
      'Enable','on',...
      'TooltipString',' Edit / Change pixel values by group ',...
      'Callback','landmarksvol(''changeimage'')');
   
data.hupdatecenter=uicontrol(data.hfig,'Style','push','Tag','updatecenter',...
      'Position',[ 0.0392857142857143 0.188095238095238 0.2875 0.0690476190476191 ],...
      'String','Update center',...
      'Enable','off',...
      'TooltipString',' Move to the integer closer to the center of the sphere ',...
      'Callback','landmarksvol(''updatecenter'')');

  data.hplotsphere=uicontrol(data.hfig,'Style','push','Tag','plotsphere',...
      'Position',[ 0.0392857142857143 0.266666666666667 0.2875 0.0690476190476191 ],...
      'String','Plot Sphere',...
      'Enable','off',...
      'TooltipString',' Plot best fitting sphere ',...
      'Callback','landmarksvol(''plotsphere'')');

  data.hselectslp=uicontrol(data.hfig,'Style','push','Tag','selectslp',...
      'Position',[ 0.0392857142857143 0.345 0.2875 0.0690476190476191 ],...
      'String','Solution Points Grid',...
      'Enable','off',...
      'TooltipString',' Select solution points in the gray matter  ',...
      'Callback','landmarksvol(''selectslp'')');
   
  pinta_data_puntoxyz(data)
   % para que no cree problemas si el cursor esta en la figura al comenzar
   set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
      'Interruptible','on');      

  % saving data
   set(data.hfig,'UserData',data)
   
case 'fileout'
   data=get(gcf,'Userdata');
   [filename, pathname] = uiputfile([data.filein,'.txt'], 'Save As');
   if filename~=0,
      data.fileout=[pathname,filename];
      if exist(data.fileout)==2,
         delete(data.fileout);
      end
      set(data.hfileout,'String',[pathname(1:3),'...\',filename])
      set(data.hsavepoint,'Enable','on')
      set(data.hfig,'UserData',data);
   end   
   
case 'savepoint'
   data=get(gcf,'Userdata');
   fid=fopen(data.fileout,'at');
   if isempty(data.linearindex),
      fprintf(fid,'%8d  %8d  %8d\n',[data.x,data.y,data.z])
   else % the volume position and the location in the slp file.
      otro=sub2ind([data.nx,data.ny,data.nz],data.x,data.y,data.z);
      pp=find(data.linearindex==otro);
      fprintf(fid,'%8d  %8d %8d %8d\n',[data.x,data.y,data.z,pp])
   end
   fclose(fid);
   
case 'readslp'
   data=get(gcf,'Userdata');
   [filename,pathname]=uigetfile('*.sp?',' File of solution points ');
   if filename~=0,
      data.vol=[];    data.linearindex=[];
      set(data.hfig,'WindowButtonMotionFcn','','Pointer','watch')
      data.filein=[pathname,filename];
      nada=load(data.filein);
      [data.vol,data.linearindex]=slp2vol(nada,[256,160,256]); clear nada %!!!!!
      data=nuevovolumen(data,argu);
      data=actualiza_sliders(data);
      % salvando data y pintando
      pinta_data_puntoxyz(data);
      set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
      data.fileout='No file selected';
      set(data.hfileout,'String',data.fileout)
      set(data.hsaveimage,'Enable','off')
      set(data.hsavepoint,'Enable','off')
      set(data.hfilein,'String',[pathname(1:3),'...\',filename])
      set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
         'Interruptible','on','Pointer','arrow');      
      set(data.hfig,'UserData',data);
   end
   
case 'readfld'
   data=get(gcf,'Userdata');
   [data.vol,pathname,filename]=fld2vol;
   if filename~=0,
      set(data.hfig,'WindowButtonMotionFcn','','Pointer','watch')
      data.linearindex=[];
      data.filein=[pathname,filename];
      data=nuevovolumen(data,argu);
      data=actualiza_sliders(data);
      % salvando data y pintando
      pinta_data_puntoxyz(data)
      set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
      data.fileout='No file selected';
      set(data.hfileout,'String',data.fileout)
      set(data.hsaveimage,'Enable','off')
      set(data.hsavepoint,'Enable','off')
      set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
         'Interruptible','on','Pointer','arrow');      
      set(data.hfilein,'String',[pathname(1:3),'...\',filename])
      set(data.hfig,'UserData',data);
   end
   
case 'readbin'
   data=get(gcf,'Userdata');
   [filename,pathname]=uigetfile('*.*',' Binary file containing the volume. Size(Vol)=[1..Nx,1..Ny,1..Nz ] ');
   if filename~=0,
      data.vol=[]; data.linearindex=[];
      set(data.hfig,'WindowButtonMotionFcn','','Pointer','watch')
      
      data.filein=[pathname,filename];
      prompt={'Number of slices in x (Nx):','Number of slices in y (Ny):',...
            'Number of slices in z (Nz):', 'File type :'};
      def={'0','0','0','uint8 or uint16 or uint32'};
      dlgTitle='*** Volume Dimensions [Nx,Ny,Nz]****';
      lineNo=1;
      answer=inputdlg(prompt,dlgTitle,lineNo,def);
      data.vol=Read_Binary_Vol4(data.filein,str2num(char(answer(1))),str2num(char(answer(2))),str2num(char(answer(3))),1,deblank(char(answer(4))));
      data=nuevovolumen(data,argu);
      data=actualiza_sliders(data);
      % salvando data y pintando
      pinta_data_puntoxyz(data);
      set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
      data.fileout='No file selected';
      set(data.hfileout,'String',data.fileout)
      set(data.hsaveimage,'Enable','off')
      set(data.hsavepoint,'Enable','off')
      set(data.hfilein,'String',[pathname(1:3),'...\',filename])
      set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
         'Interruptible','on','Pointer','arrow');      
      set(data.hfig,'UserData',data);
   end
   
case 'readana'
   data=get(gcf,'Userdata');
   [filename,pathname]=uigetfile('*.img',' Analize image file ');
   if filename~=0,
      data.vol=[];   data.linearindex=[];
      set(data.hfig,'WindowButtonMotionFcn','','Pointer','watch')
      data.filein=[pathname,filename]; kk=length(data.filein);
      [data.vol,data.header]=Read_Ana_File(data.filein(1:kk-4));
                                                 
      data=nuevovolumen(data,argu);
      data=actualiza_sliders(data);
      status=ana_header_io(data.header,-6);
      if status==12345,
          temp=ana_header_io(data.header,-5);
          data.center=temp(1:3);
          data.radio=temp(4);
          set(data.hplotsphere,'Enable','on');
          set(data.hupdatecenter,'Enable','on');
          %set(data.hselectslp,'Enable','on');
      end
      set(data.hselectslp,'Enable','on');
      % salvando data y pintando
      pinta_data_puntoxyz(data);
      set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
      set(data.hfilein,'String',[pathname(1:3),'...\',filename])
      data.fileout='No file selected';
      set(data.hfileout,'String',data.fileout)
      set(data.hsaveimage,'Enable','off')
      set(data.hsavepoint,'Enable','off')
      set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
         'Interruptible','on','Pointer','arrow');      
      set(data.hfig,'UserData',data);
   end
   
case 'movimiento', % flecha si esta fuera de las 3 figuras
   data=get(gcf,'Userdata');
   [Inx,Iny,Inz]=DondeEsta(data);
   
   if Inx,
      set(data.hfig,'Pointer','fullcrosshair',...
         'WindowButtonUpFcn','landmarksvol(''mouseselectx'')','Interruptible','on')
   elseif Iny
      set(data.hfig,'Pointer','fullcrosshair',...
         'WindowButtonUpFcn','landmarksvol(''mouseselecty'')','Interruptible','off')
   elseif Inz
       set(data.hfig,'Pointer','fullcrosshair',...
         'WindowButtonUpFcn','landmarksvol(''mouseselectz'')','Interruptible','off')
   else
      set(data.hfig,'Pointer','arrow','WindowButtonUpFcn','')
   end   
   
case 'mouseselectx',
   data=get(gcf,'Userdata');
   set(data.hfig,'Pointer','watch');
   cp=get(data.hxcte,'CurrentPoint');
   
   data.z=round(cp(2,1)); data=controla_mouse_point(data,3);
   data.y=round(cp(2,2)); data=controla_mouse_point(data,2);
   
   % salvando data y pintando
   data=new_sliders_pos(data);
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   pinta_data_puntoxyz(data)
   set(data.hfig,'Pointer','arrow');
   set(data.hfig,'UserData',data);
   
case 'mouseselecty',
   data=get(gcf,'Userdata');
   set(data.hfig,'Pointer','watch');
   cp=get(data.hycte,'CurrentPoint');
   
   data.z=round(cp(2,1)); data=controla_mouse_point(data,3);
   data.x=round(cp(2,2)); data=controla_mouse_point(data,1);
   
   % salvando data y pintando
   data=new_sliders_pos(data);
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   pinta_data_puntoxyz(data)
   set(data.hfig,'Pointer','arrow');
   set(data.hfig,'UserData',data);
   
case 'mouseselectz',
   data=get(gcf,'Userdata');
   set(data.hfig,'Pointer','watch');
   cp=get(data.hzcte,'CurrentPoint');
   
   data.y=round(cp(2,1)); data=controla_mouse_point(data,2);
   data.x=round(cp(2,2)); data=controla_mouse_point(data,1);
   
   % salvando data y pintando
   data=new_sliders_pos(data);
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   pinta_data_puntoxyz(data)
   set(data.hfig,'Pointer','arrow');
   set(data.hfig,'UserData',data);
   
case 'cambiax',
   data=get(gcf,'Userdata');
   set(data.hfig,'Pointer','watch');
   data.x=round(get(data.hsliderx,'Value'));
   set(data.hsliderx,'Value',data.x)
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   % salvando data y pintando
   pinta_data_puntoxyz(data)
   set(data.hfig,'Pointer','arrow');
   set(data.hfig,'UserData',data);
   
case 'cambiay',
   data=get(gcf,'Userdata');
   set(data.hfig,'Pointer','watch');
   data.y=round(get(data.hslidery,'Value'));
   set(data.hslidery,'Value',data.y)
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
 % salvando data y pintando
   pinta_data_puntoxyz(data)
   set(data.hfig,'Pointer','arrow');
   set(data.hfig,'UserData',data);

case 'cambiaz',
   data=get(gcf,'Userdata');
   set(data.hfig,'Pointer','watch');
   data.z=round(get(data.hsliderz,'Value'));
   set(data.hsliderz,'Value',data.z)
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   % salvando data y pintando
   pinta_data_puntoxyz(data)
   set(data.hfig,'Pointer','arrow');
   set(data.hfig,'UserData',data);
   
case 'editval',
   data=get(gcf,'Userdata');
   valor=str2num(get(data.hedit,'String')),
   if ~isempty(valor) & ~isnan(valor) & length(valor)==1,
      data.vol(data.x,data.y,data.z)=valor;
      data.volchanged=1;
      set(data.hsaveimage,'Enable','on')
      data.min=min(data.vol(:));
      data.max=max(data.vol(:));
   else
     set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   end
  set(data.hfig,'UserData',data);
  
case 'saveimage'  %Salvar usando Ana_file
   data=get(gcf,'Userdata');
   set(data.hfig,'WindowButtonMotionFcn','');
   kk=length(data.filein);
   [filename,pathname]=uiputfile(data.filein, 'Save Transformed Image as :');
   if filename~=0,
      data.filein=[pathname,filename];  kk=length(data.filein);
      set(data.hfilein,'String',[pathname(1:3),'...\',filename])
      % analize format in action
      if isempty(data.header), data.header=IsAnalize(data); end
      Save_Ana_file(data.vol,data.header,data.filein(1:kk-4));
      % end of analize stuff
      set(data.hsaveimage,'Enable','off')
      data.volchanged=0;
   else
      myerror(' Image cannot be stored / filename');
   end
   set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
      'Interruptible','on');      
   set(data.hfig,'UserData',data); 
   
case 'changeimage'
    data=get(gcf,'Userdata');
    set(data.hfig,'WindowButtonMotionFcn','');
    dlgTitle=' Editing pixel values ';
    prompt={' Change all pixels with value : ',' To the new pixel value :'};
    def={num2str(double(data.vol(data.x,data.y,data.z))),'0'};
    lineNo=1;
    dimelo=inputdlg(prompt,dlgTitle,lineNo,def); 
    if ~isempty(dimelo),
        orival=str2num(deblank(char(dimelo(1))));
        kk=find(data.vol(:)==orival);
        newval=str2num(deblank(char(dimelo(2))));
        temp=data.vol(data.x,data.y,data.z);
        if (orival>=data.min) & (orival<=data.max) & ~isempty(kk)...
                & (newval>=0) & (newval<=65535) & abs(orival-newval)~=0,
            data.vol(kk)=newval;
            data.volchanged=1;
            data.min=min([data.min,newval]);
            data.max=max([data.max,newval]);
            if temp==orival, % update the "edit current value" value
                set(data.hedit,'String',char(dimelo(2)));
            end
            set(data.hsaveimage,'Enable','on')
        end
    end
    set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
        'Interruptible','on');      
    set(data.hfig,'UserData',data); 
    
case 'roiselect'
   data=get(gcf,'Userdata');
   set(data.hroiselect,'Enable','off');
   set(data.hfig,'WindowButtonMotionFcn','');
   kk=waitforbuttonpress;
   [Inx,Iny,Inz]=DondeEsta(data);
   if Inx,
      [z,y]=myroipoly(data.hxcte); [z,y]
      y(y<1)=1; y(y>data.ny)=data.ny;
      z(z<1)=1; z(z>data.nz)=data.nz;
      set(data.hroiselect,'Enable','on');
      set(data.hfig,'UserData',data);
      BW=roipoly(squeeze(data.vol(data.x,:,:)),z,y);
      [zz,yy]=find(BW); clear BW;
      xx=ones(size(yy))*data.x;
   elseif Iny,
      [z,x]=myroipoly(data.hycte);
      x(x<1)=1; x(x>data.nx)=data.nx;
      z(z<1)=1; z(z>data.ny)=data.nz;
      set(data.hroiselect,'Enable','on');
      set(data.hfig,'UserData',data);
      BW=roipoly(squeeze(data.vol(:,data.y,:)),z,x);
      [zz,xx]=find(BW); clear BW
      yy=ones(size(xx))*data.y;
   elseif Inz,
      [y,x]=myroipoly(data.hzcte);
      x(x<1)=1; x(x>data.nx)=data.nx;
      y(y<1)=1; y(y>data.ny)=data.ny;
      set(data.hroiselect,'Enable','on');
      set(data.hfig,'UserData',data);
      Idata=getimage(data.hzcte);clear Idata
      BW=roipoly(squeeze(data.vol(:,:,data.z)),y,x);
      [xx,yy]=find(BW); clear BW
      zz=ones(size(xx))*data.z;
   end
   [filename,pathname]=uiputfile(data.fileout, 'Save ROI in file ( If file already exists just APPENDS!!)');
   if filename~=0,
      data.fileout=[pathname,filename];
      set(data.hfileout,'String',[pathname(1:3),'...\',filename])
      set(data.hsavepoint,'Enable','on')
      fid=fopen(data.fileout,'at');
      if isempty(data.linearindex),
         temp=[xx(:),yy(:),zz(:)];
         for kk=1:length(xx), fprintf(fid,'%8d  %8d  %8d\n',temp(kk,:)); end
      else % the volume position and the location in the slp file.
         otro=sub2ind([data.nx,data.ny,data.nz],xx,yy,zz);
         for kk=1:length(xx),
            pp=find(data.linearindex==otro(kk));
            if isempty(pp), pp=-1111; end
            temp=[xx(kk),yy(kk),zz(kk),pp];
            fprintf(fid,'%8d  %8d %8d %8d\n',temp); 
         end
      end
      fclose(fid);
   else
      myerror(' ROI cannot be stored / filename');
   end
   set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
      'Interruptible','on');      
   set(data.hfig,'UserData',data); 
   
case 'plotsphere'
    data=get(gcf,'Userdata');
    t=linspace(0,2*pi,100); 
    % in x
    radio=sqrt( data.radio^2 - (data.x-data.center(1) )^2 );
    y=radio*sin(t)+data.center(2); z=radio*cos(t)+data.center(3);
    axes(data.hxcte),
    hold on, plot(z,y,'.-'), plot(data.center(3),data.center(2),'*'), hold off
    
    % in y
    radio=sqrt( data.radio^2 - (data.y-data.center(2) )^2 );
    y=radio*sin(t)+data.center(1); z=radio*cos(t)+data.center(3);
    axes(data.hycte),
    hold on, plot(z,y,'.-'), plot(data.center(3),data.center(1),'*'),  hold off

    % in z
    radio=sqrt( data.radio^2 - (data.y-data.center(3) )^2 );
    y=radio*sin(t)+data.center(1); z=radio*cos(t)+data.center(2);
    axes(data.hzcte),
    hold on, plot(z,y,'.-'), plot(data.center(2),data.center(1),'*'), hold off

    
    
case 'updatecenter'
   data=get(gcf,'Userdata');
   data.x=round(data.center(1)); data.y=round(data.center(2)); data.z=round(data.center(3)); 
   set(data.hsliderx,'Value',data.x)
   set(data.hslidery,'Value',data.y)
   set(data.hsliderz,'Value',data.z)
   set(data.hedit,'String',num2str(double(data.vol(data.x,data.y,data.z))));
   pinta_data_puntoxyz(data)
   set(data.hfig,'UserData',data); 
   
    
case 'selectslp',
    data=get(gcf,'UserData');
    set(data.hfig,'WindowButtonMotionFcn','');
    deltaint=ana_header_io(data.header,-2);
    Np=sum(data.vol(:)==13); 
    mystring=['Points Interdistance (slices)', blanks(2),'Dx =',num2str(deltaint(1)),blanks(2),'Dy =',num2str(deltaint(2)),blanks(2),'Dz =',num2str(deltaint(3))];
    myhelp=helpdlg({['Number of solution points =',num2str(Np)],mystring},'Current Grid of solution');
    
    dlgTitle=' (Re-) Computation of the Grid of Solution Points ';
    prompt={[' Interdistance  in X direction (current value =', num2str(deltaint(1)),' )' ],...
            [' Interdistance  in Y direction (current value =', num2str(deltaint(2)),' )' ],...
            [' Interdistance  in Z direction (current value =', num2str(deltaint(3)),' )' ],...
            [' Defining isolated points. Number of neighbors smaller than:']};
    def={num2str(deltaint(1)),num2str(deltaint(2)),num2str(deltaint(3)),'1 to 6'};
    lineNo=1;
    otravez=1; 
    while otravez,
        estabien=0;
        dimelo=inputdlg(prompt,dlgTitle,lineNo,def); 
        %delta=deltaint;
        if ~isempty(dimelo),
            delta(1)=round(str2num(deblank(char(dimelo(1)))));
            delta(2)=round(str2num(deblank(char(dimelo(2)))));
            delta(3)=round(str2num(deblank(char(dimelo(3)))));
            Nmin=round(str2num(deblank(char(dimelo(4)))));
            estabien=isinrang(delta(1),[1,round(data.nx/2)]) &...
                     isinrang(delta(2),[1,round(data.ny/2)]) &...
                     isinrang(delta(3),[1,round(data.nz/2)])  &...
                     isinrang(Nmin,[1,6]);
        end
        otravez=(~isempty(dimelo) ) &  (~estabien );
    end
    if ~isempty(dimelo),
        deltaint=delta,
        VolOut=select_slp_in_vol(data.vol,deltaint(1),deltaint(2),deltaint(3),Nmin);
        Np=sum(VolOut(:)==13); 
        button = questdlg(['There are ',num2str(Np),' solution points. Do you accept this grid?' ],'New Grid Computation Results','Yes','No','No');
        if strcmp(button,'Yes')
            data.vol=VolOut;
            data.header=ana_header_io(data.header,2,deltaint);
            data.volchanged=1;
            set(data.hsaveimage,'Enable','on')
            data.min=min(data.vol(:));
            data.max=max(data.vol(:));
            set(data.hfig,'UserData',data); 
        end
    end
    close(myhelp)
    set(data.hfig,'WindowButtonMotionFcn','landmarksvol(''movimiento'')',...
        'Interruptible','on');      
    
otherwise
    myerror('Option not available')
end

function pinta_data_puntox(data)
axes(data.hxcte),
h=imshow(squeeze(data.vol(data.x,:,:)),[data.min, data.max]); hold on
%colormap(jet(256))
%h=imagesc(squeeze(data.vol(data.x,:,:)),[data.min, data.max]);
%set(h,'EraseMode','background');
scatter(data.z,data.y,'+','r','filled'), hold off
set(data.hxcte,'Xtick',[],'Ytick',[],'Ztick',[],'Box','on')   
title(['X = ',num2str(data.x)])

function pinta_data_puntoy(data)
axes(data.hycte),
h=imshow(squeeze(data.vol(:,data.y,:)),[data.min, data.max]); hold on
%colormap(jet(256))
%h=imagesc(squeeze(data.vol(:,data.y,:)),[data.min, data.max]);
%set(h,'EraseMode','background');
scatter(data.z,data.x,'+','r','filled'), hold off
set(data.hycte,'Xtick',[],'Ytick',[],'Ztick',[],'Box','on')   
title(['Y = ',num2str(data.y)])

function pinta_data_puntoz(data)
axes(data.hzcte)
h=imshow(squeeze(data.vol(:,:,data.z)),[data.min, data.max]); hold on
%colormap(jet(256))
%h=imagesc(squeeze(data.vol(:,:,data.z)),[data.min, data.max]);
%set(h,'EraseMode','background');
scatter(data.y,data.x,'+','r','filled'), hold off
set(data.hzcte,'Xtick',[],'Ytick',[],'Ztick',[],'Box','on')   
title(['Z = ',num2str(data.z)])

function pinta_data_puntoxyz(data)
pinta_data_puntox(data)
pinta_data_puntoy(data)
pinta_data_puntoz(data)
drawnow

function data=nuevovolumen(data,argu)
[data.nx,data.ny,data.nz]=size(data.vol);
data.x=round(data.nx*0.5);
data.y=round(data.ny*0.5);
data.z=round(data.nz*0.5);
data.center=[data.x,data.y,data.z];
data.max=max(data.vol(:));
data.min=min(data.vol(:));
data.fileout='No file selected';
data.volchanged=0;


function data=actualiza_sliders(data)
set(data.hsliderx,'Min',1,'Max',data.nx,'SliderStep',[1/data.nx, 0.25],'Value',data.x)
set(data.hslidery,'Min',1,'Max',data.ny,'SliderStep',[1/data.ny, 0.25],'Value',data.y)
set(data.hsliderz,'Min',1,'Max',data.nz,'SliderStep',[1/data.nz, 0.25],'Value',data.z)

function data=new_sliders_pos(data)
set(data.hsliderx,'Value',data.x)
set(data.hslidery,'Value',data.y)
set(data.hsliderz,'Value',data.z)

function data=controla_mouse_point(data, cual)
switch cual,
case 1,
   if data.x<1, data.x=1; end
   if data.x>data.nx, data.x=data.nx; end
case 2,
   if data.y<1, data.y=1; end
   if data.y>data.ny, data.y=data.ny; end
case 3,
   if data.z<1, data.z=1; end
   if data.z>data.nz, data.z=data.nz; end
end

function [xi,yi]=myroipoly(axishandle)
[xi,yi] = getline(axishandle,'closed');

function [Inx,Iny,Inz]=DondeEsta(data);
% posicion de los ejes
posx=get(data.hxcte,'Position');
posy=get(data.hycte,'Position');
posz=get(data.hzcte,'Position');
% posicion del cursor
cp=get(data.hfig,'CurrentPoint');
x = cp(1,1);
y = cp(1,2);
Inx=(x>posx(1) & x<posx(1)+posx(3) & y>posx(2) & y<posx(2)+posx(4));
Iny=(x>posy(1) & x<posy(1)+posy(3) & y>posy(2) & y<posy(2)+posy(4));
Inz=(x>posz(1) & x<posz(1)+posz(3) & y>posz(2) & y<posz(2)+posz(4));
% over

function header=IsAnalize(data)
% creates a header to save in Analize format 

%[PATH,NAME,EXT,VER] = FILEPARTS(data.filein);
%filename=[PATH,'\',NAME];
%if exist([PATH,'\',NAME,'.img'])==2 & exist([PATH,'\',NAME,'.hdr'])==2, % Analize
   %[header,formato]=Read_Ana_Header([PATH,'\',NAME,'.hdr']);
   %else
   header=Init_Ana_Header;
   header.dime.dim(1)=4;  % number of dimensions
   header.dime.dim(2)=data.nx;    % pixels in X
   header.dime.dim(3)=data.ny;    % pixels in Y
   header.dime.dim(4)=data.nz;    % pixels in Z slices
   header.dime.dim(5)=1;    % number of time frames
   if data.max<=65535 & data.max>=256, 
      header.dime.bitpix=16;  % use uint16
   elseif data.max<=255, 
      header.dime.bitpix=8;    % use uint16
   else
      myerror( ' Unable to detect number of bit per pixel ');
   end
   %end

   
   function local_selectslp(pregunta)
   %scrsz = get(0,'ScreenSize');
   % the handle of the second figure 
   fh=figure('MenuBar','None','Units','Normalized','NumberTitle','off',...
       'DefaultUicontrolUnits','Normalized','Name',[blanks(20),'Volume & Grid Information'],...
       'Position',[0.4600    0.0534    0.3555    0.5391]);

   % copy the data to the new figure data
   set(fh,'UserData',data);
   
   % the Volume Information Part       
   %hvolinf=uicontrol(fh,'Style','frame','Position',[-0.0027472527472527475, 0.0024154589371980675,1,1]);
   hvolTitle=uicontrol(fh,'Style','text','String','Volume Information',...
       'Position',[-0.0027472527472527475,0.9323671497584539,1,0.06521739130434782],...
       'ForeGroundColor','Red');
   
   mystring=['Number of slices', blanks(5),'Nx =',num2str(data.nx),blanks(5),'Ny =',num2str(data.ny),blanks(5),'Nz =',num2str(data.nz)];
   hvolDim=uicontrol(fh,'Style','text','String',mystring,...
       'Position',[-0.0027472527472527475,0.8743961352657005,1,0.06521739130434782],...
       'ForeGroundColor','White','Fontsize',9,'FontWeight','bold');

   mystring=['Pixel sizes ',blanks(5), 'Sx =',num2str(data.header.dime.pixdim(2)),blanks(5),...
             'Sy =',num2str(data.header.dime.pixdim(3)),blanks(5),...
             'Sz =',num2str(data.header.dime.pixdim(4)),blanks(5), 'in ',data.header.dime.vox_units];
   hvolDim=uicontrol(fh,'Style','text','String',mystring,...
       'Position',[-0.0027472527472527475,0.8115942028985507,1,0.06521739130434782],...
       'ForeGroundColor','White','Fontsize',9,'FontWeight','bold');
   
   
   % the Grid Information Part       
   %hgridinf=uicontrol(fh,'Style','frame');
   hgridTitle=uicontrol(fh,'Style','text','String','Grid Information',...
       'Position',[-0.0027472527472527475,0.6666666666666666,1,0.06521739130434782],...
       'ForeGroundColor','Red');

   mystring=['Points Interdistance ', blanks(5),'Dx =',num2str(data.nx),blanks(5),'Dy =',num2str(data.ny),blanks(5),'Dz =',num2str(data.nz)];
   hgridDim=uicontrol(fh,'Style','text','String',mystring,...
       'Position',[-0.0027472527472527475,0.6038647342995169,1,0.06521739130434782],...
       'ForeGroundColor','White','Fontsize',9,'FontWeight','bold');
 
   mystring=[blanks(5), 'Number of solution points = ',num2str(sum(data.vol(:)==13)) ];
   hgridNslp=uicontrol(fh,'Style','text','String',mystring,...
       'Position',[-0.0027472527472527475,0.538647342995169,1,0.06521739130434782],...
       'ForeGroundColor','White','Fontsize',9,'FontWeight','bold');
   
   mystring=[blanks(5), 'Press here to (re) compute the Grid of solution points' ];
   hgridGrid=uicontrol(fh,'Style','pushbutton','String',mystring,...
       'Position',[-0.0027472527472527475,0.4492753623188406,1,0.06521739130434782],...
       'ForeGroundColor','White','Fontsize',9,'FontWeight','bold', 'callback','landmarksvol(''computegrid'')');
   uiwait(fh)
   
   
function local_compute_grid
data=get(gcf,'UserData');
deltaint=ana_header_io(data.header,-2);
dlgTitle=' (Re-) Computation of the Grid of Solution Points ';
prompt={[' Interdistance  in X direction (current value =', num2str(deltaint(1)),' )' ],...
        [' Interdistance  in Y direction (current value =', num2str(deltaint(2)),' )' ],...
        [' Interdistance  in Z direction (current value =', num2str(deltaint(3)),' )' ]};
def={num2str(deltaint(1)),num2str(deltaint(2)),num2str(deltaint(3))};
lineNo=1;
otravez=1; 
while otravez,
    estabien=0;
    dimelo=inputdlg(prompt,dlgTitle,lineNo,def); 
    %delta=deltaint;
    if ~isempty(dimelo),
        delta(1)=round(str2num(deblank(char(dimelo(1)))));
        delta(2)=round(str2num(deblank(char(dimelo(2)))));
        delta(3)=round(str2num(deblank(char(dimelo(3)))));
        estabien=isinrang(delta(1),[1,round(data.nx/2)]) & isinrang(delta(2),[1,round(data.ny/2)]) & isinrang(delta(3),[1,round(data.nz/2)]);
    end
    otravez=(~isempty(dimelo) ) &  (~estabien );
end
if ~isempty(dimelo),
    deltaint=delta;
    VolOut=select_slp_in_vol(data.vol,deltaint(1),deltaint(2),deltaint(3));
    data.vol=VolOut;
    data.volchanged=1;
    % copy to data to Main figure
    set(data.hfig,'UserData',data); 
end   
   