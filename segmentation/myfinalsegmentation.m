function nvol = myfinalsegmentation


% Implicit Values

BACKGROUND  =  0;  % not used elsewhere below
ENDGRAY     = 20;  % not used elsewhere below


GRAY        = 11;  % Gray Matter will go from 11 to 20;
GRAYBORDER  = 12;
WHITE       = 21;
SKIN        = 31;
BGPARAMETER = 70;


% BET parameter. Smaller values lead to more brain and strange things
BRAINSTRIP  = .6;
BETDIR      = 'C:\Program Files\MRIcro '; % Implicit directory to localize the bet.exe



% check for existence of bet.exe, available in the MRICRO distribution

% should check the operating platform and then check for either bet.exe
% (windows) or bet (linux).

if exist([BETDIR(1:end-1),'.exe'])~=2,
  [betFilename, betPath] = uigetfile('bet.hdr', 'Pick location of the bet.exe file');
  if isequal(betFilename,0) | isequal(betPath,0)
    disp('User pressed cancel')
  else
    disp(['User selected ', fullfile(betPath, betFilename)])
  end
  BETDIR = fullfile(betPath, betFilename(1:end-4))
  if exist([BETDIR,'.exe'])~=2,
    errStr = ['Segmentation impossible without the bet.exe.\n',...
        'Go to web site www.mricro.com, and download/install mricro.'];
    error(errStr);
  end;
end;


% Select the MRI volume to process
[filename, pathname] = uigetfile('*.hdr', 'Pick the original hdr-file');
if isequal(filename,0) | isequal(pathname,0)
  disp('User pressed cancel')
else
  disp(['User selected ', fullfile(pathname, filename)])
end

ofilename = fullfile(pathname, filename); 
bfilename = fullfile(pathname, ['segmented']); 

% replacing this with avw_read function; DLW 10/2003
%[Vols,head]=Read_Ana_File(ofilename(1:end-4));
%sizeh=head.dime.dim(2:4);
avw = avw_read(ofilename);
sizeh = avw.hdr.dime.dim(2:4);


%SKULL STRIPPING BY MEANS OF BET. It will create a file on the same
%directory called segmented that will be deleted afterwards

skullStripped = 0;
while ~skullStripped,
  
  h = helpdlg('BET based brain segmentation in progress','Skull stripping using BET');
  set(h,'Pointer','watch');
  t = [BETDIR,pathname,filesep,filename,' ',pathname,filesep,'segmented',' -s -f ',num2str(BRAINSTRIP)];
  status = dos([t]),
  set(h,'Pointer','arrow');
  close(h);
  
  %[Volb,head]=Read_Ana_File(bfilename);
  %avw.img=Volb; avw.hdr=head;
  avw = avw_read(bfilename);
  hh=avw_view(avw) 
  uiwait(hh)
  
  ans = questdlg('Do you like the Brain Segmentation Result?','Brain Segmentation Quality')
  switch ans,
    case 'No', 
      skullStripped = 0;
      prompt = {['Enter new brain strip parameter. Increasing it gives less brain and surroundings']};
      answer = inputdlg(prompt,'Brain Strip Parameter',1,num2str(BRAINSTRIP));
      BRAINSTRIP = answer{1},
    case 'Cancel',
      error('You do not like the program. What a pity')
    case 'Yes',
      skullStripped = 1;
  end % switch
end


% SELECTING GRAY AND WHITE MATTER FROM SKULL STRIPPED IMAGE

nVol = avw.img;
nVol(:) = 0;

h = helpdlg('Gray White Matter selection','Segmentation');
set(h,'Pointer','watch');

bgt=min([max(max(Vols(1,:,:))),max(max(Vols(sizeh(1),:,:))),max(max(Vols(:,1,:))),max(max(Vols(:,sizeh(2),:))),max(max(Vols(:,:,1))),max(max(Vols(:,:,sizeh(3))))]),
cmax=[];
cmin=[];

[ScalpBoundary,BrainBoundary]=GetBoundaries(Volb,Vols,BGPARAMETER);
for i=10:3:sizeh(3)-10,
  I=Volb(:,:,i); c=improfile(I,repmat(sizeh(1)/2,sizeh(1),1),[1:1:sizeh(1)]);
  cmax=[cmax;c(maxima(c))];
  cmin=[cmin;c(maxima(-c))];
end;    

mmin=mean(cmin);
mmax=mean(cmax);
imminimaparameter=fix(mmax)-0.25*std(cmax);
otrop=(min(cmin)+2*std(cmin));

set(h,'Pointer','arrow');
close(h);


bb=0;
while ~bb
  kk=imextendedmin(Volb,imminimaparameter);  
  
  avw.img=kk<1;
  
  hh=avw_view(avw) 
  uiwait(hh)
  ans=questdlg('Do you like the White Segmentation Result?','Brain Segmentation Quality')
  switch ans,
    case 'No', 
      bb=0;
      prompt={['Enter new parameter. Increasing it let to thinner white matter']};
      def={num2str(imminimaparameter)};
      dlgTitle='Brain Strip Parameter';
      lineNo=1;
      answer=inputdlg(prompt,dlgTitle,lineNo,def);
      imminimaparameter=str2num(answer{1}),
      
    case 'Cancel',
      myerror('You do not like my program. What a pity')
    case 'Yes',
      bb=1;
  end % switch
end


white=Volb; white(kk==1)=0;  white=white(:);
gray=Volb>otrop; gray(kk==0)=0; gray=gray(:);
ScalpBoundary=ScalpBoundary(:);
BrainBoundary=BrainBoundary(:);
nVol=nVol(:); nVol(:)=0;

nVol(gray==1)=GRAY;
nVol(white>0)=WHITE;
nVol(BrainBoundary>0)=GRAYBORDER;
nVol(ScalpBoundary>0)=SKIN;

nVol=reshape(nVol,size(Volb));

avw.img=nVol; avw.hdr=head;

hh=avw_view(avw) 
uiwait(hh)

delete([pathname, 'segmented*.*']);

avw_write(avw,[ofilename,'.seg']);
%Save_Ana_file(nVol,head,[ofilename,'.seg'])


return





function [ScalpBoundary,BrainBoundary]=GetBoundaries(Volb,Vols,BGPARAMETER);

sizeh=size(Vols);
ll=(Volb>0); [N,X]=hist(double(Volb(ll))); X(1)=2000;
Vols(Vols<X(1))=0;
%Volb(Volb<X(1))=0;
h = waitbar(0,'Please wait.Simultaneous processing of scalp and brain extractor..');
for i=1:sizeh(3), 
  waitbar(i/sizeh(3),h)
  I = medfilt2(squeeze(Vols(:,:,i)),[4,4]); 
  Ifill = imfill(I,'holes'); 
  
  Ib=(Volb(:,:,i));
  
  [B,L] = bwboundaries(Ifill>BGPARAMETER,'noholes');
  [Bs,L] = bwboundaries(ll(:,:,i)>0,'noholes');
  nv=zeros(size(I));
  allb=[];
  if ~isempty(B),
    for k = 1:length(B)
      boundary = B{k};
      allb=[allb;boundary];
    end
    ind=sub2ind(size(I),allb(:,1),allb(:,2)); nv=nv(:); nv(ind)=1; nv=reshape(nv,size(I));
  end;  
  ScalpBoundary(:,:,i)=nv;
  
  nv=zeros(size(I));
  allb=[];
  if ~isempty(Bs),
    for k = 1:length(Bs)
      boundary = Bs{k};
      allb=[allb;boundary];
    end
    ind=sub2ind(size(I),allb(:,1),allb(:,2)); nv=nv(:); nv(ind)=1; nv=reshape(nv,size(I));    
  end
  BrainBoundary(:,:,i)=nv;   
  
end
close(h);
% Voy en la otra direccion
h = waitbar(0,'Please wait.Simultaneous processing of scalp and brain extractor..');

for i=1:sizeh(2), 
  waitbar(i/sizeh(2),h)
  I = medfilt2(squeeze(Vols(:,i,:)),[4,4]); 
  Ifill = imfill(I,'holes'); 
  
  Ib=squeeze(Volb(:,i,:));
  
  [B,L] = bwboundaries(Ifill>BGPARAMETER,'noholes');
  [Bs,L] = bwboundaries(Ib>X(1),'noholes');
  nv=zeros(size(I));
  allb=[];
  if ~isempty(B),
    for k = 1:length(B)
      boundary = B{k};
      allb=[allb;boundary];
    end
    ind=sub2ind(size(I),allb(:,1),allb(:,2)); nv=nv(:); nv(ind)=1; nv=reshape(nv,size(I));
  end;  
  sb=squeeze(ScalpBoundary(:,i,:)); sb=sb(:); nv=nv(:); sb(nv>0)=1; sb=reshape(sb,size(I));
  ScalpBoundary(:,i,:)=sb;
  
  nv=zeros(size(I));
  allb=[];
  if ~isempty(Bs),
    for k = 1:length(Bs)
      boundary = Bs{k};
      allb=[allb;boundary];
    end
    ind=sub2ind(size(I),allb(:,1),allb(:,2)); nv=nv(:); nv(ind)=1; nv=reshape(nv,size(I));    
  end
  sb=squeeze(BrainBoundary(:,i,:)); sb=sb(:); nv=nv(:); sb(nv>0)=1; sb=reshape(sb,size(I));
  BrainBoundary(:,i,:)=sb;
  
  
end
close(h);
