function poi(slice_no,selection);

% File/function name:  poi.m   (Pixels Of Interest)
%         Written by:  Tim DeMonte
%               Date:  Mar 27, 2001
%            Version:  1.0
%        Description:  A tool for studying Regions of Interest (ROI)
%                      defined by selected groups of pixels in 2D
%                      images taken from multi-slice (3D) data sets.
%

global image Jmax Jtotal doublemap
global ROI_1 ROI_1sum ROI_1no ROI_1mean ROI_1std ROI_1min ROI_1max
global ROI_1percent
global h1 h2
global rad1 rad2 rad3
global editbox1 editbox2 editbox3 editbox4 editbox5 editbox6 editbox7 editbox8
if nargin<2,
   selection=0;
   clear global ROI_1 ROI_1mean ROI_1std ROI_1min ROI_1max ROI_1percent
end

if selection==0,
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Modify these lines to read in a particular image file               %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   load Jxyz.mat;
   image=Jz(:,:,slice_no);
   Jtotal=sum(sum(image(:,:)));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Modify this line to specify the maximum range of data to display    %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Jmax=15;
   
   ROI_1mean=0;
   ROI_1std=0;
   ROI_1min=0;
   ROI_1max=0;
   ROI_1sum=0;
   ROI_1no=0;
   dfc=get(0,'defaultfigurecolor');
   gui_dim=[200,100,820,625];
   figure('name',...
      'POI (Pixels of Interest)',...
      'numbertitle','off','position',gui_dim,'menubar','none');
   h1=imagesc(image(:,:),[-Jmax Jmax]);
   axis square
   pwr=0.5; % exponent for colormap (1 is linear)
   for m=1:128
      if m<64
         doublemap(m,1)=((m-1)./64).^pwr;
         doublemap(m,2)=((m-1)./64).^pwr;
         doublemap(m,3)=1;
      elseif m>64
         doublemap(m,1)=1;
         doublemap(m,2)=((128-m)./64).^pwr;
         doublemap(m,3)=((128-m)./64).^pwr;
      else
         doublemap(m,1)=1;
         doublemap(m,2)=1;
         doublemap(m,3)=1;
      end
   end
   colormap(doublemap)
   colorbar
   set(h1,'buttondownfcn','global slice_no,poi(slice_no,1),clear slice_no');
   
   % Text
   uicontrol('style','text','position',[5,440,80,20],...
      'string','Mean of ROI:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   uicontrol('style','text','position',[5,390,80,20],...
      'string','Std Dev of ROI:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   uicontrol('style','text','position',[5,90,80,20],...
      'string','Last Pick Value:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   uicontrol('style','text','position',[5,340,80,20],...
      'string','Min of ROI:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   uicontrol('style','text','position',[5,290,80,20],...
      'string','Max of ROI:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   uicontrol('style','text','position',[5,510,90,20],...
      'string','Sum of ROI:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   uicontrol('style','text','position',[5,210,80,20],...
      'string','No. of Pixels:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   uicontrol('style','text','position',[5,160,80,20],...
      'string','% of Total:','backgroundcolor',dfc,'foregroundcolor','black',...
      'fontname','ms sans serif','fontsize',8,'fontunits','points',...
      'horizontalalignment','left');
   
   % Edit Boxes
   editbox1=uicontrol('style','edit','position',[5,420,90,20],...
      'string',num2str(ROI_1mean),'horizontalalignment','left',...
      'backgroundcolor','white');
   editbox2=uicontrol('style','edit','position',[5,370,90,20],...
      'string',num2str(ROI_1std),'horizontalalignment','left',...
      'backgroundcolor','white');
   editbox3=uicontrol('style','edit','position',[5,70,90,20],...
      'string','','horizontalalignment','left',...
      'backgroundcolor','white');
   editbox4=uicontrol('style','edit','position',[5,320,90,20],...
      'string',num2str(ROI_1min),'horizontalalignment','left',...
      'backgroundcolor','white');
   editbox5=uicontrol('style','edit','position',[5,270,90,20],...
      'string',num2str(ROI_1max),'horizontalalignment','left',...
      'backgroundcolor','white');
   editbox6=uicontrol('style','edit','position',[5,490,90,20],...
      'string',num2str(ROI_1sum),'horizontalalignment','left',...
      'backgroundcolor','white');
   editbox7=uicontrol('style','edit','position',[5,190,90,20],...
      'string',num2str(ROI_1no),'horizontalalignment','left',...
      'backgroundcolor','white');
   editbox8=uicontrol('style','edit','position',[5,140,90,20],...
      'string','','horizontalalignment','left',...
      'backgroundcolor','white');
   
   % Radio Buttons
   rad1=uicontrol('style','radiobutton','position',[5,585,90,20],...
      'string','Pick/Clear','min',0,'max',1,'value',1);
   rad2=uicontrol('style','radiobutton','position',[5,560,90,20],...
      'string','Big Block','min',0,'max',1,'value',0);
   rad3=uicontrol('style','radiobutton','position',[5,535,90,20],...
      'string','Zoom ON/OFF','min',0,'max',1,'value',0,...
      'callback','global slice_no,poi(slice_no,2),clear slice_no');
   
elseif selection==1,
   pos=round(get(gca,'currentpoint'));
   x=[pos(1,1)-0.5 pos(1,1)-0.5;pos(1,1)-0.5 pos(1,1)+0.5;pos(1,1)+0.5 pos(1,1)+0.5;...
         pos(1,1)+0.5 pos(1,1)-0.5;pos(1,1)-0.5 pos(1,1)];
   y=[pos(1,2)-0.5 pos(1,2)+0.5;pos(1,2)+0.5 pos(1,2)+0.5;pos(1,2)+0.5 pos(1,2)-0.5;...
         pos(1,2)-0.5 pos(1,2)-0.5;pos(1,2)-0.5 pos(1,2)+0.5];
   if(get(rad1,'value')),   % Pick
      if(get(rad2,'value')),  % Big Block (i.e. 3x3 array)
         h3=line(x,y);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)==ROI_1(:,2));
            b=find(pos(1,2)==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
         end
         %%% End of Check ... %%%
         h3=line(x,y-1);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)==ROI_1(:,2));
            b=find(pos(1,2)-1==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)) pos(1,1) pos(1,2)-1];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)) pos(1,1) pos(1,2)-1];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)) pos(1,1) pos(1,2)-1];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)) pos(1,1) pos(1,2)-1];
         end
         %%% End of Check ... %%%
         
         h3=line(x+1,y-1);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)+1==ROI_1(:,2));
            b=find(pos(1,2)-1==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)+1) pos(1,1)+1 pos(1,2)-1];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)+1) pos(1,1)+1 pos(1,2)-1];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)+1) pos(1,1)+1 pos(1,2)-1];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)+1) pos(1,1)+1 pos(1,2)-1];
         end
         %%% End of Check ... %%%
         
         h3=line(x-1,y);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)-1==ROI_1(:,2));
            b=find(pos(1,2)==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)-1) pos(1,1)-1 pos(1,2)];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)-1) pos(1,1)-1 pos(1,2)];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2),pos(1,1)-1) pos(1,1)-1 pos(1,2)];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2),pos(1,1)-1) pos(1,1)-1 pos(1,2)];
         end
         %%% End of Check ... %%%
         
         h3=line(x+1,y);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)+1==ROI_1(:,2));
            b=find(pos(1,2)==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)+1) pos(1,1)+1 pos(1,2)];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)+1) pos(1,1)+1 pos(1,2)];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2),pos(1,1)+1) pos(1,1)+1 pos(1,2)];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2),pos(1,1)+1) pos(1,1)+1 pos(1,2)];
         end
         %%% End of Check ... %%%
         
         h3=line(x-1,y+1);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)-1==ROI_1(:,2));
            b=find(pos(1,2)+1==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)-1) pos(1,1)-1 pos(1,2)+1];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)-1) pos(1,1)-1 pos(1,2)+1];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)-1) pos(1,1)-1 pos(1,2)+1];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)-1) pos(1,1)-1 pos(1,2)+1];
         end
         %%% End of Check ... %%%
         
         h3=line(x,y+1);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)==ROI_1(:,2));
            b=find(pos(1,2)+1==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)) pos(1,1) pos(1,2)+1];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)) pos(1,1) pos(1,2)+1];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)) pos(1,1) pos(1,2)+1];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)) pos(1,1) pos(1,2)+1];
         end
         %%% End of Check ... %%%
         
         h3=line(x+1,y+1);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)+1==ROI_1(:,2));
            b=find(pos(1,2)+1==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)+1) pos(1,1)+1 pos(1,2)+1];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)+1) pos(1,1)+1 pos(1,2)+1];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)+1) pos(1,1)+1 pos(1,2)+1];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2)+1,pos(1,1)+1) pos(1,1)+1 pos(1,2)+1];
         end
         %%% End of Check ... %%%
         
         h3=line(x-1,y-1);
         set(h3,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)-1==ROI_1(:,2));
            b=find(pos(1,2)-1==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)-1) pos(1,1)-1 pos(1,2)-1];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)-1) pos(1,1)-1 pos(1,2)-1];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)-1) pos(1,1)-1 pos(1,2)-1];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2)-1,pos(1,1)-1) pos(1,1)-1 pos(1,2)-1];
         end
         %%% End of Check ... %%%
         
         ROI_1mean=mean(ROI_1(:,1));
         ROI_1std=std(ROI_1(:,1));
         ROI_1min=min(ROI_1(:,1));
         ROI_1max=max(ROI_1(:,1));
         ROI_1sum=sum(ROI_1(:,1));
         [ROI_1no,n]=size(ROI_1);
         ROI_1percent=(ROI_1sum./Jtotal).*100;
         set(editbox1,'string',num2str(ROI_1mean));
         set(editbox2,'string',num2str(ROI_1std));
         set(editbox3,'string',num2str(image(pos(1,2),pos(1,1))));
         set(editbox4,'string',num2str(ROI_1min));
         set(editbox5,'string',num2str(ROI_1max));
         set(editbox6,'string',num2str(ROI_1sum));
         set(editbox7,'string',num2str(ROI_1no));
         set(editbox8,'string',num2str(ROI_1percent));
      else   % single pixel
         h2=line(x,y);
         set(h2,'color','g');
         %%% Check if pixel is already selected %%%
         if ~isempty(ROI_1),
            a=find(pos(1,1)==ROI_1(:,2));
            b=find(pos(1,2)==ROI_1(:,3));
            if isempty(a),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
            elseif isempty(b),
               ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
            else
               c=0;
               for n=1:length(a),
                  if find(a(n)==b(:)),
                     c=a(n);
                  end
               end
               if c==0,
                  ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
               end
            end
         else
            ROI_1=[ROI_1;image(pos(1,2),pos(1,1)) pos(1,1) pos(1,2)];
         end
         %%% End of Check ... %%%
         ROI_1mean=mean(ROI_1(:,1));
         ROI_1std=std(ROI_1(:,1));
         ROI_1min=min(ROI_1(:,1));
         ROI_1max=max(ROI_1(:,1));
         ROI_1sum=sum(ROI_1(:,1));
         [ROI_1no,n]=size(ROI_1);
         ROI_1percent=(ROI_1sum./Jtotal).*100;
         set(editbox1,'string',num2str(ROI_1mean));
         set(editbox2,'string',num2str(ROI_1std));
         set(editbox3,'string',num2str(image(pos(1,2),pos(1,1))));
         set(editbox4,'string',num2str(ROI_1min));
         set(editbox5,'string',num2str(ROI_1max));
         set(editbox6,'string',num2str(ROI_1sum));
         set(editbox7,'string',num2str(ROI_1no));
         set(editbox8,'string',num2str(ROI_1percent));
      end
   else   % Clear
      pix=image(pos(1,2),pos(1,1));
      color=round(((pix+Jmax)./(2.*Jmax)).*128);
      if color>128,
         color=128;
      elseif color<1,
         color=1;
      end
      if(get(rad2,'value')),  % Clear Big Block
         h3=line(x,y);
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)==ROI_1(:,2));
         b=find(pos(1,2)==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x,y-1);
         pix=image(pos(1,2)-1,pos(1,1));
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)==ROI_1(:,2));
         b=find(pos(1,2)-1==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x+1,y-1);
         pix=image(pos(1,2)-1,pos(1,1)+1);
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)+1==ROI_1(:,2));
         b=find(pos(1,2)-1==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x-1,y);
         pix=image(pos(1,2),pos(1,1)-1);
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)-1==ROI_1(:,2));
         b=find(pos(1,2)==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x+1,y);
         pix=image(pos(1,2),pos(1,1)+1);
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)+1==ROI_1(:,2));
         b=find(pos(1,2)==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x-1,y+1);
         pix=image(pos(1,2)+1,pos(1,1)-1);
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)-1==ROI_1(:,2));
         b=find(pos(1,2)+1==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x,y+1);
         pix=image(pos(1,2)+1,pos(1,1));
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)==ROI_1(:,2));
         b=find(pos(1,2)+1==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x+1,y+1);
         pix=image(pos(1,2)+1,pos(1,1)+1);
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)+1==ROI_1(:,2));
         b=find(pos(1,2)+1==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         h3=line(x-1,y-1);
         pix=image(pos(1,2)-1,pos(1,1)-1);
         color=round(((pix+Jmax)./(2.*Jmax)).*128);
         if color>128,
            color=128;
         elseif color<1,
            color=1;
         end
         set(h3,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)-1==ROI_1(:,2));
         b=find(pos(1,2)-1==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         ROI_1mean=mean(ROI_1(:,1));
         ROI_1std=std(ROI_1(:,1));
         ROI_1min=min(ROI_1(:,1));
         ROI_1max=max(ROI_1(:,1));
         ROI_1sum=sum(ROI_1(:,1));
         [ROI_1no,n]=size(ROI_1);
         ROI_1percent=(ROI_1sum./Jtotal).*100;
         set(editbox1,'string',num2str(ROI_1mean));
         set(editbox2,'string',num2str(ROI_1std));
         set(editbox3,'string',num2str(image(pos(1,2),pos(1,1))));
         set(editbox4,'string',num2str(ROI_1min));
         set(editbox5,'string',num2str(ROI_1max));
         set(editbox6,'string',num2str(ROI_1sum));
         set(editbox7,'string',num2str(ROI_1no));
         set(editbox8,'string',num2str(ROI_1percent));
      else  % Clear single pixel
         h2=line(x,y);
         set(h2,'color',doublemap(color,:));
         %%% Remove a pixel from ROI list %%%         
         a=find(pos(1,1)==ROI_1(:,2));
         b=find(pos(1,2)==ROI_1(:,3));
         if ~isempty(a) & ~isempty(b),
            c=0;
            for n=1:length(a),
               if find(a(n)==b(:)),
                  c=a(n);
               end
            end
            if c,
               d=size(ROI_1);
               ROI_1(c,1:3)=ROI_1(d(1),1:3);
               ROI_1=ROI_1(1:(d(1)-1),:);
            end
         end
         %%% End of Remove ... %%%
         ROI_1mean=mean(ROI_1(:,1));
         ROI_1std=std(ROI_1(:,1));
         ROI_1min=min(ROI_1(:,1));
         ROI_1max=max(ROI_1(:,1));
         ROI_1sum=sum(ROI_1(:,1));
         [ROI_1no,n]=size(ROI_1);
         ROI_1percent=(ROI_1sum./Jtotal).*100;
         set(editbox1,'string',num2str(ROI_1mean));
         set(editbox2,'string',num2str(ROI_1std));
         set(editbox3,'string',num2str(image(pos(1,2),pos(1,1))));
         set(editbox4,'string',num2str(ROI_1min));
         set(editbox5,'string',num2str(ROI_1max));
         set(editbox6,'string',num2str(ROI_1sum));
         set(editbox7,'string',num2str(ROI_1no));
         set(editbox8,'string',num2str(ROI_1percent));
      end
   end
elseif selection==2,
   if(get(rad3,'value')),
      zoom
   else
      zoom off
   end
end