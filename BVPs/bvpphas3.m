function stop=bvpphas3(t,y,flag)

%   BVPPHAS3 output function for BVP and BVP2
%
%   Plots a 3-dimensional phase portrait of the current approximation
%
%   When the string 'bvpphas3' is passed to a BVP solver as the 'OutputFcn'
%   property, i.e. options = bvpset('OutputFcn','bvpphas3'), the solver calls
%   BVPPHAS3(T,Y,FLAG) after every iteration (See BVPSET for details).
%   To plot only particular components, specify their indices in the
%   'OutputSel' property passed to the BVP solver(again BVPSET for details).
%   
%   See also BVPPLOT  BVPPHAS2  BVPSET  BVP  BVP2

%  Copyright (c) 1999 Guenter Kneisl
%                     University of Technology, Vienna
%                     e9425595@fbma.tuwien.ac.at

if isempty(flag)
   usrdata=get(gcf,'userdata');
   stop=usrdata(1);
   li=usrdata(2);

   set(li,'xdata',y(:,1),'ydata',y(:,2),'zdata',y(:,3));
   drawnow;
else
   if flag=='init'
      f=figure(gcf);
      rotate3d on;
      li=line(y(:,1),y(:,2),y(:,3),'erasemode','normal','Linestyle','-');
   
      h = findobj(f,'Tag','stop');
     
      if isempty(h)
        stop = 0;
        pos = get(0,'DefaultUicontrolPosition');
        pos(1) = pos(1) - 15;
        pos(2) = pos(2) - 15;
        str = 'usrdata=get(gcf,''userdata''); stop=1; set(gcf,''UserData'',[stop usrdata(2)]);';
        uicontrol( ...
             'Style','push', ...
             'String','Stop', ...
             'Position',pos, ...
             'Callback',str, ...
             'Tag','stop');
      else
         set(h,'Visible','on');            % make sure it's visible
         stop = get(f,'UserData');
      end
      set(f,'UserData',[0 li]);
      pause;
   else  % flag=='done'
      stop=0;
   end
end

