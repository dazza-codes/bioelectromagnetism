function stop=bvpplot(t,y,flag)

%   BVPPLOT output function for BVP and BVP2
%
%   Plots the current approximation vector Y versus T.
%
%   When the string 'bvpplot' is passed to a BVP solver as the 'OutputFcn'
%   property, i.e. options = bvpset('OutputFcn','bvpplot'), the solver calls
%   BVPPLOT(T,Y,FLAG) after every iteration (See BVPSET for details).
%   To plot only particular components, specify their indices in the
%   'OutputSel' property passed to the BVP solver (again BVPSET for details).
%   BVPPLOT is the default output function of the solvers when they are
%   called without output arguments.
%   
%   See also BVPPHAS2  BVPPHAS3  BVPSET  BVP  BVP2

%  Copyright (c) 1999 Guenter Kneisl
%                     University of Technology, Vienna
%                     e9425595@fbma.tuwien.ac.at

if isempty(flag)
   usrdata=get(gcf,'userdata');
   stop=usrdata{1};
   li  =usrdata{2};
   
   for i=1:length(li)
      set(li(i),'ydata',y(:,i));
   end
   
   drawnow;
else
   if flag=='init'
      f=figure(gcf);
      rotate3d on;
      li=line(t,y,'erasemode','normal','Linestyle','-');
   
      h = findobj(f,'Tag','stop');
     
      if isempty(h)
        stop = 0;
        pos = get(0,'DefaultUicontrolPosition');
        pos(1) = pos(1) - 15;
        pos(2) = pos(2) - 15;
        str = 'usrdata=get(gcf,''userdata''); stop=1; set(gcf,''UserData'',{stop usrdata{2}});';
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
      set(f,'UserData',{0 li});
      pause;
   else  % flag=='done'
      stop=0;
   end
end

