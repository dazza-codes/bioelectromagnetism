function stop=bvpgunplot(t,y,flag)

% Output function for BVP used by BVPDEMO
%
% See also BVPPHAS2, BVPPHAS3, BVPSET, BVP


n=y(:,1);
w=y(:,2);
E=y(:,3);

if isempty(flag)
   usrdata=get(gcf,'userdata');
   stop=usrdata(1);
   li=usrdata(2);
   set(li,'xdata',E,'ydata',n,'zdata',w);
   drawnow;
else   
   if flag=='init'
      f=figure(gcf);
      rotate3d on;
      li=line(E,n,w,'erasemode','normal','Linestyle','-','marker','o');
      
      h = findobj(f,'Tag','stop');
     
      if isempty(h)
        stop = 0;
        pos = get(0,'DefaultUicontrolPosition');
        pos(1) = pos(1) - 15;
        pos(2) = pos(2) - 15;
        %str = 'stop=get(gcf,''userdata''); stop=1; set(gcf,''UserData'',stop);';
        str = 'usdata=get(gcf,''UserData'');stop=1;set(gcf,''UserData'',[stop usdata(2)]);';

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
   else   % flag=='done'
      stop=0;
   end
end

