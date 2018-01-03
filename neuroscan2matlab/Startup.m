%
%STARTUP Startup for utilities.
%	STARTUP Setup some needed globals and looks for utilities
%	toolboxes installed.
%
%DIAGNOSTICS
%	Used global variables are DIRS, DIRSS and PCWIN.
%
%EXAMPLES
%	startup

%Mention source when using or modifying these Shareware tools

%JVIR,
%JVIR,Jussi.Virkkala@occuphealth.fi

%JVIR, 4-Feb-1999 Changed folder names.
%JVIR, 2-Feb-1999 Modified for PCWIN Matlab 5.2.

%	J.Virkkala 5-Apr-95
%	J.Virkkala 4-May-95 Started to make to work with IBM RS.

%JVIR,
close all
clear all 
clear globals

set(0,'defaultfigurepapertype','a4')

%JVIR,
global BORDERX BORDERY BORDERB
BORDERY=100;
BORDERX=100;
BORDERB=100;

% directory separator
clear all
global DIRS DIRSS PCWIN
if strcmp(computer,'PCWIN')==1,
  DIRS='\';
  DIRSS=';';
  PCWIN=1;
else
  DIRS='/';
  DIRSS=':';
  PCWIN=0;
end
	% looking existing utilities
p=path;
op=cd;
i1=fliplr(find(path==DIRSS));
i2=fliplr(find(path==DIRS));
p=p(i1(1)-1:i2(1)-1);

%JVIR, Matlab and toolbox in different drives
%JVIR, Currently only scanutil
p(1)=op(1);
s=['scanutil';'aumautil';'edsautil'];
for i=1:size(s,1),
  tn=s(i,:);
  eval(['cd ' p DIRS tn],'tn=[];');
  if ~isempty(tn),
	% path is different between PCWIN ja IBM RS
    path(path,[p DIRS tn]);
  end
end 

%JVIR, path can include spaces
eval(['cd(''' op ''');']);
disp(' ');
ver
path(path,cd);

eval('setscan','');
eval('setauma','');
eval('setedsa','');

%END OF STARTUP