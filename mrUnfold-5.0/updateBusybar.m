function newProgStat=updateBusybar(handle,currentProgressStatus)
%function progressStatus=updateProgressbar(handle,currentProgressStatus)
%GUIs can have a text box that acts as a progress bar. This routine updates 
% one by adding a new character to the box and zeroing it at the end.
posit=get(handle,'Position');

boxLen=fix(posit(3));
newProgStat=mod(currentProgressStatus+1,boxLen);
set(handle,'String',repmat('|',1,newProgStat));
