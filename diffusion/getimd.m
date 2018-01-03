function [nvox,nim,imd]=getimd(deck)
%
%
%
%  [nvox,nim,imd]=getimd(deck);
%
%
%
%  ER, DRCMR, Wed Jan  7 20:49:23 MET 1998
%



n=ndims(deck);
imd=size(deck);
imd(end+1:4)=1;
switch n
case 2  
  nvox=prod(imd);
  nim=1;
case 3
  nvox=prod(imd(1:2));
  nim=imd(3);
case 4
  nvox=prod(imd(1:3));
  nim=imd(4);
end  


