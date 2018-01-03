%------------------------------------------------------------------- 

function filename=ZeroPad(prefix,npos,fileno) 

% function filename=ZeroPad(path,npos,fileno).m 
% prefix e.g. path and start of filename (string) 
% npos total no of positions (int) 
% fileno file number (int) 
% 
% by gh may 1999 

if fileno >= 10^(npos-1) %no zeropadding needed 
  filename = [prefix num2str(fileno)]; 
else 
  for i=1:npos-1 
    if (fileno < (10^i)) & (exist('filename') == 0) 
       for k=1:npos-i 
         zeropd(k) = '0'; 
       end; 
       filename = [prefix zeropd num2str(fileno)]; 
    end 
  end 
end 

