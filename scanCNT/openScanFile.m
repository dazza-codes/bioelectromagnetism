% fid = openScanFile('filename')
%               
%   filename    <-  Filename string
%   fid         <-  File identifier
%
function fid = openScanFile(filename)

fid = fopen(filename,'r','l');