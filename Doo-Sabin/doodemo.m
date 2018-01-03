function doodemo
% DOODEMO  demonstrates the usage of DOO,PGONDISP,PGONTRACE
% See m-files for details.



% Setting paths
installpath=fileparts(which(mfilename));  % path of mfile this m-file

addpath(...                               % add paths to subfolders
       [installpath '\Colorfiles'],... 
       [installpath '\Data']);



% Definition of the polyhedral mesh to start with,
% first the points, then the polygons

pt=...          
    [-2 -2 -2;  % outer cube
      2 -2 -2;
      2  2 -2;
     -2  2 -2;
     -2 -2  2;
      2 -2  2;
      2  2  2;
     -2  2  2;
     -1 -1 -1;  % inner cube
      1 -1 -1;
      1  1 -1;
     -1  1 -1;
     -1 -1  1;
      1 -1  1;
      1  1  1;
     -1  1  1;
     -1 -2 -1;  % -y pierce
      1 -2 -1; 
      1 -2  1;
     -1 -2  1;
     -1 -1  2;  % +z pierce
      1 -1  2;
      1  1  2;
     -1  1  2; 
      2 -1 -1;  % +x pierce
      2  1 -1; 
      2  1  1;
      2 -1  1;
     -2 -1 -1;  % -x pierce  
     -2  1 -1;
     -2  1  1;    
     -2 -1  1;
     -1  2 -1;  % +y pierce    
      1  2 -1;
      1  2  1;
     -1  2  1;
     -1 -1 -2;  % -z pierce
      1 -1 -2;
      1  1 -2;
     -1  1 -2;];
   
pgon=...
   {[2 1 17 18]  % -y face
    [1 5 20 17]
    [5 6 19 20]
    [6 2 18 19]
    [2 6 28 25]  % +x face
    [6 7 27 28]
    [7 3 26 27]
    [3 2 25 26]
    [6 5 21 22]  % +z face
    [5 8 24 21]
    [8 7 23 24]
    [7 6 22 23]
    [5 1 29 32]   % -x face
    [8 5 32 31]
    [4 8 31 30]
    [1 4 30 29]
    [7 8 36 35]   % +y face
    [8 4 33 36]
    [4 3 34 33]
    [3 7 35 34]
    [1 2 38 37]   % -z face
    [2 3 39 38]
    [3 4 40 39]
    [4 1 37 40]
    [17 20 13  9] % -y pierce
    [20 19 14 13]
    [19 18 10 14]
    [18 17  9 10]
    [21 24 16 13] % +z pierce
    [24 23 15 16]
    [23 22 14 15]
    [22 21 13 14]
    [25 28 14 10] % +x pierce
    [28 27 15 14]
    [27 26 11 15]
    [26 25 10 11]
    [ 9 13 32 29] % -x pierce
    [13 16 31 32]
    [16 12 30 31]
    [12  9 29 30]
    [16 15 35 36] % +y pierce
    [15 11 34 35]
    [11 12 33 34]
    [12 16 36 33]
    [10  9 37 38] % -z pierce
    [11 10 38 39]
    [12 11 39 40]
    [9  12 40 37]};
 
 
fprintf('\n This may take some time. Please be patient. \n'); 
fprintf('\n Now calculating the Doo-Steps using DOO.\n')

[pt1,pgon1]=doo(pt,pgon,1,'orientation','negative','stats','off');
[pt2,pgon2]=doo(pt,pgon,2,'orientation','negative','stats','off');
[pt3,pgon3]=doo(pt,pgon,3,'orientation','negative','stats','off');

fprintf('\n Now displaying polygons using PGONDISP. \n');

colordef('black');
fig=figure('units','normalized','position',[0 0.4 0.49 0.54]);

subplot(2,2,1);
pgondisp(pt,pgon);
title('Initial Data','color',[1 0.7 0.5]);
drawnow;

subplot(2,2,2);
pgondisp(pt1,pgon1);
title('Doo Step 1','color',[1 0.7 0.5]);
drawnow;

subplot(2,2,3);
pgondisp(pt2,pgon2);
title('Doo Step 2','color',[1 0.7 0.5]);
drawnow;

subplot(2,2,4);
pgondisp(pt3,pgon3);
title('Doo Step 3','color',[1 0.7 0.5]);
drawnow;

fprintf('\n Now raytracing polygons using PGONTRACE. \n');

fig=figure('units','normalized','position',[0.5 0.4 0.49 0.54]);
drawnow;
pgontrace(pt3,pgon3,'sine');
li=light('pos',[-0.4 -0.2 -1]); 
 
 
 
