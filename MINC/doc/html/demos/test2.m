h=openimage('/local/matlab/toolbox/local/examples/yates_19445.mnc');
h2 = newimage('new.mnc',[0 15], ...
     '/local/matlab/toolbox/local/examples/yates_19445.mnc');

ftimes = getimageinfo (h,'MidFrameTimes');

img = zeros (16384,1);

for j=1:15
   PET = getimages(h,j,1:21,PET);
   for i=1:128:16257;
      line = PET(i:i+127,:);
      img (i:i+127) = ntrapz(ftimes, line')';
   end;

   putimages(h2,img,j);
   disp('Generated image.');
end;

closeimage(h);
closeimage(h2);
