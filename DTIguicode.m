tic



getDicomsgui           %uses my script to load an entire Dicom directory.

disp('Calculating tensor...')

noise=evalin('base','noise');                                   %noise threshold
doFAmap=evalin('base','doFAmap');                               %write DICOMs?
dotADC=evalin('base','dotADC');
dolambda1=evalin('base','dolambda1');
dolambda2=evalin('base','dolambda2');
dolambda3=evalin('base','dolambda3');
docm=evalin('base','docm');
counter=evalin('base','counter'); %counter comes from getDicomsgui
extension=evalin('base','extension');


slices=counter/7;


FAmaps=[];
tADCmaps=[];
lambda1_maps=[];
lambda2_maps=[];
lambda3_maps=[];
cmaps=[];
xyzmaps=[];







warning off MATLAB:divideByZero                                 %supress matlab warning


for numslice=1:slices

bo=tensorvol(:,:,numslice,1);%bo=bo';                            %transposes the images to correct for getDicoms
b1=tensorvol(:,:,numslice,2);%b1=b1'; 
b2=tensorvol(:,:,numslice,3);%b2=b2'; 
b3=tensorvol(:,:,numslice,4);%b3=b3'; 
b4=tensorvol(:,:,numslice,5);%b4=b4'; 
b5=tensorvol(:,:,numslice,6);%b5=b5'; 
b6=tensorvol(:,:,numslice,7);%b6=b6'; 


ADC101=log(im2double(b1)./im2double(bo))./(-1000);              %Calculate ADC maps for each direction
ADC_101=log(im2double(b2)./im2double(bo))./(-1000);
ADC011=log(im2double(b3)./im2double(bo))./(-1000);
ADC01_1=log(im2double(b4)./im2double(bo))./(-1000);
ADC110=log(im2double(b5)./im2double(bo))./(-1000);
ADC_110=log(im2double(b6)./im2double(bo))./(-1000);

transform=[1 1 0 2 0 0;1 0 1 0 2 0;0 1 1 0 0 2;1 1 0 -2 0 0;1 0 1 0 -2 0;0 1 1 0 0 -2];  %the transformation matrix to obtain tensor elements

xx=zeros(256);                                                  %Create empty tensor elements
yy=zeros(256);
zz=zeros(256);
xy=zeros(256);
xz=zeros(256);
yz=zeros(256);
lambda1=zeros(256);                                             %create empty images
lambda2=zeros(256);
lambda3=zeros(256);
r=zeros(256);
g=zeros(256);
b=zeros(256);
FAmap=zeros(256);
tADC=zeros(256);

for i=1:256                                                     %Fill in tensor elements pixel by pixel. 
    for j=1:256
        if (bo(i,j)>noise)                                      %noise threshold.  Change this to zero if you want the whole image to be calculated (takes longer)
            ADCm=[ADC110(i,j);ADC101(i,j);ADC011(i,j);ADC_110(i,j);ADC_101(i,j);ADC01_1(i,j)]; %see Hedehus, http://www-radiology.stanford.edu/majh/
            ADCe=inv(transform)*ADCm;
            xx(i,j)=ADCe(1,1);                                  %get tensor elements
            yy(i,j)=ADCe(2,1);
            zz(i,j)=ADCe(3,1);
            xy(i,j)=ADCe(4,1);
            xz(i,j)=ADCe(5,1);
            yz(i,j)=ADCe(6,1);
            ten=[xx(i,j) xy(i,j) xz(i,j);xy(i,j) yy(i,j) yz(i,j);xz(i,j) yz(i,j) zz(i,j)];  %the tensor itself
            [V,D]=eig(ten);
            D=eig(ten);                                         %get eigenvalues
            D=abs(D);                                           
            E=D;                                                %get the eigenvalues before re-ordering
            D=sort(D);                                          %sort the eigenvalues (upwards)
            e1=D(3,1);e1=e1*2;
            e2=D(2,1);e2=e2*2;
            e3=D(1,1);e3=e3*2;
            lambda1(i,j)=e1;
            lambda2(i,j)=e2;
            lambda3(i,j)=e3;
            %V=abs(V);       %get the x,y,z components of e3   
            largest=find(E==D(3,1));                            %find the value of the largest eigenvalue in the original matrix
            sz=length(largest);
            r(i,j)=V(1,largest(sz));                                %in order to get its corresponding eigenvector components in x, y and z
            g(i,j)=V(2,largest(sz));
            b(i,j)=V(3,largest(sz));
            trace=(e1+e2+e3)/3;
            %FA=(sqrt(3/2)*sqrt((1/3)*((e1-e2)^2+(e2-e3)^2+(e3-e1)^2)))/(sqrt(e1^2+e2^2+e3^3)); %Another formula for FA
            FA=(sqrt(3*((e1-trace)^2+(e2-trace)^2+(e3-trace)^2)))/(sqrt(2*(e1^2+e2^2+e3^2)));  %Le Bihan 2001
            tADC(i,j)=trace;
            FAmap(i,j)=FA;
        end;
    end;
end;




rgb=cat(3,r,g,b);
xyz=cat(3,r,(g*-1),b);
FAmap3=cat(3,FAmap,FAmap,FAmap);                                %needed to combine with the colormap, need 3D.
cm=abs(rgb).*FAmap3;                                            %FA weighting for colormap.  Has absolute values for vector components

%eval (['rgb' num2str(numslice) '=rgb;']);
%nameofmap=['rgb' num2str(numslice)];
%assignin('base',nameofmap,rgb);


FAmaps=cat(3,FAmaps,FAmap);
tADCmaps=cat(3,tADCmaps,tADC);
lambda1_maps=cat(3,lambda1_maps,lambda1);
lambda2_maps=cat(3,lambda2_maps,lambda2);
lambda3_maps=cat(3,lambda3_maps,lambda3);
cmaps=cat(4,cmaps,cm);
xyzmaps=cat(4,xyzmaps,xyz);






%eval (['xyz' num2str(numslice) '=rgb;']);                       %xyz is the same as cm (colormap) but without abs.  But, apparently, doesn't make a difference.     
%nameofmap=['xyz' num2str(numslice)];
%assignin('base',nameofmap,rgb);

%image (cm);axis image
%eval (['FA' num2str(numslice) '=FAmap;']);
%nameofmap=['FA' num2str(numslice)];
%assignin('base',nameofmap,FAmap);
if (doFAmap==1)
    number=num2str(numslice);                                   %FOUR LINES FOR FILE WRITING, REPEAT FOR THE OTHER IMAGES
    st1='FA';st2='.dcm';
    nameforfile=[st1 number st2];
    dicomwrite(FAmap,nameforfile);
    disp('writing FA map to DICOM file')
end
%eval (['tADC' num2str(numslice) '=tADC;']);
%nameofmap=['tADC' num2str(numslice)];
%assignin('base',nameofmap,tADC);
if (dotADC==1)
    number=num2str(numslice);                                   %FOUR LINES FOR FILE WRITING, REPEAT FOR THE OTHER IMAGES
    st1='tADC';st2='.dcm';
    nameforfile=[st1 number st2];
    dicomwrite(tADC,nameforfile);
    disp('writing tADC map to DICOM file')
end
%eval (['lambda1' num2str(numslice) '=lambda1;']);
%nameofmap=['lambda1_' num2str(numslice)];
%assignin('base',nameofmap,lambda1);
if (dolambda1==1)
    number=num2str(numslice);                                   %FOUR LINES FOR FILE WRITING, REPEAT FOR THE OTHER IMAGES
    st1='lambda1_';st2='.dcm';
    nameforfile=[st1 number st2];
    dicomwrite(lambda1,nameforfile);
    disp('writing lambda1 map to DICOM file')
end
%eval (['lambda2' num2str(numslice) '=lambda2;']);
%nameofmap=['lambda2_' num2str(numslice)];
%assignin('base',nameofmap,lambda2);
if (dolambda2==1)
    number=num2str(numslice);                                   %FOUR LINES FOR FILE WRITING, REPEAT FOR THE OTHER IMAGES
    st1='lambda2_';st2='.dcm';
    nameforfile=[st1 number st2];
    dicomwrite(lambda2,nameforfile); 
    disp('writing lambda2 map to DICOM file')
end
%eval (['lambda3' num2str(numslice) '=lambda3;']);
%nameofmap=['lambda3_' num2str(numslice)];
%assignin('base',nameofmap,lambda3);
if (dolambda3==1)
    number=num2str(numslice);                                   %FOUR LINES FOR FILE WRITING, REPEAT FOR THE OTHER IMAGES
    st1='lambda3_';st2='.dcm';
    nameforfile=[st1 number st2];
    dicomwrite(lambda3,nameforfile);
    disp('writing lambda3 map to DICOM file')
end
%eval (['cm' num2str(numslice) '=cm;']);
%nameofmap=['cm' num2str(numslice)];
%assignin('base',nameofmap,cm);
if (docm==1)
    number=num2str(numslice);                                   %FOUR LINES FOR FILE WRITING, REPEAT FOR THE OTHER IMAGES
    st1='cm';st2='.dcm';
    nameforfile=[st1 number st2];
    dicomwrite(cm,nameforfile); 
    disp('writing colormap to DICOM file')
end

mes=[num2str(numslice) ' of ' num2str(slices) ' finished'];
disp(mes)
end
 
%assignin('base','ten',ten)
clear FA FAmap tADC lambda1 lambda2 lambda3 cm FAmap3 rgb bo b1 b2 b3 b4 b5 b6 ADC101 ADC_101 ADC011;
clear ADC01_1 ADC110 ADC_110 r g b xx yy zz xy xz yz trace ten tensorvol transform ADCe ADCm D V;
assignin('base','tensorvol',0)


assignin('base','FAmaps',FAmaps);
assignin('base','tADCmaps',tADCmaps);
assignin('base','lambda1_maps',lambda1_maps);
assignin('base','lambda2_maps',lambda2_maps);
assignin('base','lambda3_maps',lambda3_maps);
assignin('base','cmaps',cmaps);
assignin('base','xyzmaps',xyzmaps);

toc
