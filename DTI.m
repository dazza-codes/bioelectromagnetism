%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%TENSORCALC2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Enzyme.  Feb, 2004%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matlab implementation of the instructions given by Dr. Maj Hedehus 
%in:  http://www-radiology.stanford.edu/majh/
%This program calculates the tensor for an MRI slice and creates maps for
%Fractional Anisotropy (FAmap), traceADC (tADC), principal eigenvector 
%(lambda1) and its orthogonals (lambda2 and lambda3), and an FA-weighted
%colormap for the principal eigenvector (cm).
%It also creates the tensor elements (xx,yy,zz,xy,yz,xz).
%
%
%HOW TO USE:
%
%1. Load the image slices to variables bo (note it is bo, letter 'o'), 
%   b1, b2, b3, b4, b5 and b6. Use the function dicomread to do this i.e: 
%       >>bo=dicomread('image.dcm');
%2. Run this script 
%       >>DTI
%3. The maps and images are created with the names as described above.  To 
%   view a particular map or image type 
%       >>imagesc(FAmap);colormap(gray);axis image;
%
%
%
%Please note this program requires the image processing toolbox with dicomread
%function.  If you do not have it, you can load raw data images with the 'load'
%function, availabe in any matlab.
%
%
%VERY IMPORTANT!!!!   This program assumes you only have one b-value and it is 
%   set to b=1000.  
%Also:  feel free to make any modifications to this code.  
%
%ENJOY!



       


ADC101=log(im2double(b1)./im2double(bo))./(-1000);         %Calculate ADC maps for each direction
ADC_101=log(im2double(b2)./im2double(bo))./(-1000);
ADC011=log(im2double(b3)./im2double(bo))./(-1000);
ADC01_1=log(im2double(b4)./im2double(bo))./(-1000);
ADC110=log(im2double(b5)./im2double(bo))./(-1000);
ADC_110=log(im2double(b6)./im2double(bo))./(-1000);

transform=[1 1 0 2 0 0;1 0 1 0 2 0;0 1 1 0 0 2;1 1 0 -2 0 0;1 0 1 0 -2 0;0 1 1 0 0 -2]  %the transformation matrix to obtain tensor elements

xx=zeros(256);      %Create empty tensor elements
yy=zeros(256);
zz=zeros(256);
xy=zeros(256);
xz=zeros(256);
yz=zeros(256);
lambda1=zeros(256);     %create empty images
lambda2=zeros(256);
lambda3=zeros(256);

r=zeros(256);
g=zeros(256);
b=zeros(256);
FAmap=zeros(256);
tADC=zeros(256);

for i=1:256     %Fill in tensor elements pixel by pixel. 
    for j=1:256
        if (bo(i,j)>60) %noise threshold.  Change this to zero if you want the whole image to be calculated (takes longer)
            ADCm=[ADC110(i,j);ADC101(i,j);ADC011(i,j);ADC_110(i,j);ADC_101(i,j);ADC01_1(i,j)]; %see Hedehus, http://www-radiology.stanford.edu/majh/
            ADCe=inv(transform)*ADCm;
            xx(i,j)=ADCe(1,1);
            yy(i,j)=ADCe(2,1);
            zz(i,j)=ADCe(3,1);
            xy(i,j)=ADCe(4,1);
            xz(i,j)=ADCe(5,1);
            yz(i,j)=ADCe(6,1);
            ten=[xx(i,j) xy(i,j) xz(i,j);xy(i,j) yy(i,j) yz(i,j);xz(i,j) yz(i,j) zz(i,j)];  %the tensor itself
            [V,D]=eig(ten);
            D=eig(ten);     %get eigenvalues
            D=abs(D);
            D=sort(D);      %sort the eigenvalues (upwards)
            e1=D(3,1);e1=e1*2;
            e2=D(2,1);e2=e2*2;
            e3=D(1,1);e3=e3*2;
            lambda1(i,j)=e1;
            lambda2(i,j)=e2;
            lambda3(i,j)=e3;
            V=abs(V);       %get the x,y,z components of e3   
            r(i,j)=V(1,3);
            g(i,j)=V(2,3);
            b(i,j)=V(3,3);
            trace=(e1+e2+e3)/3;
            %FA=(sqrt(3/2)*sqrt((1/3)*((e1-e2)^2+(e2-e3)^2+(e3-e1)^2)))/(sqrt(e1^2+e2^2+e3^3)); %Another formula for FA
            FA=(sqrt(3*((e1-trace)^2+(e2-trace)^2+(e3-trace)^2)))/(sqrt(2*(e1^2+e2^2+e3^2)));  %Le Bihan 2001
            tADC(i,j)=trace;
            FAmap(i,j)=FA;
        end;
    end;
end;




rgb=cat(3,r,g,b);
FAmap3=cat(3,FAmap,FAmap,FAmap); %needed to combine with the colormap, need 3D.
cm=rgb.*FAmap3; %FA weighting for colormap
figure, image (cm);axis image






        

