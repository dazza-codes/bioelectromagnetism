function[Image] = RotateImage(Image, degrees)
% RotateImage - Rotates an image by X degrees
% 
%  Image_Rotated = RotateImage(Image, degrees)
%
%  Rotates an image by the angle degrees in the 
%  CCW direction.  Degrees may be any number.  
%  The function will put degrees in the range 0 
%  to 360 degrees and then into a range of -45 to 45
%  degrees after performing elementary 90 degree rotations.
%
%  The rotation performed is the 3 Pass Separable rotation
%  using FFT-based methods to perform the skews.  Aliased spectral 
%  components are masked after EACH skew using the MaskImage function.
%
%  Normal Rotation
%
%  |x'| = | cos(a) -sin(a) ||x|
%  |y'|   | sin(a)  cos(a) ||y|
%
%  3 Pass Rotation
%
%  |x'| = | 1  -tan(a/2) || 1       0 || 1  -tan(a/2) ||x|
%  |y'|   | 0     1      || sin(a)  1 || 0     1      ||y|
%
%  This function pads images during the rotation process.  Images 
%  can be non-square, but odd dimensions will cause incorrect padding
%  Therefore an error is returned, if a dimension is odd.
%
%  Dimensions need not be powers of 2 and will not necessarily be padded
%  to a power of 2.  The Matlab fft functions are capable of non-power of 2
%  transforms.
%
%  Note: The results often display some Gibbs ringing.
%
%  Example:
%
%  load mri
%  img = double( squeeze(D(:,:,1,16)) );
%  img_rot = real( RotateImage(img,45) );
%  figure;
%  subplot(1,2,1); imagesc(img,[0 88]); axis image; 
%  subplot(1,2,2); imagesc(img_rot,[0 88]); axis image;

%
%************************************************************************
%*                    (c) Copyright 1998				                *
%*                  Biomathematics Resource 			       	        *
%*                      Mayo Foundation                              	*
%************************************************************************
%
% 11/22/98  Implemented by Edward Brian Welch, edwardbrianwelch@yahoo.com
%


% NOTE:
% This function is partially vectorized (only single loops, 
% no double loops) to gain some benefit in speed.  The tradeoff 
% is that more memory is required to hold the large matrices.
% It would be possible to completely vectorize certain parts 
% of the mfile, but that would require large amounts of memory
% and would also make the code less readable.

[xdim ydim] = size(Image);

if mod(xdim,2)==1 | mod(ydim,2)==1,
   error('RotateImage() can only rotate images with even dimensions.');
end

% Determine type of image to return
if isreal(Image),
   real_flag = 1;
else
   real_flag = 0;
end

% Put degrees into the range [0,360]
while degrees<0,
   degrees = degrees + 360;
end

while degrees>360,
   degrees = degrees - 360;
end

% Count number of elementary 90 degree rotations
% to put rotation into range [-45,45]
count=0;
while degrees>45,
   degrees = degrees - 90;
   count = count + 1;
end

Image = rot90(Image, count);

% Calculate Trigonometric Values
% Not Necessary in Matlab implementation to Negate degrees 
% so that CCW rotations are positive
theta = (degrees/360)*2*pi;
sin_theta = sin(theta);
negtan_thetadiv2 = -tan(theta/2);

%---------------------------
% FIRST SKEW
% |x1| = | 1  -tan(a/2) ||x|
% |y1| = | 0     1      ||y|
%---------------------------

% Pad image rows to double the row size
Image2 = zeros(2*xdim,ydim);
Image2( (xdim/2+1):(xdim/2+xdim),:)=Image;

% Forward FFT image rows
Image2 = fft(Image2, 2*xdim, 1);

% Calculate image's center coordinates
xno = (2*xdim-1)/2;
yno = (ydim-1)/2;

% Prepare to use the Fourier Shift Theorem
% f(x-deltax) <=> exp(-j*2*pi*deltax*k/N)*F(k)

% Initialize constant part of the exponent
% expression.  The skew is row dependent.
cons1 = (-2.0*pi/(2*xdim))*negtan_thetadiv2;

% Calculate k values (Nyquist is at x=xno)
k_array = zeros((2*xdim),1);

for x=1:(2*xdim),    
      if (x-1)<=xno,
         k_array(x) = (x-1);
      else
         k_array(x) = (x-1-2*xdim);
      end
end

for y=1:ydim,
   % Skew dependent part of expression
	cons2 = cons1 * (y-1-yno);
   
   % Calculate the angles
   angle_array = cons2*k_array;
   
   % Rotate the complex numbers by those angles
   sin_ang = sin(angle_array);
   cos_ang = cos(angle_array);
   newr = real(Image2(:,y)).*cos_ang - imag(Image2(:,y)).*sin_ang;
   newi = real(Image2(:,y)).*sin_ang + imag(Image2(:,y)).*cos_ang;
   Image2(:,y) = newr + newi*i;
end

%---------------------------
% SECOND SKEW
% |x2| = | 1       0 ||x1|
% |y2| = | sin(a)  1 ||y1|
%---------------------------

% Pad the image columns in hybrid space
Image22 = zeros(2*xdim, 2*ydim);
Image22(:,(ydim/2+1):(ydim/2+ydim) )=Image2;

% Perform a Forward FFT on the image columns
Image22 = fft(Image22, 2*ydim, 2);

% Mask aliased components from SKEW 1
Image22 = fftshift(Image22);
Image22 = MaskImage( Image22, 1, 0, negtan_thetadiv2, 1);
Image22 = fftshift(Image22);

%  Inverse FFT the image rows
Image22 = ifft(Image22, 2*xdim, 1);

% Calculate image's center coordinates
xno = (2*xdim-1)/2;
yno = (2*ydim-1)/2;

% Initialize constant part of the exponent
% expression.  The skew is column dependent.
cons1 = (-2.0*pi/(2*ydim))*sin_theta;

% Calculate k values (Nyquist is at y=yno)
k_array = zeros(1,(2*ydim));
for y=1:(2*ydim),
      if (y-1)<=yno,
      k_array(y) = (y-1);
   else
      k_array(y) = (y-1-2*ydim);
   end
end

for x=1:(2*xdim),  
   % Skew dependent part of expression
   cons2 = cons1 * (x-1-xno);
   
   % Calculate the angles
   angle_array = cons2*k_array;
   
   % Rotate the complex numbers by those angles
   sin_ang = sin(angle_array);
   cos_ang = cos(angle_array);
   newr = real(Image22(x,:)).*cos_ang - imag(Image22(x,:)).*sin_ang;
   newi = real(Image22(x,:)).*sin_ang + imag(Image22(x,:)).*cos_ang;
   Image22(x,:) = newr + newi*i;
end

%---------------------------
% THIRD SKEW
% |x3| = | 1  -tan(a/2) ||x2|
% |y3| = | 0     1      ||y2|
%---------------------------

% Forward FFT image rows
Image22 = fft(Image22, 2*xdim, 1);

% Mask aliased components from SKEW 2
Image22 = fftshift(Image22);
Image22 = MaskImage( Image22, 1, sin_theta, 0,  1);
Image22 = fftshift(Image22);

% Inverse FFT the image columns
Image22 = ifft(Image22, 2*ydim, 2);

% Crop image columns
Image2 = Image22(:,(ydim/2+1):(ydim/2+ydim) );

% Calculate image's center coordinates
xno = (2*xdim-1)/2;
yno = (ydim-1)/2;

% Initialize constant part of the exponent
% expression.  The skew is row dependent.
cons1 = (-2.0*pi/(2*xdim))*negtan_thetadiv2;

% Calculate k values (Nyquist is at x=xno)
k_array = zeros((2*xdim),1);

for x=1:(2*xdim),
      
      if (x-1)<=xno,
         k_array(x) = (x-1);
      else
         k_array(x) = (x-1-2*xdim);
      end
end

for y=1:ydim,
   % Skew dependent part of expression
	cons2 = cons1 * (y-1-yno);
   
   % Calculate the angles
   angle_array = cons2*k_array;
   
   % Rotate the complex numbers by those angles
   sin_ang = sin(angle_array);
   cos_ang = cos(angle_array);
   newr = real(Image2(:,y)).*cos_ang - imag(Image2(:,y)).*sin_ang;
   newi = real(Image2(:,y)).*sin_ang + imag(Image2(:,y)).*cos_ang;
   Image2(:,y) = newr + newi*i;
end

% Forward FFT the image columns
Image2 = fft(Image2, ydim, 2);

% Mask aliased components from SKEW 3
% The /2.0 factor is there because columns have been cropped in half
Image2 = fftshift(Image2);
Image2 = MaskImage( Image2, 1, 0, negtan_thetadiv2/2.0, 1);
Image2 = fftshift(Image2);

% Inverse 2D FFT the image
Image2 = ifft2(Image2);

% Crop the rows to obtain final result
Image = Image2( (xdim/2+1):(xdim/2+xdim),:) ;

% Return a Real image if original Image was Real
if real_flag==1,
   Image = real(Image);
end

%
%  Image_Masked = MaskImage(Image, a, b, c, d)
%
%  Masks an image according to the (x,y) coordinate
%  transform matrix:  |x'| = | a  b | * |x|
%                     |y'|   | c  d |   |y|
%
%  (x,y) positions that become (x',y') outside the boundary
%  of the original image will be set to 0.
%
%  The masking is softened at the edges by a convolution.
%

%
%************************************************************************
%*                    (c) Copyright 1998				                *
%*                  Biomathematics Resource 			       	        *
%*                      Mayo Foundation                              	*
%************************************************************************
%
% 11/22/98  Implemented by Edward Brian Welch, edwardbrianwelch@yahoo.com
%
function[Image] = MaskImage(Image, a, b, c, d)

% NOTE:
% This function is highly vectorized (no loops)
% to gain some benefit in speed.  The tradeoff is that 
% more memory is required to hold the large matrices

[xdim ydim] = size(Image);

% Calculate center of the image
xno = (xdim-1)/2;
yno = (ydim-1)/2;

x_values = 1:1:xdim;
x_values = x_values - 1 - xno;
x_mat = repmat(x_values',1,ydim);

y_values = 1:1:ydim;
y_values = y_values - 1 - yno;
y_mat = repmat(y_values,xdim,1);

% Calculate new x values for this column    
new_x_mat = a*x_mat + b*y_mat;

% Calculate new y values for this column
new_y_mat = c*x_mat + d*y_mat;

% Points whos x or y position is greater than
% absolute value of xno or yno respectively
% should be masked (set to zero).
new_x_mat = abs(new_x_mat) - xno;
new_y_mat = abs(new_y_mat) - yno;

% Only valid points will have a non-zero value after 
% the operation below
new_x_mat = sign(new_x_mat)-1;
new_y_mat = sign(new_y_mat)-1;

% ANDing the matrices shows good points
mask = new_x_mat & new_y_mat;

% Soften mask edges
Nx=ceil(xdim*.2);
Ny=ceil(ydim*.2);
mask = conv2(mask,ones(Nx,Ny)/(Nx*Ny),'same');

Image = Image.*mask;
  

        
