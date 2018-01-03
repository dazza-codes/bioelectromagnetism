function [m] = pifft(pksp)

% PIFFT - Homodyne reconstruciton of partial MRI k-space (Spatial Frequency) data
%
%  IMG = PIFFT(KSP)
%  The input KSP is assumed to be a rectangular matrix with the number of rows
%  equaling the number of phase encodes (lines) of the image's k-space.  The 
%  first rows of the matrix should be the highest frequency phase encodes 
%  acquired so that the center of k-space is near the bottom of the matrix.
%  The image is reconstructed for the next highest power of 2 greater than the
%  number of phase encodes in the partial k-space data set.
%
%  Example:
%
%  load mri
%  img = double( squeeze(D(:,:,1,16)) );
%  ksp = fftshift( fft2( img ) );
%  pksp = ksp(1:68,:);
%  m=pifft(pksp);
%
%  figure;
%  subplot(1,3,1); imagesc(img); axis image; title('128 of 128 of k-space lines');
%  subplot(1,3,2); imagesc(m); axis image; title('68 of 128 of k-space lines');
%  subplot(1,3,3); imagesc(abs(img-m)); axis image; title('abs. difference');
%  colormap(gray);
%
%  Reference:
%
%  Noll DC, Nishimura DG, Macovski A. Homodyne Detection in Magnetic 
%  Resonance Imaging.  IEEE Trans. on Medical Imaging 1991; 10(2):154-163
%
%

%
%    Written by Edward Brian Welch (edwardbrianwelch@yahoo.com)
%    MRI Research Lab, Mayo Graduate School, November 2001
%

% ERROR CHECKING
if ndims(pksp)~=2 | ~isnumeric(pksp),
    error('PFFT operates on two-dimensional numerical data');
end

% Convert data into hybrid space by inverse FFT along FE dir.
% for better numerical accuracy.
pksp = ifft(pksp,[],2);

% Frequency Encode direction is assumed to be the longer direction
[Npp Nf] = size(pksp);

% Number of phase encodes in the full k-space
Np = 2^nextpow2(Npp+1);

% Number of high frequency phase encodes
NH = Np-Npp;

% NUmber of low frequency phase encodes
NL = Npp - NH;

% Row indices of high and low frequency lines
HFind = 1:NH;
LFind = (NH+1):(NH+NL);

% Create weighting window
%
% It multiplies all un-partnered high frequency lines by 2.0
% Low frequency lines are multiplied by a ramp that approaches 
% zero at the zero-padded edge of k-space.  The weighting
% factors of two partnered low frequency lines should have a 
% sum of 2.0
w = zeros(Np,1);
w(HFind) = 2;
rstep = 2/(NL+1);
w(LFind) = [(2-rstep):-rstep:rstep];

% Create weighted partial k-space
HFksp = zeros(Np,Nf);
HFksp(1:Npp,:)=pksp;
HFksp = HFksp.*repmat(w,1,Nf);

% Create Low Frequency k-space 
LFksp = zeros(Np,Nf);
LFksp(LFind,:) = pksp(LFind,:);

% Low frequency image
Rc = ifft(LFksp,[],1);

% Unsynchronous image
Ic = ifft(HFksp,[],1);

% Synchronous image
Is = Ic.*exp(-i*angle(Rc));

% Demodulated image
m = real(Is);

return
