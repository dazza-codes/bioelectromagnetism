function h = montage(avw,cmap)
%


fprintf('in development, sorry'); return

if ~exist('cmap','var'), cmap = 'gray'; end
if isempty('cmap'), cmap = 'gray'; end


b = a(1,1); % to inherit type 
b(1,1) = 0; % from a
b = repmat(b, [size(a,1)*size(a,3), size(a,2)*size(a,3), 1]);

for z =1:size(a,3),
    for i=1:10,
        b(rows+i*siz(1),cols+j*siz(2),1) = a(:,:,z);
    end
end

h = imshow(b,cm);
