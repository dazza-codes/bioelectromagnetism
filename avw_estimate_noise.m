function noise_mean = avw_estimate_noise(avw)



[s1,s2,s3] = size(avw.img);

slice_x0 = avw.img(1,:,:);
slice_x0_mean = mean(mean(slice_x0));

if s1 > 1,
    slice_x1 = avw.img(end,:,:);
    slice_x1_mean = mean(mean(slice_x1));
    x_mean = mean(slice_x0_mean, slice_x1_mean);
else
    x_mean = slice_x0_mean;
end
    
slice_y0 = avw.img(:,1,:);
slice_y0_mean = mean(mean(slice_y0));

if s2 > 1,
    slice_y1 = avw.img(:,end,:);
    slice_y1_mean = mean(mean(slice_y1));
    y_mean = mean(slice_y0_mean, slice_y1_mean);
else
    y_mean = slice_y0_mean;
end

slice_z0 = avw.img(:,:,1);
slice_z0_mean = mean(mean(slice_z0));

if s3 > 1,
    slice_z1 = avw.img(:,:,end);
    slice_z1_mean = mean(mean(slice_z1));
    z_mean = mean(slice_z0_mean, slice_z1_mean);
else
    z_mean = slice_z0_mean;
end

noise_mean = mean([x_mean,y_mean,z_mean]);

