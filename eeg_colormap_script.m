
% eeg_colormap_script - play with colour mapping

clear all, close all

p = eeg_toolbox_defaults;
p.colorMap.style = 'Red/Blue/White';

p = eeg_colormap(p);

colormap(p.colorMap.map);

for i = 1:size(p.colorMap.map,1)

    c = p.colorMap.map(i,:);   % define patch color

    xp = [-5  5  5 -5];
    yp = [ i  i  i+1 i+1];
    
    patch(xp,yp,c)
end

axis tight

colorbar
