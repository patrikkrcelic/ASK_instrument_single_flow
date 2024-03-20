function plot_ASK_image(datdir, cam, imnum, bck)

% Function plots ASK image 
%
% INPUTS:
%   datadir - data directory
%   cam     - camera number
%   imnum   - image number
%   bck     - background substraction number
%

sp=num2str(imnum,'%05.f');

imfile=strcat(datdir,'/ASK', num2str(cam), '/ask', num2str(cam), '_', sp,'.txt');
A1=importdata(imfile);

A1 = image_trimming(A1, bck);

h=pcolor(A1); 
axis equal;
set(h, 'EdgeColor', 'none');
colorbar;
xlim([1,256])
ylim([1,256])


end 