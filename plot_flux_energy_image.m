function plot_flux_energy_image(datdir, imnum, type)

% Function plots precitipating electron flux or energy tomographic image. 
%
% INPUTS:
%   datadir - data directory
%   imnum   - image number
%   type    - set ot 1 or 2. 1 plots precipitating electron flux and 2
%             plots precipitating electron energy
%

sp=num2str(imnum,'%05.f');

imfile=strcat(datdir, '/results/energy_and_flux/model_', sp, '.mat');

if type == 1
    load(imfile, 'Mflux');
    M = Mflux;
elseif type == 2
    load(imfile, 'Men');
    M = Men;
else
    disp('Wrong type input!')
    return
end

M2 = zeros(150,150);
MF = zeros(150,150);

if length(size(M))==3
    for i=1:150
        for j=1:150
            M2(i,j)=median(M(i,j,M(i,j,:)>0));
        end
    end
end

M2(isnan(M2)) = 0;
M2 = M2';

for i=1:150
    MF(i,:) = M2(150-i+1,:); 
end

if type == 2 
    MF = log(MF);
end

h=pcolor(MF); 
axis equal;
set(h, 'EdgeColor', 'none');
colorbar;
xlim([1,150])
ylim([1,150])


end 