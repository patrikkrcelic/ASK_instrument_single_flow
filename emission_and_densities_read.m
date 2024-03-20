function [] = emission_and_densities_read(datadir)

% Function reads in the emissions and densities from the ion chemistry
% model and saves them in .mat files for easier and faster use.
%
% INPUTS:
%    datadir     - Data directory of the event
%



en=importdata(strcat(datadir, '/en_range.txt'));
f1=fopen(strcat(datadir,'/emissions.txt'));
s=fgets(f1);

for i=1:120
    s=fgets(f1);
    t=str2num(s);
    for j=1:201
        s=fgets(f1);
        if i/2==round(i/2)
           alt(i/2,j)=str2num(s(1:20));
           N2_em(i/2,j)=str2num(s(21:40));
           OI_em(i/2,j)=str2num(s(41:60));
           OI_D_em(i/2,j)=str2num(s(61:80));
       end
    end
end

fclose(f1);

f2=fopen(strcat(datadir,'/densities.txt'));

for i=1:120
    s=fgets(f2);
    s=fgets(f2);
    for j=1:201
        s=fgets(f2);
        if i/2==round(i/2)
           N2_den(i/2,j)=str2num(s(31:60));
           O2_den(i/2,j)=str2num(s(61:90));
           O_den(i/2,j)=str2num(s(91:120));
           e_den(i/2,j)=str2num(s(121:150));
           Te(i/2,j)=str2num(s(151:end));
           
        end
    end
    s=fgets(f2);
    for j=1:201
        s=fgets(f2);
        if i/2==round(i/2)
           prod(i/2,j)=str2num(s);
        end
    end
end

fclose(f2);

save(strcat(datadir,'/emissions.mat'), 'alt', 'en', 'OI_D_em', 'OI_em', ...
    'N2_em', 'prod');
save(strcat(datadir,'/densities.mat'), 'alt','N2_den', 'O_den', 'O2_den',...
    'e_den', 'Te', 'en'); 

end






