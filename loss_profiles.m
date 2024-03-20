function [] = loss_profiles(datadir, EN)

% Function caluclates the "time of life" of the O=(2P) J=1/2 and J=3/2
% states for the average energy of the event and saves them in the file: 
%  'datadir+/loss_profiles.mat'
%
% INPUTS:
%    datadir - Data directory of the event
%    EN      - Mean energy of the event
%


load(strcat(datadir, '/densities.mat'));
I=find(en>=EN);
if length(I)<1
    I=length(en);
end

tau=zeros(201,1);
tau_1_2=zeros(201,1);
tau_3_2=zeros(201,1);

for i=1:201
    tau(i)=e_den(I(1),i)*10^(-6)*4*10^(-8)*(Te(I(1),i)/300)^(-0.5);
    tau(i)=tau(i)+e_den(I(1),i)*10^(-6)*1.5*10^(-7)*(Te(I(1),i)/300)^(-0.5);
    tau(i)=tau(i)+O_den(I(1),i)*10^(-6)*5*10^(-11);
    tau(i)=tau(i)+N2_den(I(1),i)*10^(-6)*1.8*10^(-10);
    
    tau_1_2(i)=tau(i)+5.63*10^(-2)+9.39*10^(-2)+2.34*10^(-2);
    tau_3_2(i)=tau(i)+1.07*10^(-1)+5.78*10^(-2)+5.7*10^(-2);
    
end

tau_1_2=1./tau_1_2;
tau_3_2=1./tau_3_2;

save(strcat(datadir,'/loss_profiles.mat'), 'tau_1_2', 'tau_3_2');

end