function [El, STD_El] = plot_el_field(velocities, STD, B, theta) 





STD(isnan(STD(:,1)),1)=mean(STD(~isnan(STD(:,1)),1));
STD(isnan(STD(:,2)),2)=mean(STD(~isnan(STD(:,2)),2));

J=jet(length(velocities(5:end-2,1)));

n=size(velocities);
ft=[1 3 5 3 1]/13;
VV = zeros(n(1)-4,2);
S = zeros(n(1)-4, 2);

for i=3:n(1)-2
    VV(i-2,1)=ft*velocities(i-2:i+2,1);
    VV(i-2,2)=ft*velocities(i-2:i+2,2);
    S(i-2,1)=ft*STD(i-2:i+2,1);
    S(i-2,2)=ft*STD(i-2:i+2,2);
end

velocities = VV;
STD = S;

if nargin == 4
    Rot=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];
else
    Rot=[cosd(-1.651) -sind(-1.651); sind(-1.651) cosd(-1.651)];
end

if nargin ~= 3
   B=48916*1e-9;
end

subplot(1,2,1)

for i=5:length(velocities(:,1))

    vv=velocities(i,:);
    vv=vv*Rot;
    velocities(i,:)=vv;
    
    hold on;
    plot([0,velocities(i,1)], [0, velocities(i,2)],'color', J(i-4,:));
    axis equal;
    
end
 

xlim([-4, 4]);
ylim([-4, 4]);
xlabel('v_E [km/s]');
ylabel('v_N [km/s]');
hc = colorbar;
clim([0.5, (n(1)-4)/10]);
ylabel(hc, 'time steps');


h=plot([0, mean(velocities(1:end,1))], [0, mean(velocities(1:end,2))], 'k');
set(h,'LineWidth',2);

% electric fields

El = zeros(n(1)-4, 2);
STD_El = zeros(n(1)-4, 2);

subplot(1,2,2)
for i = 5 : length(velocities(:,1))
    vv=velocities(i,:);
    E = -cross(1000*[vv(1), vv(2), 0], B*[0, 0, -1]); % V/m
    El(i,:) = -E(1:2)*1000;  % mV/m

    STD_El(i,1) = sqrt(B^2*(STD(i,2)*1000)^2) *1000;
    STD_El(i,2) = sqrt(B^2*(STD(i,1)*1000)^2) *1000;
              
    plot([0,El(i,1)], [0, El(i,2)],'color', J(i-4,:));
    axis equal;
    hold on;

end

h=plot([0, mean(El(:,1))], [0, mean(El(:,2))], 'k');
set(h,'LineWidth',2);

hc = colorbar;
colormap jet
xlim([-200, 200]);
ylim([-200, 200]);
xlabel('E_E [mV/m]');
ylabel('E_N [mV/m]');
clim([0.5, (n(1)-4)/10]);
ylabel(hc, 't [s]');
 
 
end 
 
 
 
 
 
 
 
 
 

