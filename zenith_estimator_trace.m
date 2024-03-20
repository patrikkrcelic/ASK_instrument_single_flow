function Xzen = zenith_estimator_trace(datdir, event_setup, imnum1, imnum2, imnum3, cam)

% Function estimates the magentic zenith position by modeling the 3 ASK
% images labeled with numbeers imnum1, imnum2, imnum3. Function also plots 
% errors of ray tracing for visual inspection. 
%
%
% INPUTS:
%    datdir     - data directory of the event
%    event_seup - matadata of the event
%    imnum1     - number of the image 1
%    imnum2     - number of the image 2
%    imnum3     - number of the image 3
%    cam        - camera number
%
% OUTPUT:
%    Xzen - magnetic zenith position    
%
%


% image num 1

sp=num2str(imnum1,'%05.f');

imfile=strcat(datdir,'/ASK', num2str(cam), '/ask', num2str(cam), '_', sp,'.txt');
A1=importdata(imfile);

A1 = image_trimming(A1, event_setup.bck1);

h=pcolor(cumsum(ones(256,1)),cumsum(ones(256,1)),A1); 
axis equal;
set(h, 'EdgeColor', 'none');
colorbar;
[t, u]=ginput(2);

dx=t(2)-t(1);
dy=u(2)-u(1);

aa=dy/dx;
bb=u(1)-dy/dx*t(1);

a(1)=aa*dx;
b(1)=-dx;
c(1)=bb*dx;


% image num 2

sp=num2str(imnum2,'%05.f');

imfile=strcat(datdir,'/ASK', num2str(cam), '/ask', num2str(cam), '_', sp,'.txt');
A1=importdata(imfile);

A1 = image_trimming(A1, event_setup.bck1);

h=pcolor(cumsum(ones(256,1)),cumsum(ones(256,1)),A1); 
axis equal;
set(h, 'EdgeColor', 'none');
colorbar;
[t, u]=ginput(2);

dx=t(2)-t(1);
dy=u(2)-u(1);

aa=dy/dx;
bb=u(1)-dy/dx*t(1);

a(2)=aa*dx;
b(2)=-dx;
c(2)=bb*dx;


% image num 3

sp=num2str(imnum3,'%05.f');

imfile=strcat(datdir,'/ASK', num2str(cam), '/ask', num2str(cam), '_', sp,'.txt');
A1=importdata(imfile);

A1 = image_trimming(A1, event_setup.bck1);

h=pcolor(cumsum(ones(256,1)),cumsum(ones(256,1)),A1); 
axis equal;
set(h, 'EdgeColor', 'none');
colorbar;
[t, u]=ginput(2);

dx=t(2)-t(1);
dy=u(2)-u(1);

aa=dy/dx;
bb=u(1)-dy/dx*t(1);

a(3)=aa*dx;
b(3)=-dx;
c(3)=bb*dx;



A(1,1)=a*a';
A(1,2)=a*b';
A(2,1)=A(1,2);
A(2,2)=b*b';
B(1,1)=c*a';
B(2,1)=c*b';

Xzen=-inv(A)*B;

e = zeros(256, 256);

for i=1:256
    for j=1:256
        for k=1:3
            e(i,j)=e(i,j)+(a(k)*j+b(k)*i+c(k))^2;
        end
    end
end
        
er=min(e);       
J=find(er==min(er));
I=find(e(:,J)==min(e(:,J)));

subplot(2,2,1)
h=pcolor(cumsum(ones(256,1)),cumsum(ones(256,1)),AA1); 
axis equal; 
xlim([0 256]);
ylim([0 256]);
set(h, 'EdgeColor', 'none');
hold on;
scatter(X(1),X(2),20,'k','*');
scatter(J,I,20,'r','*');
plot([0 256],[-a(1)/b(1)*0-c(1)/b(1) -a(1)/b(1)*256-c(1)/b(1)],'k');


subplot(2,2,2)
h=pcolor(cumsum(ones(256,1)),cumsum(ones(256,1)),AA2); 
axis equal;
xlim([0 256]);
ylim([0 256]);
set(h, 'EdgeColor', 'none');
hold on;
scatter(X(1),X(2),20,'k','*');
scatter(J,I,20,'r','*');
plot([0 256],[-a(2)/b(2)*0-c(2)/b(2) -a(2)/b(2)*256-c(2)/b(2)],'k');
c1=colorbar;
ylabel(c1, 'Re', 'Rotation',-90)

subplot(2,2,3)
h=pcolor(cumsum(ones(256,1)),cumsum(ones(256,1)),AA3); 
axis equal;
xlim([0 256]);
ylim([0 256]);
set(h, 'EdgeColor', 'none');
hold on;
scatter(X(1),X(2),20,'k','*');
scatter(J,I,20,'r','*');
plot([0 256],[-a(3)/b(3)*0-c(3)/b(3) -a(3)/b(3)*256-c(3)/b(3)],'k');

er=e/min(min(e));
er=er*4-4;

subplot(2,2,4)
h=pcolor(cumsum(ones(256,1)),cumsum(ones(256,1)),er); 
axis equal; 
xlim([0 256]);
ylim([0 256]);
set(h, 'EdgeColor', 'none');
c1=colorbar;
ylabel(c1, '\chi^2/\chi^2_{min}', 'Rotation',-90)

end


