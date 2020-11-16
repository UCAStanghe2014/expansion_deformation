%% By Dr. He Tang at UCAS, Beijing, 2020/11/16.
%% calculate deformation by a supposed nuclear expansion
clc;
clear;
close all;
lambda = 34.2e9;
mu = 26.6e9;
depth = 700; R= 6371d3;
volume = 2e6; % m^3
ff=(lambda+2*mu)/(3*lambda+2*mu)*volume;
ff1=ff/(R^2);
ff2=ff/(R^3);
ff3=9.8*ff/(R^3);
distance = linspace(0,200/6371,800); distance(distance>180)=[]; %degree

openid = fopen('distance.txt','w');
fprintf(openid,'%e\n',distance);
fclose(openid);

ntheta = length(distance);
greenF = zeros(ntheta,6);
greenS = zeros(ntheta,10);
for ii = 1:ntheta
  theta = distance(ii);
  % spherical Earth model
  greenS(ii,1) = theta;
  greenS(ii,2:end) = green_expansion(lambda,mu,depth,R,theta);
end
tang = greenS;
tang(:,1)=tang(:,1)*6371;
tang(:,2)=tang(:,2)*ff1;
tang(:,3)=tang(:,3)*ff1;
tang(:,4)=tang(:,4)*ff2*1e3;
tang(:,5)=tang(:,5)*ff2*1e3;
tang(:,6)=tang(:,6)*ff2*1e3;
tang(:,7)=tang(:,7)*ff1*1e8;
tang(:,8)=tang(:,8)*ff3*1e10;
tang(:,9)=tang(:,9)*ff3*1e6;
tang(:,10)=tang(:,10)*ff2*1e12;
%% plot
xmax = max(tang(:,1))/sqrt(2);
xstep = abs(tang(2,1)-tang(1,1))/2;
[x0,y0] = meshgrid(-xmax:xstep:xmax,-xmax:xstep:xmax);
z = zeros(size(x0,1),size(x0,2),9);
for ii = 1:9
z(:,:,ii) =  interp1(tang(:,1),tang(:,1+ii),x0.*x0+y0.*y0);
end
%
zlables = {'(a) u_r (m)','(b) u_{\theta} (m)','(c) e_{rr} (\times10^{-3})','(d) e_{\theta\theta} (\times10^{-3})','(e) e_{\phi\phi} (\times10^{-3})','(f) N (\times10^{-8} m)','(g) {\Delta}g (\times10^{-10} m/s^2)','(h) {\delta}g (\times10^{-6} m/s^2)','(i) \xi (\times10^{-12}rad)'};
zrange=[-0.5,1;...
    0,0.6;...
    -1.5,0;...
    -0.5,1.5;...
    0,1.5;...
    6.5,7.5;...
    0,1.5;...
    -3,0;...
    -3,0];
%
figure(1); clf
set(gcf,'Color','White')
for ii = 1:9
subplot(3,3,ii);
surf(x0,y0,squeeze(z(:,:,ii)),'edgecolor','none','facecolor','interp')
bar = colorbar;
caxis(zrange(ii,:));
colormap jet
bar.TickLength=0.02;
bar.TickDirection='out';
bar.LineWidth=1.5;
title(zlables(ii))
set(gca,... 
'Box','off',... 
'LooseInset',get(gca,'TightInset')* 1.5,... 
'TickDir','in',... 
'XMinorTick','off',... 
'YMinorTick','off',... 
'TickLength' ,[.02 .02],... 
'LineWidth',1,... 
'XGrid','off',... 
'YGrid','off' ,... 
'FontSize',18);
xlim([-20,20]);
ylim([-20,20]);
end
