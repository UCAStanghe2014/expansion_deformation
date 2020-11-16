%% By Dr. He Tang at UCAS,Beijing,2020.
%% Calculate Green's functions by our asymptotic expressons
lambda = 34.2e9; mu = 26.6e9; %PREM model
depth = 32; R= 6371;
f1=1e3;
f2=1e6;

distance = 10.^(linspace(-3,2.5,40))'; distance(distance>180)=[];
openid = fopen('distance.txt','w');
fprintf(openid,'%e\n',distance);
fclose(openid);
ntheta = length(distance);
greenF = zeros(ntheta,6);
greenS = zeros(ntheta,10);
for ii = 1:ntheta
  theta = distance(ii);  
  % half-space model
  deforF = mogi_green(lambda,mu,depth,R,theta);
  greenF(ii,1) = theta;
  greenF(ii,2:end) = deforF;
  % spherical Earth model
  greenS(ii,1) = theta;
  greenS(ii,2:end) = green_expansion(lambda,mu,depth,R,theta);
end
% Okada's results code available at https://www.bosai.go.jp/e/dc3d.html
okada=load('green_okada.txt');

%% displacement and strain
mogi_u_up  = -greenF(:,3)/f1;
mogi_u_h   =  greenF(:,2)/f1;
mogi_area  =  greenF(:,4)/f2+greenF(:,5)/f2;
mogi_e_upup= -greenF(:,6)/f2;
mogig = [mogi_u_up,mogi_u_h,mogi_area,mogi_e_upup];

tang_u_up  = greenS(:,2)/f1;
tang_u_h   = greenS(:,3)/f1;
tang_area  = greenS(:,5)/f2+greenS(:,6)/f2;
tang_e_upup= greenS(:,4)/f2;
tangg = [tang_u_up,tang_u_h,tang_area,tang_e_upup];

okada_u_up  = okada(:,4);
okada_u_h   = okada(:,2);
okada_area  = okada(:,5)+okada(:,9);
okada_e_upup= okada(:,13);
okadag = [okada_u_up,okada_u_h,okada_area,okada_e_upup];

%% Plot
zlables = {'(a) u_r','(b) u_{\theta}','(c) e_{area}','(d) e_{rr}'};

figure(1); clf
set(gcf,'Color','White')
for ii = 1:4
subplot(2,2,ii)
semilogx(distance,tangg(:,ii),'r',...
         distance,mogig(:,ii),'kx',...
         okada(:,1),okadag(:,ii),'bo',...
         'MarkerSize',4,...
         'LineWidth',1);
title(zlables(ii))
xlim([1e-3,180])
set(gca,'xtick',10.^(-3:2))
set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'tickdir','out');
grid on
if(ii==4)
    legend('This study','Mogi (1958)','Okada (1985)','Location','southeast')
end
if (ii>2)
    xlabel('Angular distance \theta (\circ)');
end
end

figure(2)
set(gcf,'Color','White')
semilogx(distance,-tangg(:,4)./(tangg(:,3)*(lambda/(lambda+2*mu))),'r',...
         okada(:,1),-mogig(:,4)./(mogig(:,3)*(lambda/(lambda+2*mu))),'bo',...
         'MarkerSize',4,...
         'LineWidth',1);
xlim([1e-3,180])
set(gca,'xtick',10.^(-3:2))
set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'tickdir','out');
grid on
ylim([-3,4])
legend('This study','Mogi (1985)','Location','southeast')
xlabel('Angular distance \theta (\circ)');
ylabel('Ratio of direct and indirect vertical linear strain')
