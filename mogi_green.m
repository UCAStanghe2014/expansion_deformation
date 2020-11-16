function green = mogi_green(la,mu,depth,R,theta)
% By Dr. He Tang at UCAS, Beijing, 2020/11/16.
% displacement,strain in Cylindrical system (r,theta,z)
r = R*(theta/180*pi); % �����ϵľ���
d = depth;
ratio=R^2*((3*la+2*mu)/(la+mu)/(2*pi));
% displacement
ur =  r/(sqrt(d^2+r^2))^3;
uz = -d/(sqrt(d^2+r^2))^3;
% strain
er2 = R*(d^2-2*r^2)/(sqrt(d^2+r^2))^5; %% ˮƽ������
etheta2 = R/(sqrt(d^2+r^2))^3;         %% ˮƽ������
ez2 = R*(la*(r^2+d^2-3*d)+3*(d-1)*d*mu)/(sqrt(d^2+r^2))^5/(la+2*mu);  %% ��ֱ���� ���ܴ�Ľ���!
% ez2 = R*(-2*d^2+r^2)/(sqrt(d^2+r^2))^5*(la/(la+2*mu)); %��ȷ��ֵ
green(1)=ur*ratio;      % UdS/R^2 is omitted
green(2)=uz*ratio;      % UdS/R^2 is omitted
green(3)=er2*ratio;     % UdS/R^3 is omitted
green(4)=etheta2*ratio; % UdS/R^3 is omitted
green(5)=ez2*ratio;     % UdS/R^3 is omitted
end