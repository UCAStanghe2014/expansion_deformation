function [A,B,C]=ABC_sum(theta,ep)
% By Dr. He Tang at UCAS, Beijing, 2020/11/16.
% ep = epsilon = rs/a
c = cosd(theta);
s = sind(theta);
w = sqrt(1+ep^2-2*ep*c);

A030 = (-ep^4*(c+ep)+2*ep^3*(5-c^2)+c*ep^2*(c^2-9)+ep*(5*c^2-4)+c)/(w^7);

A020 = ((ep-c)*ep^2+(c^2-2)*ep+c)/(w^5);
A010 = (c-ep)/(w^3);
A000 = (1-w)/(ep*w);
A0m0 = log(2/(1-ep*c+w))/ep;

B010 = -s*(1+ep*c-2*ep^2)/(w^5);
B000 = -s/(w^3);
B0m0 = -s*(1+w)/(w*(1-ep*c+w));

C010 = -((9-4*c^2)*ep^3+c*(c^2-4)*ep^2+(5*c^2-6)*ep-2*c*ep^4+c)/(w^7);
C000 = (3*ep*s^2-c*w^2)/(w^5);
C0m0 = (1/4)*w^(-3)*(1+w-c*ep)^(-2)*(-c*(8*(1+w)+(35+16*w)*ep^2+4*ep^4)+...
       ep*(2*(10+8*w+(7+2*w)*ep^2)+(8+4*w+6*ep^2)*cosd(2*theta)-ep*cosd(3*theta)));

A = [A030;A020;A010;A000;A0m0];
B = [B010;B000;B0m0];
C = [C010;C000;C0m0]; 
   

end