function green=green_expansion(la,mu,depth,R,theta)
% By Dr. He Tang at UCAS, Beijing, 2020/11/16.
ep = (R-depth)/R;
ratio = la/(la+2*mu);

HLK = hlk_expansion(la,mu,depth,R);
[A,B,C] = ABC_sum(theta,ep);

A030=A(1);
A020=A(2);
A010=A(3);
A000=A(4);
A0m0=A(5);

% B010=B(1);
B000=B(2);
B0m0=B(3);

% C010=C(1);
C000=C(2);
C0m0=C(3);

h1=HLK(1,1);
h0=HLK(1,2);
hm=HLK(1,3);

l0=HLK(2,2);
lm=HLK(2,3);

k0=HLK(3,2);
km=HLK(3,3);

 green(1) = h1*A010 + h0*A000 + hm*A0m0;      %u_r
 green(2) =           l0*B000 + lm*B0m0;      %u_theta

 green(3) = l0*A020 + (lm+l0-2*h1)*A010 + ...
            (lm-2*h0)*A000 - 2*hm*A0m0;        %e_rr
 green(3) = ratio*green(3);
 green(4) = h1*A010 + h0*A000 + hm*A0m0 + ...
            l0*C000 + lm*C0m0;                 %e_thetatheta
 green(5) = h1*A010 + h0*A000 + hm*A0m0 + ...
            (l0*B000 + lm*B0m0)*cotd(theta);   % e_phiphi

 green(6) = k0*A000 + km*A0m0;                 % N geoid
 green(7) = k0*A010 + (k0+km)*A000 + km*A0m0;  % gravity change at space-fixed point
 green(8) = (k0-2*h1)*A010 + (km+k0-2*h0)*A000 + (km-2*hm)*A0m0;
                                               % gravity change tat deformed surface
 green(9) = k0*B000 + km*B0m0;                 % deflection or local tilt change

end 