clear;clc

Parameters
r  = [0.71,  0, 0.71] * AU;
v  = [0   , -41.34, 0];

df = 0.01;
f1 = 2.68;
f2 = f1 + df;
mu = sun.mu;
a  = norm(r)/(2 - norm(r)*dot(v,v)/mu);
e  = norm(((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu);

diff = @(n) f2 - f1 - ...
       acos((1+e*cos(f2))/sqrt(1+2*e*cos(f2)+e^2)) + ...
       acos((1+e*cos(f1))/sqrt(1+2*e*cos(f1)+e^2)) - ...
       atan(n/sqrt(mu*(1+2*e*cos(f1)+e^2)/a/(1-e^2))) - ...
       atan(n/sqrt(mu*(1+2*e*cos(f2)+e^2)/a/(1-e^2)));
   
noise_lim = fzero(diff,3)