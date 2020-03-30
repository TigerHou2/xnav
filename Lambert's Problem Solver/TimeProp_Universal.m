function [r1,v1] = TimeProp_Universal(r0,v0,mu,dt)
%TIMEPROP_UNIVERSAL Propagates an orbit over time.
%   Given the initial position and velocity vectors of an orbiting
%   spacecraft and the gravitational parameter of the central body, 
%   this function calculates the position and velocity vectors of the
%   spacecraft at another time.
%
%   r0: initial position vector (km)
%   v0: initial velocity vector (km)
%   mu: standard gravitational parameter of the central body (km^3/s^2)
%   dt: propagation time (days)
%
%   Arguments r0, v0, and mu can be given in any units, but dt must be
%   given in days, which is then converted internally to seconds. 
%   The orbit can be propagated backwards in time.

% if no propagation happens, simply return the input position and velocity
% this is necessary to prevent division by zero errors
if dt == 0
    r1 = r0;
    v1 = v0;
    return
end

% semi-major axis (1.53 or p.51)
a = norm(r0)/(2 - norm(r0)*dot(v0,v0)/mu);
% eccentricity (3.2)
e_vec = ((dot(v0,v0)-mu/norm(r0))*r0 - dot(r0,v0)*v0)/mu;
e = norm(e_vec);

dT = dt * 24 * 3600;

syms x

alpha = 2/norm(r0) - norm(v0)^2/mu;

temp1 = alpha*x^2;

r = r0 * ( 1 - x^2/norm(r0)*C(temp1)) + ...
    v0 * (dT - x^3/sqrt(mu)*S(temp1));
v = r0 * x*sqrt(mu)/norm(r)/norm(r0) * (temp1*S(temp1)-1) + ...
    v0 * ( 1 - x^2/norm(r) *C(temp1));

s0 = dot(r0,v0)/sqrt(mu);

F = s0*x^2*C(temp1)+(1-norm(r0)*alpha)*x^3*S(temp1)+norm(r0)*x-sqrt(mu)*dT;
Fp = s0*x*(1-temp1*S(temp1))+(1-norm(r0)*alpha)*x^2*C(temp1)+norm(r0);
Fpp = diff(Fp, x);

rp = a * (1-e);
x0 = 0;
x1 = mu*dT^2/(rp *(double(subs(F,x,sqrt(mu)*dT/rp))+sqrt(mu)*dT));
n = 6;

while abs(x1-x0) > 0.001
    x0 = x1;
    F0 = double(subs(F, x, x0));
    F1 = double(subs(Fp, x, x0));
    F2 = double(subs(Fpp, x, x0));
    x1 = x0-n*F0/(F1+sign(F1)*sqrt(abs((n-1)^2*F1^2-n*(n-1)*F0*F2)));
end

r1 = double(subs(r,x,x1));
v1 = double(subs(v,x,x1));

end

function [out] = C(y)
%C(y) Prussing & Conway eq. 2.40a

out = (sign(y)/2+0.5) * ((1 - cos(sqrt(y))) / y) - ...
      (sign(y)/2-0.5) * (( cosh(sqrt(-y)) - 1 ) / (-y));

end

function [out] = S(y)
%C(y) Prussing & Conway eq. 2.40b

out = (sign(y)/2+0.5) * ((sqrt(y) - sin(sqrt(y))) / sqrt(y^3)) - ...
      (sign(y)/2-0.5) * (( sinh(sqrt(-y)) - sqrt(-y)) / sqrt(-y^3));

end
