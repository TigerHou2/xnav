function [r1,v1] = TimeProp(r0,v0,mu,dt)
%TIMEPROP Summary of this function goes here
%   Detailed explanation goes here

% semi-major axis (1.53 or p.51)
a = norm(r0)/(2 - norm(r0)*dot(v0,v0)/mu);
% eccentricity (3.2)
e_vec = ((dot(v0,v0)-mu/norm(r0))*r0 - dot(r0,v0)*v0)/mu;
e = norm(e_vec);


if e <= 1

    % eccentric anomaly at t0 (2.12)
    % NOTE the small subtraction at the end of the equation. This was done to
    %   prevent getting imaginary results due to floating point accuracy
    E0 = acos((-norm(r0)/a+1)/e - 10^(-15));

    % (2.7)
    M = E0 - e*sin(E0);
    E1 = 0;
    E2 = M + e/2;

    syms E
    % transcendental Kepler's equation (2.27)
    F = E - e*sin(E) - M - sqrt(mu/a^3)*dt*24*3600;
    Fp = diff(F,E);

    while abs(E1-E2) > 0.001
        E1 = E2;
        F0 = double(subs(F,E,E1));
        F1 = double(subs(Fp,E,E1));
        E2 = E1 - F0/F1;
    end

    E1 = E2;

    % (2.26, 2.28, 2.29)
    f = 1 - a/norm(r0) * (1 - cos(E1 - E0));
    g = dt*24*3600 - sqrt(a^3/mu)*((E1-E0)-sin(E1-E0));

    % (2.17)
    r1 = f*r0 + g*v0;

    fdot = -sqrt(mu*a)/norm(r1)/norm(r0)*sin(E1-E0);
    gdot = 1 - a/norm(r1)*(1-cos(E1-E0));

    % (2.18)
    v1 = fdot*r0 + gdot*v0;

else

    % hyperbolic anomaly
    H0 = acosh((1-norm(r0)/a)/e);

    % (2.7)
    M = e*sinh(H0) - H0;
    H1 = 0;
    H2 = M;

    syms H
    % transcendental Kepler's equation (2.27)
    F = M - e*sinh(H) + H - sqrt(mu/-a^3)*dt*24*3600;
    Fp = diff(F,H);

    while abs(H1-H2) > 0.001
        H1 = H2;
        F0 = double(subs(F,H,H1));
        F1 = double(subs(Fp,H,H1));
        H2 = H1 - F0/F1;
    end

    H1 = H2

    % (2.26, 2.28, 2.29)
    f = 1 - a/norm(r0) * (1 - cos(H1 - H0));
    g = dt*24*3600 - sqrt(-a^3/mu)*(sinh(H1-H0)-(H1-H0));

    % (2.17)
    r1 = f*r0 + g*v0;

    fdot = -sqrt(-mu*a)/norm(r1)/norm(r0)*sin(H1-H0);
    gdot = 1 - a/norm(r1)*(1-cos(H1-H0));

    % (2.18)
    v1 = fdot*r0 + gdot*v0;
end

end
