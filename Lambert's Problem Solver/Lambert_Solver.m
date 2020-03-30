function [a1, v1, v2, alpha_0, beta_0] = Lambert_Solver(r1,r2,mu,dt,tol)
%SOLVER Summary of this function goes here
%   Detailed explanation goes here

% chord length
c = abs(norm(r1-r2));
% semiperimeter of space triangle
s = (c + norm(r1) + norm(r2)) / 2;
% minimum semi-major axis
am = s / 2;
% angle
theta = acos(dot(r1,r2)/norm(r1)/norm(r2));

syms a
alpha = 2 * asin(sqrt(s/2/a));
beta  = 2 * asin(sqrt((s-c)/2/a));
s_alpha = sqrt(s*(2*a-s))/a;
s_beta  = sqrt((s-c)*(2*a-s+c))/a;

eqn1 = a^1.5 * (alpha-beta-s_alpha+s_beta) - dt*sqrt(mu);
eqn2 = diff(eqn1, a);
eqn3 = diff(eqn2, a);
n = 4;

fplot(eqn1)

a0 = 0;
a1 = (norm(r1)+norm(r2))/2;
while abs(a1-a0)>tol
%     clc
    a0 = a1;
    F0 = double(subs(eqn1, a, a0));
    F1 = double(subs(eqn2, a, a0));
    F2 = double(subs(eqn3, a, a0));
    a1 = n*F0 / (F1 + sign(F1)*sqrt(abs((n-1)^2*F1^2 - n*(n-1)*F0*F2)));
    disp(a1);
    disp(abs(a1-a0));
end

theta = atan2(norm(cross(r1,r2)), dot(r1,r2));
if 0 <= theta && theta < pi
    beta_0 = double(subs(beta, a, a1));
else
    beta_0 = -double(subs(beta, a, a1));
end

tm = double(subs(a^1.5*(alpha-beta-sin(alpha)+sin(beta)),a,am))/sqrt(mu);
if dt > tm
    alpha_0 = 2*pi - double(subs(alpha, a, a1));
else
    alpha_0 = double(subs(alpha, a, a1));
end

u1 = r1 / norm(r1);
u2 = r2 / norm(r2);
uc = (r2-r1) / c;

A = sqrt(mu/4/a1)*cot(alpha_0/2);
B = sqrt(mu/4/a1)*cot(beta_0/2);

v1 = (B + A)*uc + (B - A)*u1;
v2 = (B + A)*uc - (B - A)*u2;

end

