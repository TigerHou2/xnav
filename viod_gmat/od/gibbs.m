function [V1] = gibbs(R,mu)
%GIBBS Solves the Gibb's problem (OD from 3 velocity vectors).
%   Only provides the first velocity measurement.
%
% Author:
%   Tiger Hou
%
% Arguments:
%   R - km  , N-by-3 position vector
%
% Outputs:
%   V - km/s, N-by-3 velocity vector
%
% Source:
%   https://smallsats.org/2013/01/26/gibbs-method-of-preliminary-orbit-determination/

R = vclean(R);
idx = ceil(size(R,1)/2);

R1 = R(1,:);
R2 = R(idx,:);
R3 = R(end,:);

r1 = norm(R1); r2 = norm(R2); r3 = norm(R3);
% Vector cross products
R12 = cross(R1,R2);
R23 = cross(R2,R3);
R31 = cross(R3,R1);

Nv  = r1*R23 + r2*R31 + r3*R12;
Dv  = R23 + R31 + R12;
Sv  = (r2-r3)*R1 + (r3-r1)*R2 + (r1-r2)*R3;

N   = norm(Nv);
D   = norm(Dv);
% Velocity vector
V1 = (mu/(N*D))^0.5*(cross(Dv,R1)/r1 + Sv);

end

