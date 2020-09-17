%% hodoHyp.m
function [r] = hodoHyp(N,mu)
%HODOHYP Solves the three-velocity initial orbit determination problem.
%
% Author:
%   Tiger Hou
%
% Note:
%   This version uses the hyperaccurate algebraic fit described here:
%   https://people.cas.uab.edu/~mosya/cl/AC1c.pdf
%   and converted to C++ here:
%   https://github.com/SohranEliassi/Circle-Fitting-Hyper-Fit
%
% Description:
%   The spacecraft's orbit is calculated using the velocity hodograph from
%   three velocity vectors collected at different points in orbit.
%
% Limitations:
%   - Assumes the measurements are provided in order.
%   - Assumes 50%+ of adjacent measurements are less than half an orbit
%   away from each other
%
% Arguments:
%   N:      [km/s]      N-by-3 matrix containing N velocity vectors
%   mu:     [km^3/s^2]  gravitational parameter of central body
%
% References:
%   [1] - Geometric Solutions for Problems in Velocity-Based Orbit Determination

% use singular value decomposition to find orbit plane normal
[~,~,V] = svd(N,0);
k = V(:,end);

% check if orbit plane normal is in the right direction
% we assume more than half of the measurements are taken less than half an
% orbit apart from its adjacent measurements.
kEst = zeros(1,3);
for i = 1:size(N,1)-1
    kEst = kEst + cross(N(i,:),N(i+1,:));
end
if dot(k,kEst') < 0
    k = -k;
end

% define rotation matrix
ux = cross(N(1,:),k') / norm(cross(N(1,:),k')); % [1] eqn.4
uy = cross(k',ux);
T = [ux; uy; k']; % [1] eqn.5

% transform velocity from inertial frame to orbit frame
N2d = (T * N')'; % [1] eqn.3

[a,b,R] = hyperfit(N2d(:,1:2)');

% find center of hodograph
c = T' * [a; b; 0]; % [1] eqn.11

% find eccentricity vector
e = cross(c,k) / R; % [1] eqn.18

% find perpendicular unit vector
upp = (N-c') ./ vecnorm(N-c',2,2); % [1] eqn.19

% find parallel unit vector
ull = cross(upp,repmat(k',size(N,1),1),2); % [1] eqn.20

% compute perpendicular velocity components [1] eqn.23
Npp = nan(size(N));
for i = 1:size(N,1)
    Npp(i,:) = ( (upp(i,:)'*upp(i,:)) * N(i,:)' )';
end

% compute position magnitudes [1] eqn.29
rho = mu * vecnorm(e'+ull,2,2) ./ vecnorm(N,2,2) ./ vecnorm(Npp,2,2);

% compute position vectors [1] eqn.30
r = rho .* ull;

end %hodoHyp.m

function [x,y,R] = hyperfit(dat)
%HYPERFIT Summary of this function goes here
%   Detailed explanation goes here

% compute moments
means = mean(dat,2);
dat_c = dat - means;
dat_z = dat_c(1,:).^2 + dat_c(2,:).^2;
Mxx = mean(dat_c(1,:).*dat_c(1,:));
Mxy = mean(dat_c(1,:).*dat_c(2,:));
Myy = mean(dat_c(2,:).*dat_c(2,:));
Mxz = mean(dat_c(1,:).*dat_z(1,:));
Myz = mean(dat_c(2,:).*dat_z(1,:));
Mzz = mean(dat_z(1,:).*dat_z(1,:));

% compute coefficients of characteristic polynomial
Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
Var_z = Mzz - Mz*Mz;

A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
A22 = A2 + A2;

% root finding using Newton's algorithm
x = 0;
y = A0;
while true
    Dy = A1 + x*(A22 + 16*x^2);
    xnew = x - y/Dy;
    if xnew == x || isinf(xnew) || isnan(xnew), break; end
    ynew = A0 + xnew*(A1 + xnew*(A2 + 4*xnew^2));
    if abs(ynew) >= abs(y), break; end
    x = xnew; y = ynew;
end

% calculate parameters of resulting circle
DET = x*x - x*Mz + Cov_xy;
Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;

x = Xcenter + means(1);
y = Ycenter + means(2);
R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);

end