function [optDiff,V,optOut] = rrFun_sine(OPT,obsv,pulsar,mu,time,debug)
%RRFUN_SINE Is the objective function for range-rate hodograph fitting.
%
% Author:
%   Tiger Hou
%
% Note: 
%   This version uses the fact that range-rate vs. true anomaly is a sine
%   wave to perform RROD. 
%
% Arguments:
%   OPT    - 1x3 array, optimization variable
%            |- OPT(1) - eccentricity
%            |- OPT(2) - period
%            |- OPT(3) - time since periapsis
%   obsv   - MxN array of range-rate measurements
%            |- M is the number of pulsars (>=3)
%            |- N is the number of measurements per pulsar (>=3)
%   pulsar - 3xM matrix of pulsar directions
%            |- M is the number of pulsars
%   mu     - gravitational parameter of central body
%   time   - MxN array of observation time stamps since first obsv
%            |- M is the number of pulsars
%            |- N is the number of measurements per pulsar
%
% Outputs:
%   v_diff - range-rate difference between guess and true values
%   V      - velocities corresponding to hodograph defined by C

% take data from optimization variables
g_e = OPT(1); % eccentricity
g_period = OPT(2); % period
g_timepe = OPT(3); % time since periapsis

% get number of pulsars
M = size(pulsar,2);

% convert time of flight to true anomaly
g_tpArray = time - g_timepe;
g_M = mod( g_tpArray/g_period*2*pi, 2*pi );
g_E = kepler(g_M,g_e);
if g_e < 1, g_f = 2 * atan( tan(g_E/2) * sqrt((1+g_e)/(1-g_e)) );
else,       g_f = 2 * atan(tanh(g_E/2) * sqrt((g_e+1)/(g_e-1)) ); end

% preallocate matrix of range-rates derived from sine waves
numHodoPoints = 6;
angles = linspace(-pi/2,pi/2,numHodoPoints);
rrFromSine = nan(M,numHodoPoints);

% fit range-rate and true anomaly data to sine waves for each pulsar
for i = 1:M
    % find number of observations for this pulsar
    N = sum(~isnan(obsv(i,:)));
    % construct matrices for solving linear system for sine wave
    % --- here we break down [[ y = amp * sin(f-pha) + off ]]
    %     as [[ y = amp * sin(f)cos(pha) - amp * cos(f)sin(pha) + off ]]
    %     where sin(f) and cos(f) are known
    A = [sin(g_f(i,:))' -cos(g_f(i,:))' ones(N,1)];
    b = obsv(i,:)';
    % solve the least squares linear system
    x = A\b;
    pha = atan(x(2)/x(1)); % phase
    amp = x(2) / sin(pha); % amplitude
    off = x(3);            % offset
    % get points on the sine wave to use later to define hodograph
    rrFromSine(i,:) = amp * sin(angles-pha) + off;
    % debugging
    if exist('debug','var')
        figure(1)
        hold on
        lb = min(min(g_f(:)),0);
        rb = max(max(g_f(:)),2*pi);
        fplot(@(x) amp * sin(x-pha) + off, [lb,rb])
        hold off
    end
end

% calculate the intersection of planes defined by range rates from sines
% --- create 3xQ matrix where
%     Q is the number of intersections used to create hodograph circle
hodoPoints = nan(3,numHodoPoints);
for k = 1:numHodoPoints
    % construct linear system to solve for intersection of planes
    % --- where the plane is defined by the pulsar vector [a0,b0,c0]
    %     and the point [x0,y0,z0]
    %     so the plane has a0(x-x0) + b0(y-y0) + c0(z-z0) = 0
    % --- transforming to a linear system, we have
    %     [a1  b1  c1]  [x]     [ a1*x1 + b1*y1 + c1*z1 ]
    %     [a2  b2  c3]  [y]  =  [ a2*x2 + b2*y2 + c2*z2 ]
    %     [a3  b3  c3]  [z]     [ a3*x3 + b3*y3 + c3*z3 ]
    A = pulsar';
    b = sum( rrFromSine(:,k).*pulsar'.*pulsar'./vecnorm(pulsar',2,2), 2 );
    hodoPoints(:,k) = A\b;
end

% use singular value decomposition to find hodograph plane normal
[~,~,V] = svd(hodoPoints',0);
normal = V(:,end);
% --- check if the sign is correct by comparing against cross product of
% velocities at f = -90 and f = 0
if normal'*cross(hodoPoints(:,1),hodoPoints(:,2)) < 0
    normal = -normal;
end

% find the transformation matrix to bring the hodograph points to 2D
u1 = hodoPoints(:,1) / norm(hodoPoints(:,1));
u2 = cross(normal,u1);
T = [u1';u2';normal'];

% transform hodograph points to 2D, then fit circle
hodo2D = T * hodoPoints; % third row should be zero
[C,R,residual] = fitcircle(hodo2D(1:2,:));

% convert hodograph center back to 3D
C = T' * [C;0];

% calculate eccentricity and period
e = norm(C) / R;
vp = R - norm(C);
a = mu / vp^2 * (1+e^2-2*e)/(1-e^2);
period = 2*pi * sqrt(abs(a^3/mu));

if exist('debug','var')
    % output full velocity measurements
    V = nan(size(obsv,1),size(obsv,2),3);
    radius = C*R / norm(C);
    for i = 1:M
        for j = 1:N
            df = g_f(i,j);
            V(i,j,:) = rotVec(radius,normal,df) + C;
        end
    end
else
    V = nan;
end

% finalizing outputs
optDiff = [e,period,residual] - [g_e,g_period,0];
weight = [50 1 10];
scale = 100;
optDiff = optDiff .* weight / norm(weight) * scale;
optOut  = [e,period,g_timepe];

% if the time since periapse is out of bounds (i.e. < 0 or > period)
% then we add a penalty for being a periodic solution that is out of range
if mod(g_timepe,g_period) ~= g_timepe
    optDiff = optDiff * (abs(mod(g_timepe,g_period)-g_timepe)/g_period+1);
end

end

function Evect = kepler(Mvect,e)

Evect = nan(size(Mvect));

for i = 1:length(Mvect(:))
    
M = Mvect(i);

if M == 0
    Evect(i) = 0;
    continue
end

d = 1;
n = 5;

if e < 1
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(E-e*sin(E)-M) / ...
            (1-e*cos(E) + sign(1-e*cos(E)) ...
                        * sqrt(abs((n-1)^2*(1-e*cos(E))^2 ...
                                   -n*(n-1)*(E-e*sin(E)-M)*(e*sin(E)))));
        d = abs(E-e*sin(E)-M);
    end
else
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(e*sinh(E)-E-M) / ...
            (e*cosh(E)-1 + sign(e*cosh(E)-1) ...
                         * sqrt(abs((n-1)^2*(e*cosh(E)-1)^2 ...
                                    -n*(n-1)*(e*sinh(E)-E-M)*(e*sinh(E)))));
        d = abs(e*sinh(E)-E-M);
    end
end

Evect(i) = E;

end

end