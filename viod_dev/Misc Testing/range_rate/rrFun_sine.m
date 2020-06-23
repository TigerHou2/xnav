function [optDiff,V] = rrFun_sine(OPT,obsv,pulsar,mu,time)
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
g_tpArray = time + g_timepe;
g_M = mod( g_tpArray/g_period*2*pi, 2*pi );
g_E = kepler(g_M,g_e);
if g_e < 1, g_f = 2 * atan( tan(g_E/2) * sqrt((1+g_e)/(1-g_e)) );
else,       g_f = 2 * atan(tanh(g_E/2) * sqrt((g_e+1)/(g_e-1)) ); end

% preallocate matrix of range-rates derived from sine waves
rrFromSine = nan(M,3);

% fit range-rate and true anomaly data to sine waves for each pulsar
for i = 1:M
    % find number of observations for this pulsar
    N = sum(~isnan(obsv(i,:)));
    % construct matrices for solving linear system for sine wave
    A = [cos(g_f(i,:))' sin(g_f(i,:))' ones(N,1)];
    b = obsv(i,:)';
    % solve the least squares linear system
    x = A\b;
    amp = sqrt(x(1)^2+x(2)^2); % amplitude
    off = x(3);                % offset
    pha = atan(x(1)/x(2));     % phase
    % get three points on the sine wave to use later to define circle
    % --- we simply choose f = -90, 0, 90
    rrFromSine(i,:) = [ amp * sin(-pi/2+pha) + off,...
                        amp * sin( 0   +pha) + off,...
                        amp * sin( pi/2+pha) + off];
end

% calculate the intersection of planes defined by range rates from sines
hodoPoints = nan(3); % 3xQ matrix, where Q is the number of intersections
for k = 1:3
    % construct linear system to solve for intersection of planes
    A = pulsar';
    b = sum( rrFromSine(:,k).*pulsar'./vecnorm(pulsar',2,2), 2 );
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
[C,R] = fitcircle(hodo2D(1:2,:));

% convert hodograph center back to 3D
C = T' * [C;0];

% calculate eccentricity, period, and time since periapsis
e = norm(C) / R;
vp = R - norm(C);
a = mu / vp^2 * (1+e^2+2*e)/(1-e^2);
period = 2*pi * sqrt(abs(a^3/mu));
timepe = min(g_f(:)) / 2*pi * period;

% finalizing outputs
V = nan;
optDiff = [e,period,timepe] - [g_e,g_period,g_timepe];

end

function Evect = kepler(Mvect,e)

Evect = nan(size(Mvect));

for i = 1:length(Mvect(:))
    
M = Mvect(i);

if M == 0
    E = 0;
    continue
end

d = 1;
n = 5;

if e < 1
    E = mod(M+e/2,2*pi);
    while d > 1e-6
        E = E - n*(E-e*sin(E)-M) / ...
            (1-e*cos(E) + sign(1-e*cos(E)) ...
                        * sqrt(abs((n-1)^2*(1-e*cos(E))^2 ...
                                   -n*(n-1)*(E-e*sin(E)-M)*(e*sin(E)))));
        d = abs(E-e*sin(E)-M);
    end
else
    E = mod(M+e/2,2*pi);
    while d > 1e-6
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