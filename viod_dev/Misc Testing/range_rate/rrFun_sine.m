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
g_Moffset = OPT(3); % time since periapsis

% get number of pulsars
M = size(pulsar,2);

% convert time of flight to true anomaly
g_M = mod( time/g_period*2*pi-g_Moffset, 2*pi );
g_E = kepler(g_M,g_e);
if g_e < 1, g_f = 2 * atan( tan(g_E/2) * sqrt((1+g_e)/(1-g_e)) );
else,       g_f = 2 * atan(tanh(g_E/2) * sqrt((g_e+1)/(g_e-1)) ); end

% preallocate matrix of range-rates derived from sine waves
numHodoPoints = 6;
angles = linspace(-pi/2,pi/2,numHodoPoints);
rrFromSine = nan(M,numHodoPoints);

% fit range-rate and true anomaly data to sine waves for each pulsar
fVal = 0;
for i = 1:M
    % find number of observations for this pulsar
    N = sum(~isnan(obsv(i,:)));
    % construct matrices for solving linear system for sine wave
    % --- here we break down [[ y = amp * sin(f-pha) + off ]]
    %     as [[ y = amp * sin(f)cos(pha) - amp * cos(f)sin(pha) + off ]]
    %     where sin(f) and cos(f) are known
    if g_e ~= 0
        A = [-sin(g_f(i,:))'/g_e cos(g_f(i,:))'/g_e+1];
        b = obsv(i,:)';
        % solve the least squares linear system
        x = A\b;
        off = x(2);        % offset
        pha = atan(off/x(1)); % phase
        amp = -off/g_e/sin(pha); % amplitude
        fVal = fVal + norm(A*x-b);
    else
        A = [sin(g_f(i,:))' -cos(g_f(i,:))'];
        b = obsv(i,:)';
        % solve the least squares linear system
        x = A\b;
        pha = atan(x(2)/x(1)); % phase
        amp = x(2) / sin(pha); % amplitude
        off = 0;            % offset
        fVal = fVal + norm(A*x-b);
    end
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

% if exist('debug','var')
%     figure(9)
%     scatter3(hodoPoints(1,:),hodoPoints(2,:),hodoPoints(3,:))
%     axis equal
% end

%     % here, we have a problem: the circle found by sine waves may not all
%     % contain the origin, which is required for the hodograph
%     % for now, we will simply shift the circle.
%     % --- use vectors from adjacent points on circle to find normal
%     hodoAdj = hodoPoints(:,2:end) - hodoPoints(:,1:end-1);
%     [~,~,V] = svd(hodoAdj',0);
%     normal = V(:,end);
%     % --- check if the sign is correct by comparing against cross product of
%     % adjacent vectors in hodoAdj
%     if normal'*cross(hodoAdj(:,1),hodoAdj(:,2)) < 0
%         normal = -normal;
%     end
%     % --- now, we will offset hodoPoints along the normal
%     offset = mean(normal'*hodoPoints);
%     hodoPoints = hodoPoints - offset;

    % use singular value decomposition to find hodograph plane normal
    [~,~,V] = svd(hodoPoints',0);
    normal = V(:,end);
    % --- check if the sign is correct by comparing against cross product of
    % velocities at adjecent true anomalies
    if normal'*cross(hodoPoints(:,1),hodoPoints(:,2)) < 0
        normal = -normal;
    end

% if exist('debug','var')
%     figure(10)
%     scatter3(hodoPoints(1,:),hodoPoints(2,:),hodoPoints(3,:))
%     hold on
%     quiver3(0,0,0,normal(1),normal(2),normal(3))
%     hold off
%     axis equal
% end

% find the transformation matrix to bring the hodograph points to 2D
u1 = hodoPoints(:,1) / norm(hodoPoints(:,1));
u2 = cross(normal,u1);
T = [u1';u2';normal'];

% transform hodograph points to 2D
hodo2D = T * hodoPoints; % third row should be zero if no error

% construct linear system [1] eqn.9
A = 2*hodo2D';
A(:,3) = -1;
B = (hodo2D').^2;
B = sum(B(:,1:2),2);
x = A\B;
% find radius of hodograph
R = sqrt(x(1)^2 + x(2)^2 - x(3));
% find center of hodograph
C = T' * [x(1); x(2); 0];
% find residual
residual = norm(A*x-B);
% add radius error to residual
g_sma = ((g_period/2/pi)^2*mu)^(1/3);
g_R = (sqrt(mu*(1+g_e)/g_sma/(1-g_e))+sqrt(mu*(1-g_e)/g_sma/(1+g_e)))/2;
residual = residual + (R-g_R)^2;

% % circle fitting with known radius
% g_sma = ((g_period/2/pi)^2*mu)^(1/3);
% R = ( sqrt(mu*(1+g_e)/g_sma/(1-g_e)) + sqrt(mu*(1-g_e)/g_sma/(1+g_e)) )/2;
% A = hodo2D';
% A(:,3) = [];
% fun = @(x) norm(vecnorm(A-x,2,2).^2-R^2);
% options = optimoptions('fminunc','Display','none');
% C = fminunc(fun,mean(A),options);
% residual = fun(C) / 1000;
% C = T' * [C,0]';

if exist('debug','var')
    figure(11)
    scatter3(hodo2D(1,:),hodo2D(2,:),hodo2D(3,:))
    axis equal
end

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
optDiff = [e,period,residual+fVal] - [g_e,g_period,0];
weight = [5 1 2];
scale = 100;
optDiff = optDiff .* weight / norm(weight) * scale;
optOut  = [e,period,g_Moffset];

% if the time since periapse is out of bounds (i.e. < 0 or > period)
% then we add a penalty for being a periodic solution that is out of range
if mod(g_Moffset,2*pi) ~= g_Moffset
    optDiff = optDiff * (abs(mod(g_Moffset,2*pi)-g_Moffset)/2/pi+1);
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