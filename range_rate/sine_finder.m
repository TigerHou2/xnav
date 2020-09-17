function fVal = sine_finder(OPT,obsv,pulsar,mu,time,debug)
%SINE_FINDER Attempts to find the best fitting sine wave for RROD.
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
        fVal = fVal + abs(A*x-b);
    else
        A = [sin(g_f(i,:))' -cos(g_f(i,:))'];
        b = obsv(i,:)';
        % solve the least squares linear system
        x = A\b;
        pha = atan(x(2)/x(1)); % phase
        amp = x(2) / sin(pha); % amplitude
        off = 0;            % offset
        fVal = fVal + abs(A*x-b);
    end
    % get points on the sine wave to use later to define hodograph
    rrFromSine(i,:) = amp * sin(angles-pha) + off;
    % debugging
    if exist('debug','var')
        colors = ['r','g','b'];
        figure(1)
        hold on
        lb = min(min(g_f(:)),0);
        rb = max(max(g_f(:)),2*pi);
        fplot(@(x) amp * sin(x-pha) + off, [lb,rb], colors(i))
        hold off
    end
end

end

%% function definitions
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