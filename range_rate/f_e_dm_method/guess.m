function fVal = guess(OPT,obsv,pulsar,mu,time,debug)
%GUESS returns the sine wave fitting result for RROD.
%
% Author:
%   Tiger Hou
%
% Note:
%   This version guesses the initial true anomaly, eccentricity, and mean
%   anomaly span of the data.
%
% Arguments:
%   OPT    - 1x3 array, optimization variable
%            |- OPT(1) - starting true anomaly
%            |- OPT(2) - eccentricity
%            |- OPT(3) - mean anomaly span
%   obsv   - pxq array of range-rate measurements
%            |- p is the number of pulsars (>=3)
%            |- q is the number of measurements per pulsar (>=3)
%   pulsar - 3xp matrix of pulsar directions
%            |- p is the number of pulsars
%   mu     - gravitational parameter of central body
%   time   - MxN array of observation time stamps since first obsv
%            |- p is the number of pulsars
%            |- q is the number of measurements per pulsar
%

f0 = OPT(1);
e = OPT(2);
dM = OPT(3);

if e < 1, E0 = 2* atan(sqrt((1-e)/(1+e))*tan(f0/2));
else,     E0 = 2*atanh(sqrt((e-1)/(e+1))*tan(f0/2)); end

if e < 1, M0 = E0 - e*sin(E0);
else,     M0 = e*sinh(E0) - E0; end

dt = max(time(:))-min(time(:));
period = dt*2*pi/dM;
M_obsv = (time-min(time(:)))/period*2*pi + M0;
E_obsv = kepler(M_obsv,e);
if e < 1, f_obsv = 2 * atan( tan(E_obsv/2) * sqrt((1+e)/(1-e)) );
else,     f_obsv = 2 * atan(tanh(E_obsv/2) * sqrt((e+1)/(e-1)) ); end

% number of pulsars
p = size(pulsar,2);

% preallocate matrix of range-rates derived from sine waves
numHodoPoints = 6;
angles = linspace(-pi/2,pi/2,numHodoPoints);
rrFromSine = nan(p,numHodoPoints);

% fit range-rate and true anomaly data to sine waves for each pulsar
fVal = 0;

for i = 1:p
    % find number of observations for this pulsar
    q = sum(~isnan(obsv(i,:)));
    % construct matrices for solving linear system for sine wave
    % --- here we break down [[ y = amp * sin(f-pha) + off ]]
    %     as [[ y = amp * sin(f)cos(pha) - amp * cos(f)sin(pha) + off ]]
    %     where sin(f) and cos(f) are known
    
    A = [-sin(f_obsv(i,:))' cos(f_obsv(i,:))'+e];
    b = obsv(i,:)';
    % solve the least squares linear system
    x = A\b;
    off = x(2) * e;        % offset
    pha = atan(x(2)/x(1)); % phase
    amp = -x(2) / sin(pha); % amplitude

    % get points on the sine wave to use later to define hodograph
    rrFromSine(i,:) = amp * sin(angles-pha) + off;
    
    % sine fitting error
    fVal = fVal + abs(A*x-b);
    
    % debugging
    if exist('debug','var')
        if strcmp(debug,'Debug')
            colors = {'--r','--g','--b'};
        else
            colors = {'r','g','b'};
        end
        figure(1)
        hold on
        lb = min(min(f_obsv(:)),0);
        rb = max(max(f_obsv(:)),2*pi);
        lb = min(f_obsv(:));
        rb = max(f_obsv(:));
        fplot(@(x) amp * sin(x-pha) + off, [lb,rb], colors{i})
        hold off
    end
end

if exist('debug','var')
    figure(1)
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
        E = E ...
          - n*(E-e*sin(E)-M) ...
             /(1-e*cos(E)+sign(1-e*cos(E)) ...
                         *sqrt(abs((n-1)^2*(1-e*cos(E))^2 ...
                                  -(n-1)*n*(E-e*sin(E)-M)*(e*sin(E)))));
        d = abs(E-e*sin(E)-M);
    end
else
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E ...
          - n*(e*sinh(E)-E-M) ...
             /(e*cosh(E)-1+sign(e*cosh(E)-1) ...
                          *sqrt(abs((n-1)^2*(e*cosh(E)-1)^2 ...
                                   -(n-1)*n*(e*sinh(E)-E-M)*(e*sinh(E)))));
        d = abs(e*sinh(E)-E-M);
    end
end

Evect(i) = E;

end
end