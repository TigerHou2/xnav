%% plotOrbitParams.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Generates plots for velocity initial orbit determination (VIOD) error
%   with respect to variations in orbital parameters, including semi-major
%   axis, eccentricity, and true anomaly.
%
% Notes:
%   In the case of varying orbital parameters, we choose to hold spacecraft
%   parameters constant. This means all simulations will be conducted with
%   equal noise, number of measurements, and duration of observation.
%   Specifically, we will use the following:
%       - Noise = 1e-6 DU/TU
%       - Number of measurements = 3
%       - Duration of observation = 0.1 orbital period

%% initialization

close all hidden
clear;clc;init

% load orbital parameter config case studies
[smaVect,eccVect,taVect,noise,names] = load_orbit_cases();

% define gravitational parameter as 1 in canonical units
mu = 1; % DU^3/TU^2

% define orbital parameters
SMA  = 1e6; % DU
ECC  = 0.1; % nd
INC  = 45;  % deg
RAAN = 30;  % deg
AOP  = 30;  % deg
TA   = 0;   % deg

% convert to radians and combine parameters
orbitParams = [SMA,ECC,deg2rad([INC,RAAN,AOP,TA])];

% spacecraft parameters (fixed)
numObsv = 3; % nd, number of measurements
duration = 0.1; % nd, measurement duration as fraction of the orbit period

% Monte Carlo settings
numSims = 100;
rngSeed = 1;

%% iterate through all cases

% plotting formats
linewidth = 1;
origFormat = '-k+';
hodoFormat = '-r+';

for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
    
% initialize position error matrix
errOrig = nan(numSims,length(taVect));
errHodo = nan(numSims,length(taVect));

% reset random number generator
rng(rngSeed);

% perform Monte Carlo sim
for n = 1:numSims
noiseVect = randn(numObsv,3);
noiseVect = noiseVect ./ vecnorm(noiseVect,2,2) * noise;

for k = 1:length(taVect)
    TA = taVect(k);
    % load new case
    orbitParams(1) = SMA;
    orbitParams(2) = ECC;
    orbitParams(6) = TA;
    % initialize empty velocity matrix
    v = nan(numObsv,3);
    % calculate mean anomaly
    period = 2*pi*sqrt(SMA^3/mu);
    E = 2 * atan(sqrt((1-ECC)/(1+ECC))*tan(TA/2));
    M = E - ECC*sin(E);
    % determine measurement locations
    Mvect = M + linspace(0,duration,numObsv) * 2*pi;
    Evect = kepler(Mvect,ECC);
    fVect = 2 * atan(sqrt((1+ECC)/(1-ECC))*tan(Evect/2));
    % store ground truth position vector at starting true anomaly
    rRef = Get_Orb_Vects(orbitParams,mu);
    % fetch all measurements
    for p = 1:numObsv
        orbitParams(6) = fVect(p);
        [~,v(p,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    % do orbit determination
    % --- note the scaling dor the original method: this is because
    % --- precision issues arise when using canonical units. 
    rOrig = viod((v+noiseVect)*1e4,mu*1e12)/1e4;
    rHodo = hodo(v+noiseVect,mu);
    % we choose to compare the position estimate at the first measurement
    rOrig = rOrig(1,:)';
    rHodo = rHodo(1,:)';
    % store error data
    errOrig(n,k) = norm(rOrig-rRef);
    errHodo(n,k) = norm(rHodo-rRef);
end %taVect

end %numSims

% plot results for each case
figure;
taDeg = rad2deg(taVect);
hold on
% scatter(taMat(:),errOrig(:),scFormat)
errorbar(taDeg,mean(errOrig)/SMA,...
                std(errOrig)/SMA,origFormat,'LineWidth',linewidth)
errorbar(taDeg,mean(errHodo)/SMA,...
                std(errHodo)/SMA,hodoFormat,'LineWidth',linewidth)
hold off
setgrid
expand(0.05,0.05,0.05,0.05)
legend('Energy Method','Hodograph Method','Location','NorthWest')
xlabel('True Anomaly, deg')
ylabel('Position Error, fraction of SMA')
title([ names{i} ':   ' ...
        'SMA = ' num2str(SMA,4) ', ' ...
        'ECC = ' num2str(ECC,4)])
latexify

end %eccVect
end %smaVect