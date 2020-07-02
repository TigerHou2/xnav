%% plotSpcraftParams.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Generates plots for velocity initial orbit determination (VIOD) error
%   with respect to variations in spacecraft parameters, including
%   measurement noise and measurement scheduling.
%
% Notes:
%   In the case of varying spacecraft parameters, we choose to hold orbital
%   parameters constant. This means all simulations will be conducted with
%   equal semi-major axis, eccentricity, and true anomaly. The orbital
%   parameters chosen are as follows:
%       - SMA = 1e5 DU
%       - ECC = 0.5
%       - TA  = 20 deg

%% initialization

close all hidden
clear;clc;init

% load spacecraft parameter config case studies
[noiseVect,durVect,obsVect] = load_spcraft_cases();

% define gravitational parameter as 1 in canonical units
mu = 1; % DU^3/TU^2

% define orbital parameters
SMA  = 1e5; % DU
ECC  = 0.5; % nd
INC  = 45;  % deg
RAAN = 30;  % deg
AOP  = 30;  % deg
TA   = 20;  % deg

% convert to radians and combine parameters
orbitParams = [SMA,ECC,deg2rad([INC,RAAN,AOP,TA])];

% Monte Carlo settings
numSims = 100;
rngSeed = 1;

%% iterate through all cases

% plotting formats
linewidth = 1;
origFormat = '-k+';
hodoFormat = '-r+';

for i = 1:length(durVect)
    duration = durVect(i);
for j = 1:length(obsVect)
    numObsv = obsVect(j);
    
% initialize position error matrix
errOrig = nan(numSims,length(noiseVect));
errHodo = nan(numSims,length(noiseVect));
% reset random number generator
rng(rngSeed);
% calculate mean anomaly
period = 2*pi*sqrt(SMA^3/mu);
e = orbitParams(2);
f = orbitParams(6);
E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);
% determine measurement locations
Mvect = M + linspace(0,duration,numObsv) * 2*pi;
Evect = kepler(Mvect,e);
fVect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
% store ground truth position vector at starting true anomaly
rRef = Get_Orb_Vects(orbitParams,mu);

% perform Monte Carlo sim
for n = 1:numSims
noiseBase = randn(numObsv,3);
noiseBase = noiseBase ./ vecnorm(noiseBase,2,2) ;

for k = 1:length(noiseVect)
    % load new case
    noise = noiseBase * noiseVect(k);
    % initialize empty velocity matrix
    v = nan(numObsv,3);
    % fetch all measurements
    for p = 1:numObsv
        orbitParams(6) = fVect(p);
        [~,v(p,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    % restore true anomaly in orbit parameter
    orbitParams(6) = f;
    % do orbit determination
    % --- note the scaling dor the original method: this is because
    % --- precision issues arise when using canonical units. 
    rOrig = viod((v+noise)*1e4,mu*1e12)/1e4;
    rHodo = hodo(v+noise,mu);
    % we choose to compare the position estimate at the first measurement
    rOrig = rOrig(1,:)';
    rHodo = rHodo(1,:)';
    % store error data
    errOrig(n,k) = norm(rOrig-rRef);
    errHodo(n,k) = norm(rHodo-rRef);
end %nVect

end %numSims

% plot results for each case
figure;
noiseScale = 1e-6;
hold on
errorbar(noiseVect/noiseScale,...
         mean(errOrig)/SMA,...
          std(errOrig)/SMA,origFormat,'LineWidth',linewidth)
errorbar(noiseVect/noiseScale,...
         mean(errHodo)/SMA,...
          std(errHodo)/SMA,hodoFormat,'LineWidth',linewidth)
hold off
setgrid
expand(0.05,0.05,0.05,0.05)
legend('Energy Method','Hodograph Method','Location','NorthWest')
xlabel(['Noise, ' num2str(noiseScale) ' DU/TU'])
ylabel('Position Error, fraction of SMA')
title([ 'DUR = ' num2str(duration,4) ', ' ...
        'Observations = ' num2str(numObsv,4)])
latexify

end %obsVect
end %durVect