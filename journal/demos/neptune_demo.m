%% neptune_demo.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Predicts VIOD error for an Earth-Neptune transfer orbit by scaling
%   results from a baseline orbit using VIOD error trends.

%% Initialization
close all
clear;clc
addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
addpath('..\fcns_misc')
savePath = 'plots\';
latexify

%% Monte Carlo Simulation for Baseline Case
mu = 1;
a = 1e5;
e = 0.5;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(90);
orbitParams = [a,e,i,o,w,f];

noise = 3e-6;
dM = 0.1 * (2*pi);
numObsv = 10;
numSims = 3000;
selObsv = 1;

nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* nGauss;

errDat = nan(numSims,1);
v = nan(numObsv,3);

% find measurement positions by true anomaly
E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M0 = E0 - e*sin(E0);
M = M0 + dM;
Mvect = linspace(M0,M,numObsv);
Evect = kepler(Mvect,e);
fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
fvect = mod(fvect,2*pi);
fend = fvect(end);
df = fend-f;
df = mod(df,2*pi);

% ground truth data
for i = 1:numObsv
    orbitParams(6) = fvect(i);
    [~,v(i,:)] = Get_Orb_Vects(orbitParams,mu);
end

% get position reference
orbitParams(6) = fvect(selObsv);
rRef = Get_Orb_Vects(orbitParams,mu);

for s = 1:numSims
    nvect = ncube(:,:,s);
    r = hodoHyp(v+nvect,mu);
    r = r(selObsv,:)';
    errDat(s)  = norm(r-rRef) / norm(rRef);
end

baseline.MSE = sqrt(mean(errDat.^2));
baseline.mu = mu;
baseline.orbitParams = orbitParams;
baseline.noise = noise;
baseline.dM = dM;
baseline.numObsv = numObsv;
baseline.df = df;

disp(['Baseline Error: ' num2str(baseline.MSE*100) '%'])

%% Convert Earth-Neptune Case to Canonical Units
AU = 1.495978e11; % m

mu = 1.327e20; % m^3/s^2

a_tgt = 30.06; % Neptune
a = (1+a_tgt)/2 * AU;
e = a_tgt*2/(1+a_tgt)-1;

i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(170);

noise = 5; % m/s
    % currently we can;t predict changes in these two variables,
    %   so they are held constant
    % dM = 0.1 * (2*pi);
    % numObsv = 10;
    % numSims = 3000;
selObsv = 1;

DU = a / 1e5;
TU = sqrt(DU^3/mu);

mu = 1;
a = a / DU;
noise = noise / DU * TU;

orbitParams = [a,e,i,o,w,f];

%% Monte Carlo Simulation for Earth-Neptune Case
nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* nGauss;

errDat = nan(numSims,1);
v = nan(numObsv,3);

% find measurement positions by true anomaly
E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M0 = E0 - e*sin(E0);
M = M0 + dM;
Mvect = linspace(M0,M,numObsv);
Evect = kepler(Mvect,e);
fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
fvect = mod(fvect,2*pi);
fend = fvect(end);
df = fend-f;
df = mod(df,2*pi);

% ground truth data
for i = 1:numObsv
    orbitParams(6) = fvect(i);
    [~,v(i,:)] = Get_Orb_Vects(orbitParams,mu);
end

% get position reference
orbitParams(6) = fvect(selObsv);
rRef = Get_Orb_Vects(orbitParams,mu);

for s = 1:numSims
    nvect = ncube(:,:,s);
    r = hodoHyp(v+nvect,mu);
    r = r(selObsv,:)';
    errDat(s)  = norm(r-rRef) / norm(rRef);
end

MSE = sqrt(mean(errDat.^2));

disp(['True Error: ' num2str(MSE*100) '%'])

%% Predict Earth-Neptune Case Error using VIOD Trends
R_baseline = 1/2 * sqrt(baseline.mu/baseline.orbitParams(1)) ...
                 * 2 / sqrt(1-baseline.orbitParams(2)^2);
R_neptune  = 1/2 * sqrt(mu/a) * 2 / sqrt(1-e^2);
delta_R = (1/R_neptune) / (1/R_baseline);

noise_baseline = baseline.noise;
noise_neptune  = noise;
delta_noise = noise_neptune / noise_baseline;

df_baseline = baseline.df;
df_neptune  = df;
delta_df = (1/df_neptune^2) / (1/df_baseline^2);

MSE_predicted = baseline.MSE * delta_R * delta_noise * delta_df;

disp(['Predicted Error: ' num2str(MSE_predicted*100) '%'])