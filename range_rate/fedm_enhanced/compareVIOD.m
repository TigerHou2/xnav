%% Compare RROD vs VIOD
%
% Author: Tiger Hou
%
% This script compares RROD and VIOD performance for various orbits.

% close all
clear

addpath('../fcns_circ')
addpath('../fcns_orb')
addpath('../../journal/fcns_od')
addpath('../../journal/fcns_misc')

%% Constants

AU = 1.495978e11; % m
mu_sun = 1.327e20; % m^3/s^2
mu_earth = 3.986e14; % m^3/s^2
a_leo = (6378+410)*1e3; % m
a_neptune = 30.03 * AU;
a_mars = 1.524 * AU;
a_earth = 1.00 * AU;

%% Define Orbit

rng(4)

mu = mu_sun;
targets = [a_earth, a_mars];
a = sum(targets) / 2;
e = max(targets) / a - 1;
% e = 0.02;
i = deg2rad(51.6);
o = deg2rad(0);
w = deg2rad(0);
f0 = deg2rad(140.00);
orbitParams = [a, e, i, o, w, f0];
R = sqrt(mu/a/(1-e^2));
period = 2*pi*sqrt(a^3/mu);

noise_1sigma = 0.5;
Moffset = rand;
dM = 2*pi * 0.02;
measPeriod = period*dM/2/pi;

groupByPulsar = false; % see Section: Calculate Measurement Times
memberInterval = 0;     % after each measurement within the same group
groupInterval = 1;      % between groups
numObsvPerPulsar = 10;
numSims = 70;

%% Define Pulsars

pulsars = [ [1, 0, 0]; ... pulsar 1
            [0, 1, 1]; ... pulsar 2
            [1, 1, 4]; ... pulsar 3
%             [-1, 1, -3]; ... pulsar 4
%             [-1, -1, 2]; ... pulsar 5
          ];
pulsarMat = pulsars' ./ vecnorm(pulsars'); % transpose and normalize
pulsarRotationMatrix = rotz(rand*90) * rotx(rand*90) * roty(rand*90);
pulsarMat = pulsarRotationMatrix * pulsarMat;

numPulsars = size(pulsars,1); % counter the number of pulsars

%% Display Basic Info

[dispTimeUnit,dispTimeScale] = get_time_scale(period);
disp(['Orbit Period:       ' ...
    num2str(period/dispTimeScale) dispTimeUnit])
[dispTimeUnit,dispTimeScale] = get_time_scale(measPeriod);
disp(['Measurement Period: ' ...
    num2str(measPeriod/dispTimeScale) dispTimeUnit])

%% Calculate Measurement Times

E0 = 2*atan( tan(f0/2) * sqrt((1-e)/(1+e)) );
M0 = E0 - e*sin(E0);
Mf = M0 + dM;
if groupByPulsar % is true
    % take measurements 1, 2, 3, ... for puslar A
    % then take measurements 1, 2, 3, ... for pulsar B, ... etc.
    numMembers = numObsvPerPulsar;
    numGroups  = numPulsars;
else
    % take measurement 1 for pulsars A, B, C, ...
    % then take measurement 2 for pulsars A, B, C, ... etc.
    numMembers = numPulsars;
    numGroups  = numObsvPerPulsar;
end
% space out measurements according to member and group intervals
Mvect = [ 0 : (numGroups -1) ] * memberInterval * numMembers ...
      + [ 0 : (numMembers-1) ]'* memberInterval;
Mvect = Mvect + [ 0 : (numGroups -1) ] * groupInterval;
Mvect = Mvect / max(Mvect(:)) * dM + M0;
if groupByPulsar
    Mvect = Mvect';
end
Evect = kepler(Mvect,e);
Fvect = 2*atan( tan(Evect/2) * sqrt((1+e)/(1-e)) );


%% Monte Carlo
errorData_RROD = nan(1,numSims);
errorData_VIOD = nan(1,numSims);
errorData_VIOD_timed = nan(1,numSims);

pos_RROD = nan(3,numSims);
pos_VIOD = nan(3,numSims);
pos_VIOD_timed = nan(3,numSims);
% explanation of noise array:
% noise(i,j,k,m) where
% m = the m-th Monte Carlo sample
% i = the i-th pulsar
% j = the j-th measurement
% k = the k-th pulsar, for pulsar measurement (i,j,m)
%   this is needed because at all time instances, all pulsars have noise.
%   While RROD does not use all pulsars at all times, we need this to
%   construct the full velocity measurement for VIOD to compare.
%     Furthermore, when all pulsars have measurements taken simulteneously
%   (i.e. ~groupByPulsar && memberInterval == 0), the noise for the k-th
%   pulsar at measurement j should be the same regardless of which pulsar
%   is being used for RROD. Therefore, for (j,k,m), the noise should be
%   the same across all i values.
if ~groupByPulsar && memberInterval == 0
    noise = normrnd(0,noise_1sigma,...
                [1,numObsvPerPulsar,numPulsars,numSims]);
    noise = repmat(noise,numPulsars);
else
    noise = normrnd(0,noise_1sigma,...
                [numPulsars,numObsvPerPulsar,numPulsars,numSims]);
end
if all(noise == 0)
    warning("Warning: No perturbations applied.")
end

%%

p = gcp('nocreate');
if ~isempty(p)
    pctRunOnAll warning('off','MATLAB:singularMatrix')
    pctRunOnAll warning('off','MATLAB:nearlySingularMatrix')
    pctRunOnAll warning('off','MATLAB:rankDeficientMatrix')
end
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

p1 = pulsarMat(:,1);
p2 = pulsarMat(:,2);
p3 = pulsarMat(:,3);

% profile on
tic

parfor ii = 1:numSims
    
disp(ii)

    %% Simulate Measurements and Noise

RROD_data = nan(numPulsars,numObsvPerPulsar);
VIOD_data = nan(numGroups,3);
Tvect_RROD = (Mvect+Moffset) * sqrt(a^3/mu);
Tvect_VIOD = nan(numGroups,1);
for j = 1:numPulsars
    thisPulsar = pulsarMat(:,j);
    for k = 1:numObsvPerPulsar
        thisParam = orbitParams;
        thisParam(6) = Fvect(j,k);
        [~,v] = Get_Orb_Vects(thisParam,mu);
        v1 = p1 * (p1' * v + noise(j,k,1,ii));
        v2 = p2 * (p2' * v + noise(j,k,2,ii));
        v3 = p3 * (p3' * v + noise(j,k,3,ii));
        vn = [v1'; v2'; v3'] \ [v1'*v1; v2'*v2; v3'*v3];
        RROD_data(j,k) = thisPulsar' * v + noise(j,k,j,ii);
        
        if ~groupByPulsar && j == 1
            VIOD_data(k,:) = vn';
            Tvect_VIOD(k) = Tvect_RROD(j,k);
        elseif groupByPulsar && k == 1
            VIOD_data(j,:) = vn';
            Tvect_VIOD(k) = Tvect_RROD(j,k);
        end
    end
end

    %% Verification Data

[rTrue,~] = Get_Orb_Vects(orbitParams,mu);

    %% Orbit Determination

[rRROD,vRROD] = rrod(RROD_data,Tvect_RROD,pulsarMat,mu);
[rVIOD,vVIOD] = hodoHyp(VIOD_data,mu);
[rVIOD_timed,vVIOD_timed] = hodo_timed(VIOD_data,Tvect_VIOD,mu);
rVIOD = rVIOD(1,:)';
vVIOD = vVIOD(1,:)';
rVIOD_timed = rVIOD_timed(1,:)';
vVIOD_timed = vVIOD_timed(1,:)';

pos_RROD(:,ii) = rRROD;
pos_VIOD(:,ii) = rVIOD;
pos_VIOD_timed(:,ii) = rVIOD_timed;

    %% Compare Errors

error_RROD = norm(rRROD-rTrue) / norm(rTrue) * 100;
error_VIOD = norm(rVIOD-rTrue) / norm(rTrue) * 100;
error_VIOD_timed = norm(rVIOD_timed-rTrue) / norm(rTrue) * 100;

    %% Store Errors
errorData_RROD(ii) = error_RROD;
errorData_VIOD(ii) = error_VIOD;
errorData_VIOD_timed(ii) = error_VIOD_timed;

end

toc
% profile viewer

p = gcp('nocreate');
if ~isempty(p)
    pctRunOnAll warning('on','MATLAB:singularMatrix')
    pctRunOnAll warning('on','MATLAB:nearlySingularMatrix')
    pctRunOnAll warning('on','MATLAB:rankDeficientMatrix')
end
warning('on','MATLAB:singularMatrix')
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:rankDeficientMatrix')

%% Histogram Results

fig = figure;
[rTrue,vTrue] = Get_Orb_Vects(orbitParams,mu);

corder = get(gca,'colororder');

removeByNan = isnan(errorData_RROD) | ...
              isnan(errorData_VIOD) | ...
              isnan(errorData_VIOD_timed) ;
removeByOutlier = isoutlier(errorData_RROD) | ...
                  isoutlier(errorData_VIOD) | ...
                  isoutlier(errorData_VIOD_timed) ;
% toRemove = removeByNan | removeByOutlier;
toRemove = removeByNan;

subplot(3,1,1)
temp1 = errorData_RROD * norm(rTrue)/100/1000;
temp1(toRemove) = [];
h1 = histogram(temp1,'FaceColor',corder(1,:));
h1.Normalization = 'probability';
legend('RROD')

subplot(3,1,2)
temp2 = errorData_VIOD * norm(rTrue)/100/1000;
temp2(toRemove) = [];
h2 = histogram(temp2,'FaceColor',corder(2,:));
h2.Normalization = 'probability';
legend('VIOD')

subplot(3,1,3)
temp3 = errorData_VIOD_timed * norm(rTrue)/100/1000;
temp3(toRemove) = [];
h3 = histogram(temp3,'FaceColor',corder(3,:));
h3.Normalization = 'probability';
legend('VIOD w/ time')

scale = 1;
width = min([h1.BinWidth,h2.BinWidth,h3.BinWidth]);
h1.BinWidth = width * scale;
h2.BinWidth = width * scale;
h3.BinWidth = width * scale;

subplot(3,1,1)
xlim1 = xlim;
ylim1 = ylim;
subplot(3,1,2)
xlim2 = xlim;
ylim2 = ylim;
subplot(3,1,3)
xlim3 = xlim;
ylim3 = ylim;

limX = [0,max([xlim1(2),xlim2(2),xlim3(2)])];
limY = [0,max([ylim1(2),ylim2(2),ylim3(2)])] * 1.05;

subplot(3,1,1)
xlim(limX)
% ylim(limY)
subplot(3,1,2)
xlim(limX)
% ylim(limY)
subplot(3,1,3)
xlim(limX)
% ylim(limY)

h = axes(fig,'Visible','off');
h.Title.Visible='on';
h.XLabel.Visible='on';
h.YLabel.Visible='on';
xlabel(h,'Position Estimate Error, km')
ylabel(h,'Probability')


%% Scatter Plots

Size = 9;

figure
hold on
scatter(      pos_VIOD(1,:)-rTrue(1),      pos_VIOD(2,:)-rTrue(2), Size  , 'k', 'filled')
scatter(pos_VIOD_timed(1,:)-rTrue(1),pos_VIOD_timed(2,:)-rTrue(2), Size-1, 'g', 'filled')
scatter(      pos_RROD(1,:)-rTrue(1),      pos_RROD(2,:)-rTrue(2), Size-2, 'r', 'filled')
hold off
xlabel('X, m')
ylabel('Y, m')
axis equal
latexify(15,15)
legend('VIOD', 'VIOD w/ time', 'RROD')

figure
hold on
scatter(      pos_VIOD(2,:)-rTrue(2),      pos_VIOD(3,:)-rTrue(3), Size  , 'k', 'filled')
scatter(pos_VIOD_timed(2,:)-rTrue(2),pos_VIOD_timed(3,:)-rTrue(3), Size-1, 'g', 'filled')
scatter(      pos_RROD(2,:)-rTrue(2),      pos_RROD(3,:)-rTrue(3), Size-2, 'r', 'filled')
hold off
xlabel('Y, m')
ylabel('Z, m')
axis equal
latexify(15,15)
legend('VIOD', 'VIOD w/ time', 'RROD')

figure
hold on
scatter(      pos_VIOD(3,:)-rTrue(3),      pos_VIOD(1,:)-rTrue(1), Size  , 'k', 'filled')
scatter(pos_VIOD_timed(3,:)-rTrue(3),pos_VIOD_timed(1,:)-rTrue(1), Size-1, 'g', 'filled')
scatter(      pos_RROD(3,:)-rTrue(3),      pos_RROD(1,:)-rTrue(1), Size-2, 'r', 'filled')
hold off
xlabel('Z, m')
ylabel('X, m')
axis equal
latexify(15,15)
legend('VIOD', 'VIOD w/ time', 'RROD')


%% Function Definitions
function [dispTimeUnit,dispTimeScale] = get_time_scale(period)

if period > 2*60*60*24*365
    dispTimeUnit = ' years';
    dispTimeScale = 60*60*24*365;
elseif period > 2*60*60*24
    dispTimeUnit = ' days';
    dispTimeScale = 60*60*24;
elseif period > 2*60*60
    dispTimeUnit = ' hours';
    dispTimeScale = 60*60;
elseif period > 2*60
    dispTimeUnit = ' minutes';
    dispTimeScale = 60;
else
    dispTimeUnit = ' seconds';
    dispTimeScale = 1;
end

end