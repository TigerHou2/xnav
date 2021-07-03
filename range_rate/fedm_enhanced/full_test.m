%% f-e-dM method full test
%
% Author: Tiger Hou
%
% this script tests the f0-e-dM range-rate orbit determination function.

close all
clear

addpath('../fcns_circ')
addpath('../fcns_orb')

%% Settings
applyNoise = true;
randomNoise = true;
randomEpoch = true;
randomPulsarRotation = true;

noise = 1;
Moffset = 0.7;
pulsarRotationMatrix = rotz(25) * rotx(12) * roty(2);

%% Define Orbit

mu = 3.986e14;
a = (6378+600)*1e3;
e = 0.13423;
i = deg2rad(45);
o = deg2rad(20);
w = deg2rad(16);
f0 = deg2rad(50.62);
orbitParams = [a, e, i, o, w, f0];
R = sqrt(mu/a/(1-e^2));

%% Define Pulsars

pulsars = [ [1, 0, 0]; ... pulsar 1
            [0, 1, 1]; ... pulsar 2
            [1, 1, 4]; ... pulsar 3
            [-1, 1, -3]; ... pulsar 4
%             [-1, -1, 2]; ... pulsar 5
          ];
pulsarMat = pulsars' ./ vecnorm(pulsars'); % transpose and normalize

if (randomPulsarRotation)
    pulsarRotationMatrix = rotz(rand*90) * rotx(rand*90) * roty(rand*90);
end
pulsarMat = pulsarRotationMatrix * pulsarMat;

numPulsars = size(pulsars,1); % counter the number of pulsars

%% Simulate Measurements

groupMeasurements = true;  % measurements are taken in groups.
                            % If true, each group would contain
                            % all measurements from only one pulsar.
                            % If false, each group would contain one
                            % measurement from each pulsar.
% time interval ratios
memberInterval = 1;     % after each measurement within the same group
groupInterval = 1;      % between groups

numObsvPerPulsar = 3;
dM = 2*pi * 0.12235;
E0 = 2*atan( tan(f0/2) * sqrt((1-e)/(1+e)) );
M0 = E0 - e*sin(E0);
Mf = M0 + dM;
if groupMeasurements
    numMembers = numObsvPerPulsar;
    numGroups  = numPulsars;
else
    numMembers = numPulsars;
    numGroups  = numObsvPerPulsar;
end
Mvect = [ 0 : (numGroups -1) ] * memberInterval * numMembers ...
      + [ 0 : (numMembers-1) ]'* memberInterval;
Mvect = Mvect + [ 0 : (numGroups -1) ] * groupInterval;
Mvect = Mvect / max(Mvect(:)) * dM + M0;
if groupMeasurements
    Mvect = Mvect';
end
Evect = kepler(Mvect,e);
Fvect = 2*atan( tan(Evect/2) * sqrt((1+e)/(1-e)) );

rangeRateData = nan(numPulsars,numObsvPerPulsar);
for j = 1:numPulsars
    thisPulsar = pulsarMat(:,j);
    for k = 1:numObsvPerPulsar
        thisParam = orbitParams;
        thisParam(6) = Fvect(j,k);
        [~,v] = Get_Orb_Vects(thisParam,mu);
        rangeRateData(j,k) = thisPulsar' * v;
    end
end

if randomEpoch
    Moffset = rand;
end
Tvect = (Mvect+Moffset) * sqrt(a^3/mu);

%% Verification Data

[rTrue,vTrue] = Get_Orb_Vects(orbitParams,mu);

%% Perturbations

if ~randomNoise
    rng(1)
end
if applyNoise
    noise = normrnd(0,noise,size(rangeRateData));
else
    noise = 0;
end
if noise == 0
    warning("Warning: No perturbations applied.")
end
rangeRateData_truth = rangeRateData;
rangeRateData = rangeRateData + noise;

%% Perform RROD

trueVars = [f0, e, dM];
[rEst,vEst,guessVars] = rrod(rangeRateData,Tvect,pulsarMat,mu);