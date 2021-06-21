%% f-e-dM method test script
%
% Author: Tiger Hou
%
% Range-rate data measured in the direction of several stars are provied
% along with the time elapsed between measurements. By guessing the initial
% true anomaly, eccentricity, and the elapsed mean anomaly, range-rate data
% may be used to validate whether the guess was close to the true
% properties of the orbit.
%
%

close all
clear

addpath('../fcns_circ')
addpath('../fcns_orb')

%% Define Orbit

mu = 3.986e14;
a = (6378+600)*1e3;
e = 0.53423;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f0 = deg2rad(163.62);
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
numPulsars = size(pulsars,1); % counter the number of pulsars

%% Simulate Measurements

groupMeasurements = true;  % measurements are taken in groups.
                            % If true, each group would contain
                            % all measurements from only one pulsar.
                            % If false, each group would contain one
                            % measurement from each pulsar.
% time interval ratios
memberInterval = 1;     % after each measurement within the same group
groupInterval = 0;      % between groups

numObsvPerPulsar = 4;
dM = 2*pi * 0.1235;
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

Toffset = rand;
Tvect = Mvect * sqrt(a^3/mu) + Toffset;

%% Perturbations

noise = normrnd(0,R/1000,size(rangeRateData));
rangeRateData_truth = rangeRateData;
rangeRateData = rangeRateData + noise;

%% Test Objective Function using True Solution

trueVars = [f0, e, dM];
fVal = objFcn(trueVars, rangeRateData, Tvect, pulsarMat, mu);
disp(['True solution residual: ' num2str(fVal^2)])

%% Initial Guess and Search Radius

initGuess = [pi, 0.5, pi];
radius = [pi, 0.5, pi];
fun = @(x) objFcn(x,rangeRateData,Tvect,pulsarMat,mu);
options = optimoptions('fsolve','Display','none' ...
                               ,'MaxFunctionEvaluations',12000 ...
                               ,'StepTolerance', 1e-16 ...
                               ,'FunctionTolerance', 1e-16 ...
                               ,'MaxIterations',3000 ...
                               ,'Algorithm','levenberg-marquardt');
iteration = 0;
residual = fun(initGuess);
residualPrev = 0;

%% Guess and Search for Solution
% profile on
while 1

iteration = iteration + 1;
disp('====')
disp(['Iteration ' num2str(iteration) ':'])

lb = initGuess - radius;
ub = initGuess + radius;
lb(1) = max(lb(1),0);       % f0 >= 0
ub(1) = min(ub(1),2*pi);    % f0 <= 2*pi
lb(2) = max(lb(2),0);       % eccentricity >= 0
ub(2) = min(ub(2),0.999);   % eccentricity <= 0.999
lb(3) = max(lb(3),1e-7);    % dM >= 1e-7
res = [20,20,20];
range_f0 = linspace(lb(1),ub(1),res(1)+1);
range_e  = linspace(lb(2),ub(2),res(2)+1);
range_dM = linspace(lb(3),ub(3),res(3)+1);

initGuess = bounded_search(fun,range_f0,range_e,range_dM,res);
if abs(residual-residualPrev) / residualPrev < 5e-2
    initGuess = fsolve(fun,initGuess,options);
end
residualPrev = residual;
thisResidual = fun(initGuess)^2;
if thisResidual <= residual
    radius = radius/2;
    residual = thisResidual;
else
    radius = radius * 0.9;
end
disp(['    Guess:      ' num2str(initGuess)])
disp(['    True Soln.: ' num2str(trueVars)])
disp(['    Residual:   ' num2str(thisResidual)])

if iteration >= 30
    break
end

end
% profile viewer
