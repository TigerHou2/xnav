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

mu = 1;
a = 1;
e = 0.9489;
i = deg2rad(45);
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
          ];
pulsarMat = pulsars' ./ vecnorm(pulsars');
numPulsars = size(pulsars,1);

%% Simulate Measurements

groupMeasurements = false;

numObsvPerPulsar = 5;
dM = 2*pi * 0.3235;
E0 = 2*atan( tan(f0/2) * sqrt((1-e)/(1+e)) );
M0 = E0 - e*sin(E0);
Mf = M0 + dM;
Mvect = linspace(M0,Mf,numPulsars*numObsvPerPulsar);
if groupMeasurements
    Mvect = reshape(Mvect,numObsvPerPulsar,numPulsars)';
else
    Mvect = reshape(Mvect,numPulsars,numObsvPerPulsar);
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

%% Test Objective Function using True Solution

trueVars = [f0, e, dM];
fVal = objFcn(trueVars, rangeRateData, Tvect, pulsarMat, mu);
disp(['True solution error: ' num2str(fVal)])

%% Initial Guess and Search Radius

initGuess = [pi, 0.5, pi];
radius = [pi, 0.5, pi];
fun = @(x) objFcn(x,rangeRateData,Tvect,pulsarMat,mu);
options = optimoptions('fsolve','Display','iter' ...
                               ,'MaxFunctionEvaluations',12000 ...
                               ,'StepTolerance', 1e-16 ...
                               ,'FunctionTolerance', 1e-16 ...
                               ,'OptimalityTolerance', 1e-16 ...
                               ,'MaxIterations',3000 ...
                               ,'Algorithm','levenberg-marquardt');
iteration = 0;

%% Guess and Search for Solution

iteration = iteration + 1;
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
initGuess = fsolve(fun,initGuess,options);
if fun(initGuess) < 1e-3
    radius = radius/2;
end
disp(['Guess:      ' num2str(initGuess)])
disp(['True Soln.: ' num2str(trueVars)])