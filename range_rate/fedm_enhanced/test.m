%% f-e-dM method test script
%
% Author: Tiger Hou
%
% Range-rate data measured in the direction of several stars are provied
% along with the time elapsed between measurements. By guessing the initial
% true anomaly, eccentricity, and the elapsed mean anomaly, range-rate data
% may be used to validate whether the guess was close to the true
% properties of the orbit.

close all
clear

addpath('../fcns_circ')
addpath('../fcns_orb')

%% Settings
applyNoise = false;
randomNoise = false;
randomEpoch = false;
randomPulsarRotation = false;

noise = 1;
Moffset = 0.7;
pulsarRotationMatrix = rotz(25) * rotx(12) * roty(2);

%% Define Orbit

mu = 3.986e14;
a = (6378+600)*1e3;
e = 0.93423;
i = deg2rad(45);
o = deg2rad(0);
w = deg2rad(0);
f0 = deg2rad(50.62);
orbitParams = [a, e, i, o, w, f0];
R = sqrt(mu/a/(1-e^2));

%% Define Pulsars

pulsars = [ [1, 0, 0]; ... pulsar 1
            [0, 1, 1]; ... pulsar 2
            [1, 1, 4]; ... pulsar 3
%             [-1, 1, -3]; ... pulsar 4
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
dM = 2*pi * 0.17235;
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

%% Test Objective Function using True Solution

trueVars = [f0, e, dM];
fVal = objFcn(trueVars, rangeRateData, Tvect, pulsarMat, mu);
disp(['True solution residual: ' num2str(sum(fVal)^2)])

%% Initial Guess and Search Radius

warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

initGuess = [pi, 0.5, pi];
maxRadius = [pi+0.1, 0.5, pi];
lsq_lb = [-0.1,0,0];
lsq_ub = [2*pi+0.1,0.999,2*pi];
radius = maxRadius;
fun_scalar = @(x) sum(objFcn(x,rangeRateData,Tvect,pulsarMat,mu));
fun_vector = @(x) objFcn(x,rangeRateData,Tvect,pulsarMat,mu);
options = optimoptions('fsolve','Display','none' ...
                               ,'MaxFunctionEvaluations',8000 ...
                               ,'StepTolerance', 1e-12 ...
                               ,'FunctionTolerance', 1e-12 ...
                               ,'MaxIterations',2000 ...
                               ,'Algorithm','levenberg-marquardt');
options_lsq = optimoptions('lsqnonlin','Display','none' ...
                               ,'MaxFunctionEvaluations',400 ...
                               ,'StepTolerance', 1e-32 ...
                               ,'FunctionTolerance', 1e-16 ...
                               ,'MaxIterations',100 ...
                               ,'Algorithm','trust-region-reflective');
options_fminsearch = optimset(  'Display','none', ...
                                'TolFun',1e-32, ...
                                'TolX', 1e-16);
iteration = 0;
residual = fun_scalar(initGuess)^2;
residualBest = residual;
% res = [24,16,24];
res = [10,10,10];

%% Guess and Search for Solution
% profile on

converged = 0;

while 1

iteration = iteration + 1;
disp('====')
disp(['Iteration ' num2str(iteration) ':'])

lb = initGuess - radius;
ub = initGuess + radius;
range_f0 = linspace(lb(1),ub(1),res(1)+1);
range_e  = linspace(lb(2),ub(2),res(2)+1);
range_dM = linspace(lb(3),ub(3),res(3)+1);
range_f0(range_f0<-0.1) = [];       % f0 >= -0.1
range_f0(range_f0>=2*pi+0.1) = [];  % f0 <=  0.1 + 2*pi
range_e(range_e<0) = [];            % ecc >= 0
range_e(range_e>0.999) = [];        % ecc < 0.999
range_dM(range_dM<1e-7) = [];       % dM >= 1e-7
thisRes = [length(range_f0),length(range_e),length(range_dM)];

[initGuess,fVal] = ...
        bounded_search(fun_scalar,range_f0,range_e,range_dM,thisRes);
residual = fVal^2;
if residual / residualBest > 0.95
    [initGuess,fVal] = fsolve(fun_scalar,initGuess,options);
%     [initGuess,fVal] = fminsearch(fun_scalar,initGuess,options_fminsearch);
    residual = fVal^2;
    disp('**fsolve used**')
else
    [initGuess,fVal] = ...
        fminsearch(fun_scalar,initGuess,options_fminsearch);
    residual = fVal^2;
    disp('**fminsearch used**')
end
if abs(residual-residualBest)/residualBest < 1e-5
    radius = radius * 0.8;
    converged = converged + 1;
elseif residual < residualBest
    residualBest = residual;
    radius = radius * 0.8^(1+iteration/20);
    converged = 0;
else
    converged = converged + 1;
end

initGuess(1) = mod(initGuess(1),2*pi);
disp(['    Guess:      ' num2str(initGuess)])
disp(['    True Soln.: ' num2str(trueVars)])
disp(['    Residual:   ' num2str(residual)])

if converged >= 2
    disp('Residual converged')
    break
end

if iteration >= 20
    disp('Iteration limit exceeded')
    break
end

if norm(radius) < 1e-4
    disp('Search radius converged')
    break
end

end
% disp('====')
% disp(['Iteration FINAL:'])
% initGuess = ...
%     lsqnonlin(fun_vector,initGuess,lsq_lb,lsq_ub,options_lsq);
% initGuess = ...
%     lsqnonlin(fun_vector,initGuess,lsq_lb,lsq_ub,options_lsq);
% initGuess(1) = mod(initGuess(1),2*pi);
% residual = fun_scalar(initGuess)^2;
% disp(['    Guess:      ' num2str(initGuess)])
% disp(['    True Soln.: ' num2str(trueVars)])
% disp(['    Residual:   ' num2str(residual)])
disp('Search complete')

warning('on','MATLAB:singularMatrix')
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:rankDeficientMatrix')

% profile viewer

%% Show fVal trend between guess and true solution
step = trueVars-initGuess;
mults = linspace(0,1,1000);
y = nan(size(mults));
for i = 1:length(y)
y(i) = fun_scalar(initGuess+mults(i)*step);
end
plot(mults,y)