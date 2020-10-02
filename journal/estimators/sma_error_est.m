%% sma_error_est.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t semi-major axis length. 
%   Estimated errors are generated using analytical equations and compared
%   against Monte Carlo simulations. 

%% initialization
close all
clear;clc
addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
savePath = 'plots\';
latexify

%% setup
smaVect = linspace(1e4,1e6,15);

mu = 1;
a = 0;
e = 0.5;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(90);
    
orbitParams = [a,e,i,o,w,f];

% total duration spanned by all measurements, as fraction of orbit period
period = 0.1;
period = period * 2*pi;
% measurement noise
noise = 3e-6;
% number of measurements
numObsv = 10;
% Monte Carlo simulation size
numSims = 3000;
% select the nth observation's position error for comparison
selObsv = 1;

% line styles
MOD = 'rx:'; % model
SIM = 'ko-.'; % simulation

% prepare measurement noise
nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* nGauss;

errDat = nan(numSims,length(smaVect));
errEst = nan(1,length(smaVect));
v = nan(numObsv,3);

%% simulation
for i = 1:length(smaVect)
    % vary semi-major axis length
    sma = smaVect(i);
    orbitParams(1) = sma;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);

    % ground truth data
    for j = 1:numObsv
        orbitParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    
    % get position reference
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(:,:,s);
        r = hodoHyp(v+nvect,mu);
        r = r(selObsv,:)';
        errDat(s,i) = norm(r-rRef) / norm(rRef) * 100;
    end
    
    % estimate error

        % adj 1: error scales with sqrt(sma)
        adj1 = sqrt(sma);
    
    errEst(i) = adj1;
    
end

%% data processing & plotting
xVar = smaVect;
yRef = sqrt(mean(errDat.^2));
yVar = errEst;
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xVar,yRef,SIM,'LineWidth',1.5,'MarkerSize',7)
hold on
plot(xVar,yVar,MOD,'LineWidth',1.5,'MarkerSize',7)
hold off
legend('Simulation','Prediction','Location','Best')
xlabel('Semi-Major Axis Length, DU')
ylabel('Position MSE, \%')
latexify(20,13,18)
setgrid
expand
svnm = [savePath 'smaErr'];
print(svnm,'-dpdf','-bestfit')