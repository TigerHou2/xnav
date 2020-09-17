%% dur_err_est.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t measurement duration. 
%   Estimated errors are generated using analytical equations and compared
%   against Monte Carlo simulations. 

%% initialization
close all
clear;clc
addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
latexify

%% setup
durVect = linspace(0.05,0.9,20);
durVect = durVect * 2*pi;

mu = 1;
a = 1.56e5;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);

a = 1.29e4;
e = 0.1;
    
orbitParams = [a,e,i,o,w,f];

% measurement noise
noise = 1e-6;
% number of measurements
numObsv = 3;
% Monte Carlo simulation size
numSims = 3000;
% select the nth observation's position error for comparison
selObsv = 1;

% prepare measurement noise
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) * noise;

errDat = nan(numSims,length(durVect));
errEst = nan(1,length(durVect));
v = nan(numObsv,3);

%% simulation
for i = 1:length(durVect)
    % vary measurement duration
    period = durVect(i);
    
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
    
    % prepare SVD error estimator
    svdErrVect = nan(numSims,1);
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(:,:,s);
        r = hodo(v+nvect,mu);
        r = r(selObsv,:)';
        errDat(s,i) = norm(r-rRef) / norm(rRef) * 100;
        
        % test svd
        [~,~,V] = svd(v+nvect,0);
        k = V(:,end);
        if k'*[0;0;1]<0
            k = -k;
        end
        svdErrVect(s) = norm(k-[0;0;1]);
    end
    
    % estimate error
    df = fend-f;
    df = mod(df,2*pi);
    
        % adj 1: error scales with orbit normal error
        adj1 = mean(svdErrVect);
        
        % adj 2: error scales inversely with square of measurement span
        adj2 = 1 / df^2;
    
    errEst(i) = adj1 * adj2;
    
end

%% data processing & plotting
xVar = durVect /2/pi * 100;
yRef = mean(errDat);
yVar = errEst;
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xVar,yVar,'LineWidth',1.5)
hold on
plot(xVar,yRef,'LineWidth',1.5)
hold off
legend('Prediction','Simulation','Location','Best')
xlabel('Measurement Duration, \% of Orbit Period')
ylabel('Error Avg. \%')
set(gca,'FontSize',18)
grid on
latexify(16)