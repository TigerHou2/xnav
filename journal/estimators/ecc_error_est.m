%% ecc_error_est.m
% 
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t eccentricity. Estimated
%   errors are generated using analytical equations and compared against
%   Monte Carlo simulations. 

%% initialization
close all
clear;clc
addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
addpath('..\fcns_misc')
savePath = 'plots\';
latexify

%% setup
eccVect = linspace(0.1,0.9,23);

mu = 1;
a = 1e5;
e = 0;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(134);

orbitParams = [a,e,i,o,w,f];

% check if initial true anomaly is 0, 135, or 180. 
% If it is none of the above, the figure is not saved.
save_plot = true;
ta_accept = deg2rad([0,135,180]);
if min(abs(f-ta_accept))>1e-8
    save_plot = false;
    warning(['The initial true anomaly provided is not required for ' ...
             'the manuscript, therefore it will not be saved.'])
end

% total duration spanned by all measurements, as fraction of orbit period
period = 0.1;
period = period * 2*pi;
% select the nth observation's position error for comparison
selObsv = 1;
% number of measurements
numObsv = 10;
% measurement noise
noise = 3e-6;
% Monte Carlo simulation size
numSims = 3000;
% resolution for the error model
model_res = 1000;

% line styles
MOD = 'r:'; % model
SIM = 'ko-.'; % simulation
    
% prepare measurement noise
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) * noise;

errDat = nan(numSims,length(eccVect));
errEst = nan(1,model_res);
v = nan(numObsv,3);

%% Reference Simulation
for i = 1:length(eccVect)
    % vary eccentricity
    e = eccVect(i);
    orbitParams(2) = e;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    E = kepler(M,e);
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
end

% create reference simulation x,y vectors
xRef = eccVect;
yRef = sqrt(mean(errDat.^2));

%% Error Model

xVar = linspace(min(eccVect),max(eccVect),model_res);

for i = 1:model_res
    
    e = xVar(i);
    
    % adj 1: error scales inversely with hodograph radius
    adj1 = sqrt((1-e)/(1+e)) + sqrt((1+e)/(1-e));
    adj1 = 1 / adj1;
    
    % adj 2: error scales inversely with square of measurement span
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    E = kepler(M,e);
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);
    df = fend-f;
    df = mod(df,2*pi);
    adj2 = 1 / df^2;
    
    errEst(i) = adj1 * adj2;
    
end

% create error model x,y vectors
yVar = errEst;

%% data processing & plotting
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
offset = 0;
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xRef,yRef,SIM,'LineWidth',1,'MarkerSize',5)
hold on
plot(xVar,yVar,MOD,'LineWidth',1.5,'MarkerSize',5)
hold off
% legend('Simulation','Prediction','Location','Best')
xlabel('Eccentricity')
ylabel('Position MSE, \%')
latexify(10,10,16)
setgrid
expand

if save_plot
    svnm = [savePath 'eccErr_f=' num2str(rad2deg(f))];
    print(svnm,'-dpdf','-bestfit')
end