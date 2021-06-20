%% compareVIOD.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares the performance of two different VIOD algorithms.

%% initialization
close all
clear;clc
addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
addpath('..\fcns_misc')
latexify

%% setup

smaVect = linspace(1e4,1e6,15);
eccVect = linspace(0.5,0.99,25);
taVect  = deg2rad(linspace(150,185,25));

varVect = taVect;

if isequal(varVect,smaVect)
    idx = 1;
    xname = 'Semi-Major Axis, DU';
elseif isequal(varVect,eccVect)
    idx = 2;
    xname = 'Eccentricity';
elseif isequal(varVect,taVect)
    idx = 6;
    xname = '$f_0$, rad';
end

mu = 1;
a = 1e5;
e = 0.9356;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(180);
    
orbitParams = [a,e,i,o,w,f];

% total duration spanned by all measurements, as fraction of orbit period
period = 0.06;
period = period * 2*pi;
% measurement noise
noise = 4e-8;
% number of measurements
numObsv = 10;
% Monte Carlo simulation size
numSims = 1000;
% select the nth observation's position error for comparison
selObsv = 1;

% line styles
S1 = 'rx:'; % style 1
S2 = 'ko-.'; % style 2

% prepare measurement noise
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,noise,1,1,numSims);

e1 = nan(numSims,length(smaVect));
e2 = nan(numSims,length(smaVect));
v = nan(numObsv,3);

%% simulation
for i = 1:length(varVect)
    % vary semi-major axis length
    var = varVect(i);
    orbitParams(idx) = var;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-orbitParams(2))/(1+orbitParams(2)))*tan(orbitParams(6)/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+orbitParams(2))/(1-orbitParams(2)))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);
    
    df = rad2deg(mod(fend-orbitParams(6),2*pi));
    disp(['df = ' num2str(df) ' deg']);

    % ground truth data
    tempParams = orbitParams;
    for j = 1:numObsv
        tempParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(tempParams,mu);
    end
    [rRef,R] = hodoHyp(v,mu);
    rRef = rRef(selObsv,:)';
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(:,:,s);
        
        [r1,R1] = hodoHyp(v+nvect,mu);
        r1 = r1(selObsv,:)';
        [r2,R2] = hodoHyp_debug(v+nvect,mu);
        r2 = r2(selObsv,:)';
        
        if (R1 > R*2) || (R1 < R/2) || (R2 > R*2) || (R2 < R/2)
            continue
        end
        
        e1(s,i) = norm(r1-rRef) / norm(rRef) * 100;
        e2(s,i) = norm(r2-rRef) / norm(rRef) * 100;
        
    end
    
end

if max(sum(isnan(e1))) > 0
    disp([  'Removed ' num2str(max(sum(isnan(e1)))) ...
            ' extreme cases with large errors.'])
end

%% data processing & plotting
err1 = sqrt(mean(e1.^2,'omitnan'));
err2 = sqrt(mean(e2.^2,'omitnan'));

figure;
plot(varVect,err1,S1,'LineWidth',1.3,'MarkerSize',6)
hold on
plot(varVect,err2,S2,'LineWidth',1.3,'MarkerSize',7)
hold off
legend('RMSE$_{old}$','RMSE$_{new}$','Location','Best')
xlabel(xname)
ylabel('RMSE, \%')
latexify(10,8,15)
setgrid
expand