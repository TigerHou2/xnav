%% circ2viod_dur.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares VIOD error trend with the error in the estimate of the
%   hodograph radius. This is the key assumption upon which the entire
%   paper is built. We should see matching trends between VIOD error and
%   hodograph radius error w.r.t. the total measurement time interval.

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
durVect = linspace(0.05,0.5,20);
durVect = durVect * 2*pi;

mu = 1;
a = 1e5;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);
    
orbitParams = [a,e,i,o,w,f];

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

errDat = nan(numSims,length(durVect));
errCirc = nan(numSims,length(durVect));
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
    [~,R] = hodoHyp_debug(v,mu);
    
    % get position reference
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(:,:,s);
        [r,Rest] = hodoHyp_debug(v+nvect,mu);
        r = r(selObsv,:)';
        errDat(s,i)  = norm(r-rRef) / norm(rRef) * 100;
        errCirc(s,i) = abs(R-Rest) / abs(R) * 100;
        
        % test svd
        [~,~,V] = svd(v+nvect,0);
        k = V(:,end);
        if k'*[0;0;1]<0
            k = -k;
        end
        
        errCirc(s,i) = errCirc(s,i) * norm(k-[0;0;1])^2;
    end
end

xRef = durVect / (2*pi) * 100;
yRef = sqrt(mean(errDat.^2));

xVar = durVect / (2*pi) * 100;
yVar = sqrt(mean(errCirc.^2));

%% data processing & plotting
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xRef,yRef,SIM,'LineWidth',1,'MarkerSize',5)
hold on
plot(xVar,yVar,MOD,'LineWidth',1,'MarkerSize',5)
hold off
% legend('VIOD','Hodograph','Location','Best')
xlabel('\% of Orbit Period')
ylabel('Mean Squared Error, \%')
latexify(10,8,16)
setgrid
expand
svnm = [savePath 'durProof'];
% print(svnm,'-dpdf','-bestfit')