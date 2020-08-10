%% plotEarth.m
function plotEarth(numSims,rngSeed,noise)
%
% Author:
%   Tiger Hou
%
% Description:
%   Generates plots for velocity initial orbit determination (VIOD) error
%   in Earth orbit with perturbations from solar radiation pressure, Earth
%   oblateness, and the gravitational pull of the Moon and Sun.
%   Propagation is performed in GMAT.
%
% Arguments:
%   numSims - number of Monte Carlo simulations to run
%   rngSeed - seed for the randon number generator
%   noise   - measurement noise in m/s

%% collect ground truth data from .mat file and generate perturbations

if ~exist('numSims','var')
    numSims = 500;
    rngSeed = 1;
    noise = 3;
    disp('Warning: No inputs provided. Using default sim settings.')
end

% load cell matrix of position, velocity, and mu data
load('temp\datEarth.mat','rArray','vArray','muArray');

% get velocity data size, i.e. number of observations and dimensionality
vSize = size(vArray{1});

% initialize cell array of noise values
nArray = cell(1,numSims);
% populate noise cell array
rng(rngSeed);
for n = 1:numSims
    noiseVect = randn(vSize);
    noiseVect = noiseVect ./ vecnorm(noiseVect,2,2) * noise / 1000;
    nArray{n} = noiseVect;
end %numSims

% data logging
filepath = 'data\datEarth.mat';

%% apply noise and perform VIOD
% we want to look at several types of perturbations:
%   - measurement noise only
%   - moon + sun gravity
%   - solar radiation pressure
%   - earth oblateness
%   - everything except noise, combined
%   - everything including noise

errVect = nan(size(vArray,2),size(vArray,1)+1);

for i = 1:size(vArray,2) % iterate over LEO and GEO cases
    idx2 = i;
for j = 1:size(vArray,1)+1 % iterate over diff. perturbations
    idx1 = j;
    if idx1 > size(vArray,1)
        idx1 = size(vArray,1);
    end
    R = rArray{idx1,idx2};
    V = vArray{idx1,idx2};
    mu = muArray{idx1,idx2};
    
    if j == 1 || j == size(vArray,1)+1 % need Monte Carlo for noise
        err = nan(1,numSims);
        for s = 1:numSims
            noiseVect = nArray{s};
            r = hodo(V+noiseVect,mu);
            r = r(1,:);
            err(s) = norm(R-r)/norm(R)*100;
        end %s
        err = mean(err);
    else
        r = hodo(V,mu);
        r = r(1,:);
        err = norm(R-r)/norm(R)*100;
    end
    errVect(i,j) = err;
end %j
end %i

list = {'Noise Only','Moon \& Sun','Solar Rad. Pressure','JGM-2 4$\times$4',...
        'All Perturbations','Perturbations \& Noise'};
x = categorical(list,list); % preserve list oder
for i = 1:size(vArray,2)
    figure(i)
    bar(x,errVect(i,:))
%     set(gca,'YScale','log')
    grid on
    ylabel('Position Error \%')
    set(gca,'FontSize',16)
end

end %plotEarth.m