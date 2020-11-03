%% circ_to_viod.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares VIOD error trend with the error in the estimate of the
%   hodograph radius. This is the key assumption upon which the entire
%   paper is built. We should see matching trends between VIOD error and
%   hodograph radius error w.r.t. several variables.
%
% Details:
%   This demo obviously cannot cover all scenarios. What we will try to do
%   instead is to establish a baseline case, then vary each parameter
%   individually.

%% initialization
close all
clear;clc

addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
addpath('..\fcns_misc')
savePath = 'plots\';
latexify

%% baseline case

% mission parameters
mu = 1;
a = 1e5; % param
e = 0.5; % param
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(90); % param
orbitParams = [a,e,i,o,w,f];

% sensor parameters
period = 0.1 * (2*pi); % param
numObsv = 10; % param
noise = 3e-6; % param

% simulation parameters
selObsv = 1;
numSims = 3000;

% line styles
CIRC = 'rx:'; % error in circle fitting, from hyperfit
VIOD = 'ko-.'; % error in position MSE, from VIOD

%% parameter variations

% mission parameters
smaVect = linspace(1e4,1e6,15);
eccVect = linspace(0.1,0.9,23);
taVect = deg2rad(linspace(0,350,36));

% sensor parameters
noiseVect = linspace(5e-7,1e-5,10);
durVect = linspace(0.05,0.9,20);
obsvVect = 3:5:100;

%% semi-major axis

err_circ = nan(numSims,length(smaVect));
err_viod = nan(numSims,length(smaVect));
params = orbitParams;

for i = 1:length(smaVect)
    
    a = smaVect(i);
    params(1) = a;
    
    
end

%% eccentricity

err_circ = nan(numSims,length(eccVect));
err_viod = nan(numSims,length(eccVect));
params = orbitParams;

for i = 1:length(eccVect)
    
end

%% initial true anomaly

err_circ = nan(numSims,length(taVect));
err_viod = nan(numSims,length(taVect));
params = orbitParams;

for i = 1:length(taVect)
    
end

%% noise

err_circ = nan(numSims,length(noiseVect));
err_viod = nan(numSims,length(noiseVect));
params = orbitParams;

for i = 1:length(noiseVect)
    
end

%% measurement interval

err_circ = nan(numSims,length(durVect));
err_viod = nan(numSims,length(durVect));
params = orbitParams;

for i = 1:length(durVect)
    
end

%% number of measurements

err_circ = nan(numSims,length(obsvVect));
err_viod = nan(numSims,length(obsvVect));
params = orbitParams;

for i = 1:length(obsvVect)
    
end