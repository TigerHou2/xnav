%% load_cases_orbit.m
function [smaVect,eccVect,taVect,noise_canon,names] = load_orbit_cases()
% Author:
%   Tiger Hou
%
% Description:
%   This file contains orbital parameter vectors that are used for case
%   studies for velocity orbit determination error. This file should be
%   called by plotOrbitParams.m after initialization.

%% define noise, initialize cases

% we wish to set up canonical units such that
% noise is 1e-6 DU/TU
noise_canon = 1e-6;

noise = 3; % m/s
noise = noise / 1000; % km/s

eccVect = [0.1,0.5,0.9];
taVect  = [0, 15, 30, 45, 90, 120];
taVect  = deg2rad(taVect);
smaVect = [];
names = {};

%% LEO case
mu  = 3.986e5; % km^3/s^2
SMA = 7500; % km
smaVect(end+1) = SMA / (mu/noise^2) / noise_canon^2;
names{end+1} = 'LEO';

%% Earth-Neptune transfer case
mu = 1.327e11; % km^3/s^2
SMA = 2.3e9; % km
smaVect(end+1) = SMA / (mu/noise^2) / noise_canon^2;
names{end+1} = 'Earth-Neptune Transfer';

%% Earth-Mars transfer case
mu = 1.327e11; % km^3/s^2
SMA = 1.9e8; % km
smaVect(end+1) = SMA / (mu/noise^2) / noise_canon^2;
names{end+1} = 'Earth-Mars Transfer';

end % load_orbit_cases.m