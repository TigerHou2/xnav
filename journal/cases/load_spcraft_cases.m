%% load_spcraft_cases.m
function [noiseVect,durVect,obsVect] = load_spcraft_cases()
% Author:
%   Tiger Hou
%
% Description:
%   This file contains spacecraft parameter vectors that are used for case
%   studies for velocity orbit determination error. This file should be
%   called by plotSpcraftParams.m after initialization.

%% initialize cases

% noise values in canonical units
noiseVect = [1e-7,5e-7,1e-6,5e-6,1e-5]; % DU/TU

% measurement duration as fraction of the orbit period
durVect = [0.01 0.04 0.8]; % nd

% measurement count
obsVect = [3 10 40]; % nd

end % load_spcraft_cases.m