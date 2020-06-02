%% init.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Clears the MATLAB workspace and closes all windows. Adds relevant
%   file paths. Designed to be called from GMAT.

close all hidden
clear;clc

% change pwd to location of init.m
cd 'C:\Users\tigre\Desktop\Directory\Research\xnav\viod_gmat'

addpath(['.' filesep 'od']);
addpath(['.' filesep 'misc']);
addpath(['.' filesep 'scripts'])