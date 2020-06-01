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

addpath(['.' filesep 'od']);
addpath(['.' filesep 'misc']);