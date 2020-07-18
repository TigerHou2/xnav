%% init_gmat.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Initializes the directory to include all supporting functions and set
%   the default template for figures.
%   This initialization file should only be used by functions designed to
%   interface with GMAT, which should be located in the \gmat directory.
%

%% initialization

% change pwd to location of init_gmat.m
cd 'E:\Research\XNAV\journal\gmat'
addpath(['..' filesep 'fcns_od'])
addpath(['..' filesep 'fcns_orb'])
addpath(['..' filesep 'fcns_vis'])
addpath(['..' filesep 'fcns_misc'])
addpath(['..' filesep 'figures'])
addpath(['..' filesep 'cases'])

close all hidden
clear;clc

latexify(0)