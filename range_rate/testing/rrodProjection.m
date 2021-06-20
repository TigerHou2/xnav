%% RROD Projection Test
% This script looks at range-rate magnitude for various eccentricities, 
% true anomalies, and range-rate vectors. 
%
% Variables:
%   RA  = right ascension
%   IN  = inclination
%   f   = true anomaly
%   phi = angle between velocity vector and hodograph center vector c
%   gamma = angle between velocity vector and range-rate vector

% close all
% clear;clc

R = 1;
e = 0.7;
f = linspace(0,2*pi,2000);

RA = deg2rad(0);
IN = deg2rad(atan2d(10,1));

phi = acos( (e+cos(f)) ./ sqrt(1+e^2+2*e*cos(f)) );
v = R * sqrt(1+e^2+2*e*cos(f));
phi = phi .* (1-2*(f>=pi));

% cos_gamma = cos(phi-RA).*cos(-IN);
% rr = cos_gamma .* v;

A = R * cos(IN);
rr = A*e*cos(RA) + A*cos(f-RA);

% figure
hold on
plot(rad2deg(f),rr)
plot(rad2deg(f),v)
hold off
