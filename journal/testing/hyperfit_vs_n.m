%% hyperfit_vs_n.m
%
% Description:
%   Compares the accuracy of the radius estimate for hyperfit against the
%   number of sample points used in the fit.

close all
clear;clc

addpath('..\fcns_misc')

span = deg2rad(50);
R = 1;
c = [0,0];

n = 5:3:100;

numSims = 100000;
R_est_vect = nan(length(n),numSims);

noise_cube = randn(max(n),2,numSims) * 0.003;

for i = 1:length(n)
    angs = linspace(0,span,n(i))';
    xy = R * [cos(angs),sin(angs)] + c;
    for s = 1:numSims
        noise = noise_cube(1:n(i),:,s);
        [~,~,R_est_vect(i,s)] = hyperfit_cpp(xy + noise);
    end
end

plot(n,mean(R_est_vect,2))
xlabel('Number of Samples')
ylabel('R Estimate Error')