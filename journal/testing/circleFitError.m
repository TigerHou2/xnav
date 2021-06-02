%% circleFitError.m
% This file compares the circle fitting error of measurements taken either
% at equal separations along the circle or at an uneven distribution with
% measurements becoming more sparse on one end. 
%
% The goal is to understand how having uneven distributions of data, which
% can happen when taking hodograph measurements along an eccentric orbit,
% can affect the accuracy of the hodograph determinatino.

%% setup
close all
clear;clc

c = [0,0];
R = 1;

N = 14:2:30;
DF = deg2rad(20:5:80);
itr = 1000;

err_eq = zeros(length(N),length(DF));
err_ne = zeros(length(N),length(DF));
err_eq_local = nan(length(N),length(DF));
err_ne_local = nan(length(N),length(DF));
X = nan(length(N),length(DF));
Y = nan(length(N),length(DF));

s = 0.01; % noise, standard deviation

direction = randn(max(N),2,itr);
unit_dir = direction ./ vecnorm(direction,2,2);
noise = unit_dir .* normrnd(0,s,max(N),1,itr);

%% simulate
for k = 1:itr
    for i = 1:length(N)
        n = N(i);
        for j = 1:length(DF)
            df = DF(j);

            %% contour plot x- and y-axes
            X(i,j) = n;
            Y(i,j) = rad2deg(df);

            %% perturbation
            noise_local = noise(1:n,:,k);

            %% Equispaced
            angs = linspace(0,df,n)';
            XY = R * [cos(angs),sin(angs)] + noise_local;
            [x,y] = hyperfit(XY);
            err_temp = norm(sqrt((x-XY(:,1)).^2+(y-XY(:,2)).^2)-R);
            err_eq_local(i,j) = err_temp;

            %% Non-Equispaced
            angs = linspace(0,sqrt(df),n)'.^2;
            XY = R * [cos(angs),sin(angs)] + noise_local;
            [x,y] = hyperfit(XY);
            err_temp = norm(sqrt((x-XY(:,1)).^2+(y-XY(:,2)).^2)-R);
            err_ne_local(i,j) = err_temp;

        end
    end
    err_eq = err_eq + err_eq_local/itr;
    err_ne = err_ne + err_ne_local/itr;
end

%% plot results
close all

figure(1)
contourf(X,Y,err_eq)
xlabel('Number of Measurements')
ylabel('$\Delta f$')
colorbar

figure(2)
contourf(X,Y,err_ne)
xlabel('Number of Measurements')
ylabel('$\Delta f$')
colorbar

disp('Done!')