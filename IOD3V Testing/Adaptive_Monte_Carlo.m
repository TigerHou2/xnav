%% Adaptive_Monte_Carlo.m
%
% Author: Tiger Hou
%
% Monte Carlo simulation of two variants of the adaptive meshing IOD
% method, along with comparison to IOD based on equispaced time.


%% Initial Conditions

% clear workspace and load data
close all hidden
clear;clc

Parameters

% choose functions to run
run_predictive = false;
run_observed   = true;
run_analytic   = false;
run_equispace  = true;

% choose display options
show_orbits = false;

% choose code profiling options
do_profile = false;

% spacecraft conditions

r1 = earth.a; % departure planet
r2 = neptune.a; % arrival planet
mu = sun.mu; % gravitational parameter [km^3/s^2]

theta = pi/4;
phi   = pi/4;
r = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];

r = r * r1 * AU;
v = cross([0,0,1],r);
v = v / norm(v) * sqrt(mu*(2/norm(r)-2/(r1+r2)/AU));

noise     = 10;             % 1 sigma sensor noise [m/s]
start_day = 0.001;            % measurement start day [days]

% observation conditions

obsv_cap = 80; % max number of observation allowed
da = deg2rad(0.5);

% seed for random number generator

rng_val = 6;
rng(rng_val);

% plotting

ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit

% Monte Carlo settings

samp_size = 10;

% symbolic math preallocation

syms x

% True Position Reference

time_ref = dur/2; % apoapsis
ap_ref = TimeProp_V3(r,v,mu,time_ref,x);
v0 = TimeProp_V3(r,v,mu,start_day,x);

% code profiling

if do_profile
    profile on
end


%% Simulation (Adaptive v1, predictive time propagation)

if run_predictive

    disp('Running predictive sim...')
    tic

    errs_pred = zeros(3,samp_size);
    time_pred = zeros(samp_size,1);

    for i = 1:samp_size
        [R_est,R,V,T] = viod_adap_pred(r,v,mu,start_day,obsv_cap,da,noise,x);
        ap_pred = TimeProp_V3(R_est(end,:),V(end,:),mu,time_ref-T(end),x);
        errs_pred(:,i) = (ap_pred - ap_ref)';
        time_pred(i) = T(end)-T(1);
        
        if show_orbits
            figure(1)
            hold on
            disp_orbit(R_est,V,mu,res,dur,start_day)
            disp_obsv(R)
            hold off
        end
    end
    
    if show_orbits
        hold on
        disp_orbit(r,v,mu,res,dur,start_day,1)
        hold off
    end

    toc
    disp('Predictive sim complete!')
    disp(' ')
    
end


%% Simulation (Adaptive v2, observed time propagation)

if run_observed

    disp('Running observed sim...')
    tic

    errs_obsv = zeros(3,samp_size);
    time_obsv = zeros(samp_size,1);

    end_day = zeros(samp_size,1);

    for i = 1:samp_size
        [R_est,R,V,T] = viod_adap_obsv(r,v,mu,start_day,obsv_cap,da,noise,x);
    %     [R_est,R,V,T] = viod_adap_obsv(r,v,mu,start_day,obsv_cap,delta_ang(i),noise);
        ap_obsv = TimeProp_V3(R_est(end,:),V(end,:),mu,time_ref-T(end),x);
        errs_obsv(:,i) = (ap_obsv - ap_ref)';
        time_obsv(i) = T(end)-T(1);

        end_day(i) = T(end);
        
        if show_orbits
            figure(2)
            hold on
            disp_orbit(R_est,V,mu,res,dur,start_day)
            disp_obsv(R)
            hold off
        end
    end
    
    if show_orbits
        hold on
        disp_orbit(r,v,mu,res,dur,start_day,1)
        hold off
    end

    toc
    disp('Observed sim complete!')
    disp(' ')
    
end


%% Simulation (Adaptive v3, analytic time propagation)

if run_analytic

    disp('Running analytic sim...')
    tic

    errs_anlt = zeros(3,samp_size);
    time_anlt = zeros(samp_size,1);

    for i = 1:samp_size
        [R_est,R,V,T] = viod_adap_analytic(r,v,mu,start_day,obsv_cap,da,noise,x);
        ap_anlt = TimeProp_V3(R_est(end,:),V(end,:),mu,time_ref-T(end),x);
        errs_anlt(:,i) = (ap_anlt - ap_ref)';
        time_anlt(i) = T(end)-T(1);
        
        if show_orbits
            figure(3)
            hold on
            disp_orbit(R_est,V,mu,res,dur,start_day)
            disp_obsv(R)
            hold off
        end
    end
    
    if show_orbits
        hold on
        disp_orbit(r,v,mu,res,dur,start_day,1)
        hold off
    end

    toc
    disp('Analytic sim complete!')
    disp(' ')
    
end


%% Simulation (Equispaced)

if run_equispace

    disp('Running equispaced sim...')
    tic

    errs_eqsp = zeros(3,samp_size);
    time_eqsp = zeros(samp_size,1);

    for i = 1:samp_size
        [R_est,R,V,T] = viod_eqsp(r,v,mu,start_day,end_day(i),obsv_cap,noise,x);
        ap_eqsp = TimeProp_V3(R_est(end,:),V(end,:),mu,time_ref-T(end),x);
        errs_eqsp(:,i) = (ap_eqsp - ap_ref)';
        time_eqsp(i) = T(end)-T(1);
        
        if show_orbits
            figure(4)
            hold on
            disp_orbit(R_est,V,mu,res,dur,start_day)
            disp_obsv(R)
            hold off
        end
    end
    
    if show_orbits
        hold on
        disp_orbit(r,v,mu,res,dur,start_day,1)
        hold off
    end

    toc
    disp('Equispaced sim complete!')
    disp(' ')
    
end


%% Result Comparison

if run_predictive
    errs_norm_pred = vecnorm(errs_pred);
    disp(['======== Predictive Adaptive VIOD ========' newline ...
          '    Error Mean: ' num2str(mean(errs_norm_pred)) ' km' ...
          ' (' num2str(mean(errs_norm_pred)/AU) ' AU)' newline ...
          '   Error Stdev: ' num2str(std(errs_norm_pred)) ' km' ...
          ' (' num2str(std(errs_norm_pred)/AU) ' AU)' newline ...
          ' Avg Obsv Time: ' num2str(mean(time_pred)) ' days' newline])
end

if run_observed
    errs_norm_obsv = vecnorm(errs_obsv);
    disp(['======== Observed Adaptive VIOD ========' newline ...
          '    Error Mean: ' num2str(mean(errs_norm_obsv)) ' km' ...
          ' (' num2str(mean(errs_norm_obsv)/AU) ' AU)' newline ...
          '   Error Stdev: ' num2str(std(errs_norm_obsv)) ' km' ...
          ' (' num2str(std(errs_norm_obsv)/AU) ' AU)' newline ...
          ' Avg Obsv Time: ' num2str(mean(time_obsv)) ' days' newline])
end
  
if run_analytic
    errs_norm_anlt = vecnorm(errs_anlt);
    disp(['======== Analytic Adaptive VIOD ========' newline ...
          '    Error Mean: ' num2str(mean(errs_norm_anlt)) ' km' ...
          ' (' num2str(mean(errs_norm_anlt)/AU) ' AU)' newline ...
          '   Error Stdev: ' num2str(std(errs_norm_anlt)) ' km' ...
          ' (' num2str(std(errs_norm_anlt)/AU) ' AU)' newline ...
          ' Avg Obsv Time: ' num2str(mean(time_anlt)) ' days' newline])
end
  
if run_equispace
    errs_norm_eqsp = vecnorm(errs_eqsp);
    disp(['======== Equispace VIOD ========' newline ...
          '    Error Mean: ' num2str(mean(errs_norm_eqsp)) ' km' ...
          ' (' num2str(mean(errs_norm_eqsp)/AU) ' AU)' newline ...
          '   Error Stdev: ' num2str(std(errs_norm_eqsp)) ' km' ...
          ' (' num2str(std(errs_norm_eqsp)/AU) ' AU)' newline ...
          ' Avg Obsv Time: ' num2str(mean(time_eqsp)) ' days' newline])
end

if do_profile
    profile viewer
end