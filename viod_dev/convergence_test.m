%% convergence_test.m
%
% Author: Tiger Hou
%
% Description:
%   This script tests the convergence speed of various velocity IOD methods
%   by varying either the number of measurements or the measurement noise.
%
% Two algorithms will be compared: 
%   - equal time (mean anomaly change)
%   - equal turn (velocity vector change)
%
%
%% Initial Conditions

% clear workspace and load data
close all hidden
clear;clc

Parameters

% choose error compare method (true = use all vel, false = use final vel)
do_avg_err = true;

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
start_day = 0.001;          % measurement start day [days]

% observation conditions

obsv_ang = deg2rad(60);
obsv_cap = 20; % max number of observation allowed
da = obsv_ang / (obsv_cap-1);

% seed for random number generator

rng_val = 5;
rng(rng_val);

% plotting

ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit

% Monte Carlo settings

samp_size = 40;

% True Position Reference

time_ref = dur/2; % apoapsis
ap_ref = TimeProp_V4(r,v,mu,time_ref);
v0 = TimeProp_V4(r,v,mu,start_day);


%% Measurement Count

count_vect = [3,5,10,20,40];
da_vect = obsv_ang ./ (count_vect-1);

% ======== Equal Turn ========

errs_obsv = zeros(3,length(da_vect));
time_obsv = zeros(length(da_vect),1);

end_day = zeros(length(da_vect),1);

for i = 1:length(da_vect)
    
    for ss = 1:samp_size
    
        [R_est,R,V,T] = ...
            viod_adap_obsv(r,v,mu,start_day,count_vect(i),da_vect(i),noise);
        if do_avg_err
            ap_obsv = zeros(1,3);
            for j = 1:count_vect(i)
                ap_obsv = ap_obsv + ...
                          TimeProp_V4(R_est(j,:),V(j,:),mu,time_ref-T(j));
            end
            ap_obsv = ap_obsv / count_vect(i);
        else
            ap_obsv = TimeProp_V4(R_est(end,:),V(end,:),mu,time_ref-T(end));
        end
        
    end
    errs_obsv(:,i) = errs_obsv(:,i) + (ap_obsv - ap_ref)';
    time_obsv(i) = T(end)-T(1);

    end_day(i) = T(end);
    
end

% ======== Equal Time ========

errs_eqsp = zeros(3,length(da_vect));
time_eqsp = zeros(length(da_vect),1);

for i = 1:length(da_vect)
    
    for ss = 1:samp_size
        
        [R_est,R,V,T] = ...
            viod_eqsp(r,v,mu,start_day,end_day(i),count_vect(i),noise);
        if do_avg_err
            ap_eqsp = zeros(1,3);
            for j = 1:count_vect(i)
                ap_eqsp = ap_eqsp + ...
                          TimeProp_V4(R_est(j,:),V(j,:),mu,time_ref-T(j));
            end
            ap_eqsp = ap_eqsp / count_vect(i);
        else
            ap_eqsp = TimeProp_V4(R_est(end,:),V(end,:),mu,time_ref-T(end));
        end
        
    end
    
    errs_eqsp(:,i) = errs_eqsp(:,i) + (ap_eqsp - ap_ref)';
    time_eqsp(i) = T(end)-T(1);
    
end

% ======== Error Plot ========

figure(1)
loglog(count_vect,vecnorm(errs_obsv/samp_size,2))
hold on
loglog(count_vect,vecnorm(errs_eqsp/samp_size,2))
hold off
grid(gca,'minor')
grid on
legend('Equal Turn','Equal Time')
xlabel('Number of Measurements')
ylabel('Error at Apoapsis (km)')
 title('Error Convergence Rate')


%% Measurement Noise

noise_vect = [1e-3,1e-2,0.1,1,5,10,30];

% ======== Equal Turn ========

errs_obsv = zeros(3,length(noise_vect));
time_obsv = zeros(length(noise_vect),1);

end_day = zeros(length(noise_vect),1);

for i = 1:length(noise_vect)
    
    for ss = 1:samp_size
    
        [R_est,R,V,T] = ...
            viod_adap_obsv(r,v,mu,start_day,obsv_cap,da,noise_vect(i));
        if do_avg_err
            ap_obsv = zeros(1,3);
            for j = 1:obsv_cap
                ap_obsv = ap_obsv + ...
                          TimeProp_V4(R_est(j,:),V(j,:),mu,time_ref-T(j));
            end
            ap_obsv = ap_obsv / obsv_cap;
        else
            ap_obsv = TimeProp_V4(R_est(end,:),V(end,:),mu,time_ref-T(end));
        end
        
    end
    
    errs_obsv(:,i) = errs_obsv(:,i) + (ap_obsv - ap_ref)';
    time_obsv(i) = T(end)-T(1);

    end_day(i) = T(end);
    
end

% ======== Equal Time ========

errs_eqsp = zeros(3,length(noise_vect));
time_eqsp = zeros(length(noise_vect),1);

for i = 1:length(noise_vect)
    
    for ss = 1:samp_size
    
        [R_est,R,V,T] = ...
            viod_eqsp(r,v,mu,start_day,end_day(i),obsv_cap,noise_vect(i));
        if do_avg_err
            ap_eqsp = zeros(1,3);
            for j = 1:obsv_cap
                ap_eqsp = ap_eqsp + ...
                          TimeProp_V4(R_est(j,:),V(j,:),mu,time_ref-T(j));
            end
            ap_eqsp = ap_eqsp / obsv_cap;
        else
            ap_eqsp = TimeProp_V4(R_est(end,:),V(end,:),mu,time_ref-T(end));
        end
        
    end
    
    errs_eqsp(:,i) = errs_eqsp(:,i) + (ap_eqsp - ap_ref)';
    time_eqsp(i) = T(end)-T(1);
    
end

% ======== Error Plot ========

figure(2)
loglog(noise_vect,vecnorm(errs_obsv/obsv_cap,2))
hold on
loglog(noise_vect,vecnorm(errs_eqsp/obsv_cap,2))
hold off
grid(gca,'minor')
grid on
legend('Equal Turn','Equal Time')
xlabel('Noise Magnitude (m/s)')
ylabel('Error at Apoapsis (km)')
 title('Error Convergence Rate')

