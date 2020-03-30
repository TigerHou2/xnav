close all
clear;clc

addpath('..')
Parameters

r = [1.0,  0, 0.3] * AU;
v = [4  , 33, 0];

% change the resolution of the plot
% does not affect testing data
step = 15;
duration = 900;

% change the data collection point and separation
% will affect testing results
t1 = 100;
dt = 90;

t2 = t1 + dt;
t3 = t2 + dt;

[r1_true, v1] = TimeProp(r, v, sun.mu, t1);
[r2_true, v2] = TimeProp(r, v, sun.mu, t2);
[r3_true, v3] = TimeProp(r, v, sun.mu, t3);

% add noise to velocity 'measurements'
% assuming variance to be proportional to speed
%   don't know what the error order of magnitude is
%   so we will use this as a test to determine sensor requirements
noise_factor = 0.002;
v1 = v1 + randn(1,3) .* (noise_factor * v1) + randn(1,3) * noise_factor;
v2 = v2 + randn(1,3) .* (noise_factor * v2) + randn(1,3) * noise_factor;
v3 = v3 + randn(1,3) .* (noise_factor * v3) + randn(1,3) * noise_factor;

% solve the 3V IOD problem
[r1,r2,r3,K] = IOD3V(v1,v2,v3,sun.mu,'omega',[0,0,1],'prograde',true);

%================================
figure('Name','Orbits with Noise')
hold on

% we have been able to show r1,r2,r3 give the same orbits
% so it is no longer necessary to plot all three
a      = Propagate(r1,v1,sun.mu,step,duration,true);
a_true = Propagate( r, v,sun.mu,step,duration,true);

disp(['Semi-major axis error: ' num2str(abs(a-a_true)/a_true*100) '%'])

% r_sample = [r1;r2;r3];
% v_sample = [v1;v2;v3];
% plot3(r_sample(:,1), r_sample(:,2), r_sample(:,3), 'ro')
% quiver3(r_sample(:,1), r_sample(:,2), r_sample(:,3), ...
%         v_sample(:,1), v_sample(:,2), v_sample(:,3))
% 
% pbaspect([1 1 1])
% axis equal
% grid on
% hold off

%================================
% figure('Name','Orbits without Noise')
% hold on
% 
% Propagate(r1_true,v1,sun.mu,step,duration);
% Propagate(      r, v,sun.mu,step,duration);
% 
% pbaspect([1 1 1])
% axis equal
% grid on
% hold off