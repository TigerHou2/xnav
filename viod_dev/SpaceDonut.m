%% Initial Conditions
%       Set up the Monte Carlo environment and load data.

close all hidden
clear;clc

Parameters

% ======================================
% Orbital Parameters
% ======================================

r1 = earth.a; % departure planet
r2 = mars.a; % arrival planet
mu = sun.mu; % gravitational parameter [km^3/s^2]

theta = pi/4;
phi   = pi/4;
r = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];

r = r * r1 * AU;
v = cross([0,0,1],r);
v = v / norm(v) * sqrt(mu*(2/norm(r)-2/(r1+r2)/AU));

% additional remarks
planet = PlanetScan(r1,r2,planets);
if not(isempty(planet))
    notes = [planet.start '-' planet.end ' Transfer Orbit'];
else
    notes = input('Simulation Description', 's');
end


% ======================================
% Testing Conditions
% ======================================

samples   = 012;            % number of points on the Fibonaccia sphere
noise     = 10;             % 1 sigma sensor noise [m/s]
start_day = 020;            % measurement start day [days]
interval  = 050;            % measurement time interval [days]
v_count   = 03;             % number of velocity measurements per trial
trials    = 010;            % trials per interval (~0.2s / trial)

make_plot = false;
ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
dur       = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res       = 100;                         % # points to plot for orbit

% ======================================
% Data Settings
% ======================================

save_data = false;      % save data?
use_f     = true;       % use true anomaly instead of days on x-axis


%%
tic
params_sample = [];
points = [];

DT = interval;
reps = v_count;

t = zeros(1,reps);
dt = DT / (reps-1);
for i = 1:reps
    t(i) = start_day + (i-1)*dt;
end

R = zeros(reps,3);
V = zeros(reps,3);

for i = 1:reps
    [R(i,:),V(i,:)] = TimeProp_Universal_V2(r,v,mu,t(i));
end

[a,e,i,omg,w,f] = Get_Orb_Params(R(1,:),V(1,:),mu);
params_true = [a,norm(e),i,omg,w,f];

% create Fibonacci sphere of noise
lat = zeros(1,samples);
lon = zeros(1,samples);
for i = 1:samples
    phi = (sqrt(5)+1)/2-1;
    ga = phi * 2*pi;
    lon(i) = ga*i;
    lat(i) = asin(-1+2*i/samples);
end
[x,y,z] = sph2cart(lon,lat,noise*3/1000);
N = [x',y',z'];

hold on
for m = 1:samples
    for n = 1:samples
        for p = 1:samples
            VP = V + [N(m,:);N(n,:);N(p,:)];
            RP = IOD3V_Ext(VP,mu,'ordered',true);
            r1 = RP(1,:);
            v1 = VP(1,:);
            pos = Get_Orb_Points(r1,v1,mu,res,dur,start_day);
            points = [points; pos(:,1:3)];
            if make_plot
%                 plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue');
                scatter3(pos(:,1), pos(:,2), pos(:,3),0.5,'b')
            end
            [a,e,i,omg,w,f] = Get_Orb_Params(r1,v1,mu);
            params_sample = [params_sample; [a,norm(e),i,omg,w,f]];
        end
    end
end
view(cross(r,v))
pbaspect([1 1 1])
hold off
axis equal
toc
%%

hold on
if not(make_plot)
    for j = 1:30
        m = randi(samples);
        n = randi(samples);
        p = randi(samples);
        VP = V + [N(m,:);N(n,:);N(p,:)];
        RP = IOD3V_Ext(VP,mu,'ordered',true);
        r1 = RP(1,:);
        v1 = VP(1,:);
        pos = Get_Orb_Points(r1,v1,mu,res,dur,start_day);
        plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue');
    end
%     shp = alphaShape(points);
%     a = criticalAlpha(shp,'one-region');
%     shp = alphaShape(points,a);
%     plot(shp)
end
view(cross(r,v))
pbaspect([1 1 1])
axis equal
hold off


params_diff = params_sample - params_true;
figure

subplot(2,3,1)
boxplot(params_diff(:,1))
title('Semi-Major Axis')
grid on

subplot(2,3,2)
boxplot(params_diff(:,2))
title('Eccentricity')
grid on

subplot(2,3,3)
boxplot(params_diff(:,3))
title('Inclination')
grid on

subplot(2,3,4)
boxplot(params_diff(:,4))
title('Longitude of Ascending Node')
grid on

subplot(2,3,5)
boxplot(params_diff(:,5))
title('Argument of Periapse')
grid on

subplot(2,3,6)
boxplot(params_diff(:,6))
title('True Anomaly')
grid on