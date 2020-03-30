%% Adaptive_Meshing.m
%
% Created:  03/22/2020
% Modified: 03/25/2020
% Author:   Tiger Hou
%
% Optimizes velocity observation times to reduce influence of sensor noise
% during small true anomaly changes, attempts to generate better orbit
% estimates.


%% Initial Conditions


% clear workspace and load data

close all hidden
clear;clc

Parameters


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

obsv_cap = 20; % max number of observation allowed
delta_t = 5; % days, time between observations
delta_angle = deg2rad(6); % rad, change in inertial angle before next obsv


% seed random number generator

rng_val = 2;
rng(rng_val);


% plotting

ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit


%% Simulation (Adaptive)

tic

% perform first three observations

R =     zeros(obsv_cap,3); % position (ground truth)
R_est = zeros(obsv_cap,3); % position (IOD estimate)
V =     zeros(obsv_cap,3); % velocity (measured, w/ error)
T =     zeros(obsv_cap,1); % time since periapsis


C = combnk(1:3,2);
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
A = spalloc(comb_n*4,6,42);
B = zeros(comb_n*4,1);


% iterate for first three measurements
% assumes true orbit is not hyperbolic, so if IOD gives hyperbolic
% then the program increases the initial delta-t
% and reworks the first three measurements.
%
% we're kind of cheating the system by going back in time, 
% but on an actual mission we would just take a lot of measurements
% and use the start, midpoint, and end.

while true
    
    for i = 1:3

        % propagate forward in orbit by delta-t, get position and velocity

        T(i) = start_day + (i-1) * delta_t;
        [R(i,:),V(i,:)] = TimeProp_Universal_V2(r,v,mu,T(i));

        % add noise to simulate velocity measurements w/ error

        n_vec = randn(1,3);
        n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
        V(i,:) = V(i,:) + n_vec;

    end

    % perform IOD using velocity measurements to get IOD positions
    
    RR = IOD3V_V3(V(1:3,:),A,B,C,mu,'ordered',true);
    
    % get geometry of IOD orbit from IOD position and velocity measurements
    
    rr = RR(1,:);
    vv = V(1,:);
    ah_est = norm(rr)/(2 - norm(rr)*dot(vv,vv)/mu);
    
    % if IOD orbit is hyperbolic, redo the above with a larger delta-t
    
    if ah_est <= 0
        delta_t = delta_t + 1;
    else
        break
    end
    
end


% iterate for observation #4 and above

for i = 3:obsv_cap
    
    C = combnk(1:i,2);
    comb_n = size(C,1);
    % using sparse matrix to reduce memory consumption
    A = spalloc(comb_n*4,i*2,7*i*(i-1));
    B = zeros(comb_n*4,1);
    
    
    % perform IOD using velocity measurements to get IOD positions
    
    RR = IOD3V_V3(V(1:i,:),A,B,C,mu,'ordered',true);
    R_est(i,:) = RR(i,:);
    
    
    % get IOD orbit parameters
    
    rr = R_est(i,:);
    vv = V(i,:);
    [~,ee,~,~,~,ff] = Get_Orb_Params(rr,vv,mu);
    
    
    % calculate new inertial angle
    
    angle = acos(dot(vv,ee)/norm(vv)/norm(ee));
    
    
    % propagate to find proper dt
    
    step = 0;
    step_size = 1;
    
    diff = 0;
    diff_m1 = 0;
    
    while true
        
        step = step + step_size;
        [rt,vt] = TimeProp_Universal_V2(rr,vv,mu,step);
        [~,et,~,~,~,ft] = Get_Orb_Params(rt,vt,mu);
        angle_new = acos(dot(vt,et)/norm(vt)/norm(et));
        
        diff = abs(angle_new-angle);
        diff_m1 = diff;
        
        if diff > delta_angle
            delta_t = step;
            break
        end
        
        if abs(diff-delta_angle) > abs(diff-diff_m1)
            step_size = step_size * 1.5;
        else
            step_size = step_size * 0.7;
        end
        
    end
    
    % propagate time, collect next measurement
    
    T(i+1) = T(i) + delta_t;
    [R(i+1,:),V(i+1,:)] = TimeProp_Universal_V2(r,v,mu,T(i+1));
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    V(i+1,:) = V(i+1,:) + n_vec;
    
end

disp('Simulation complete!')

toc


%% Animate

% true orbit

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

figure(1)
hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1,'Color','Red')
plot3(0,0,0,'bo')

view(cross(r,v))
pbaspect([1 1 1])
axis equal
grid on


% estimated orbits

for i = obsv_cap:obsv_cap
    pos = Get_Orb_Points(R_est(i,:),V(i,:),mu,res,dur,start_day);
    plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue')
    pause(1)
end

scatter3(R(:,1),R(:,2),R(:,3),'MarkerEdgeColor','Blue');

hold off


%% Equispaced Comparison

rng(rng_val);

Req =     zeros(obsv_cap,3);
Veq =     zeros(obsv_cap,3);
Teq =     linspace(start_day,T(end),obsv_cap);

for i = 1:length(Teq)
    
    % propagate
    
    [Req(i,:),Veq(i,:)] = TimeProp_Universal_V2(r,v,mu,Teq(i));
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    Veq(i,:) = Veq(i,:) + n_vec;
    
end

C = combnk(1:i,2);           % permutations of two velocity vectors i,j
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
A = spalloc(comb_n*4,i*2,7*i*(i-1)); % LHS matrix of eq. 36, see line 6
B = zeros(comb_n*4,1);               % RHS matrix of eq. 36, see line 6

RReq = IOD3V_V3(Veq(1:i,:),A,B,C,mu,'ordered',true);
Req_est = RReq(1,:);

hold on
pos = Get_Orb_Points(Req_est,Veq(1,:),mu,res,dur,start_day);
plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Green')
scatter3(Req(:,1),Req(:,2),Req(:,3),'MarkerEdgeColor','Green');
hold off

legend('True Orbit','Central Body',...
       'Adaptive IOD Orbit','Adaptive Sample Points',...
       'Equispace IOD Orbit','Equispace Sample Points')


%% Accuracy Check

ap_true = TimeProp_Universal_V2(r,v,mu,dur/2);
ap_adap = TimeProp_Universal_V2(R_est(obsv_cap,:),V(obsv_cap,:),mu,dur/2-T(obsv_cap));
ap_eqsp = TimeProp_Universal_V2(Req_est,Veq(1,:),mu,dur/2-Teq(1));

diff_adap = norm(ap_true - ap_adap);
diff_eqsp = norm(ap_true - ap_eqsp);

disp([' Adaptive IOD Error: ' num2str(diff_adap) ' km'])
disp(['Equispace IOD Error: ' num2str(diff_eqsp) ' km'])

hold on
scatter3(ap_true(:,1),ap_true(:,2),ap_true(:,3),'MarkerFaceColor','Red')
scatter3(ap_adap(:,1),ap_adap(:,2),ap_adap(:,3),'MarkerFaceColor','Blue')
scatter3(ap_eqsp(:,1),ap_eqsp(:,2),ap_eqsp(:,3),'MarkerFaceColor','Green')
hold off