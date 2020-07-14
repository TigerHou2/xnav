%% Adaptive_Meshing_2.m
%
% Created:  03/24/2020
% Modified: 03/24/2020
% Author:   Tiger Hou
%
% Uses a simpler inertial angle determination method to perform
% adaptive meshing for VIOD.


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
start_day = 10000;            % measurement start day [days]


% observation conditions

obsv_cap = 20; % max number of observation allowed
delta_angle = deg2rad(3); % rad, change in inertial angle before next obsv


% seed random number generator

rng_val = 1;
rng(rng_val);


% plotting

ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit


%% Simulation (Adaptive)

tic

R = zeros(obsv_cap,3); % position (ground truth)
V = zeros(obsv_cap,3); % velocity (measured, w/ error)
T = zeros(obsv_cap,1); % time since periapsis


% iterate

next_t = start_day;

for i = 1:obsv_cap
    
    % propagate time, collect next measurement
    
    T(i) = next_t;
    [R(i,:),V(i,:)] = TimeProp_Universal_V2(r,v,mu,T(i));
    
    % store ground truth data to use for propagating
    
    rt = R(i,:);
    vt = V(i,:);
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    V(i,:) = V(i,:) + n_vec;
    vv = V(i,:);
    
    % propagate to find proper dt
    
    step = 0;
    step_size = 0.5;
    
    angle = 0;
    angle_prev = 0;
    
    while true
        
        step = step + step_size;
        [r_next,v_next] = TimeProp_Universal_V2(rt,vt,mu,step);
        
        angle = acos(dot(vv,v_next)/norm(vv)/norm(v_next));
        
        if angle > delta_angle
            next_t = T(i) + step;
            break
        end
        
        if abs(angle-delta_angle) > abs(angle-angle_prev)
            step_size = step_size * 1.5;
        else
            step_size = step_size * 0.7;
        end
        
        angle_prev = angle;
        
    end
    
end

% perform IOD

C = combnk(1:obsv_cap,2);
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
A = spalloc(comb_n*4,obsv_cap*2,7*obsv_cap*(obsv_cap-1));
B = zeros(comb_n*4,1);

R_est = IOD3V_V3(V,A,B,C,mu,'ordered',true);

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


% estimated orbit

pos = Get_Orb_Points(R_est(1,:),V(1,:),mu,res,dur,start_day);
plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue')
    
% measurement positions

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
ap_adap = TimeProp_Universal_V2(R_est(1,:),  V(1,:),mu,dur/2-T(1));
ap_eqsp = TimeProp_Universal_V2(Req_est   ,Veq(1,:),mu,dur/2-Teq(1));

diff_adap = norm(ap_true - ap_adap)
diff_eqsp = norm(ap_true - ap_eqsp)

hold on
scatter3(ap_true(:,1),ap_true(:,2),ap_true(:,3),'MarkerFaceColor','Red')
scatter3(ap_adap(:,1),ap_adap(:,2),ap_adap(:,3),'MarkerFaceColor','Blue')
scatter3(ap_eqsp(:,1),ap_eqsp(:,2),ap_eqsp(:,3),'MarkerFaceColor','Green')
hold off