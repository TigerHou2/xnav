%% Adaptive_Meshing_3.m
%
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

rng_val = 8;
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
        [R(i,:),V(i,:)] = TimeProp_V4(r,v,mu,T(i));

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
    
    
    % propagate to find proper dt
    
    step = 0;
    step_size = 1;
    
    diff = 0;
    diff_m1 = 0;
    
    while true
        
        % simulate propagating forward using current orbit estimate
        step = step + step_size;
        [rt,vt] = TimeProp_V4(rr,vv,mu,step);
        
        diff = acos(dot(vv,vt)/norm(vv)/norm(vt));
        
        if diff > delta_angle
            delta_t = step;
            break
        end
        
        if abs(diff-delta_angle) > abs(diff-diff_m1)
            step_size = step_size * 1.5;
        else
            step_size = step_size * 0.8;
        end
        
        diff_m1 = diff;
        
    end
    
    % propagate time, collect next measurement
    
    T(i+1) = T(i) + delta_t;
    [R(i+1,:),V(i+1,:)] = TimeProp_V4(r,v,mu,T(i+1));
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    V(i+1,:) = V(i+1,:) + n_vec;
    
end

disp(['Adaptive Simulation complete! (' num2str(T(obsv_cap)-T(1)) ' days)'])

toc
disp(' ')


%% Simulation (Adaptive v2)

tic

R_ad2 = zeros(obsv_cap,3); % position (ground truth)
V_ad2 = zeros(obsv_cap,3); % velocity (measured, w/ error)
T_ad2 = zeros(obsv_cap,1); % time since periapsis


% find true inertial angle change and divide evenly
[~,v_ini] = TimeProp_V4(r,v,mu,T(1));
[~,v_end] = TimeProp_V4(r,v,mu,T(obsv_cap));
delta_angle_ad2 = acos(dot(v_ini,v_end)/norm(v_ini)/norm(v_end)) / (obsv_cap-1);
% delta_angle_ad2 = delta_angle;

% iterate

next_t = start_day;

for i = 1:obsv_cap
    
    % propagate time, collect next measurement
    
    T_ad2(i) = next_t;
    [R_ad2(i,:),V_ad2(i,:)] = TimeProp_V4(r,v,mu,T_ad2(i));
    
    % store ground truth data to use for propagating
    
    rt = R_ad2(i,:);
    vt = V_ad2(i,:);
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    V_ad2(i,:) = V_ad2(i,:) + n_vec;
    vv = V_ad2(i,:);
    
    % propagate to find proper dt
    
    step = 0;
    step_size = 0.5;
    
    angle = 0;
    angle_prev = 0;
    
    while true
        
        step = step + step_size;
        [r_next,v_next] = TimeProp_V4(rt,vt,mu,step);
        
        angle = acos(dot(vv,v_next)/norm(vv)/norm(v_next));
        
        if angle > delta_angle_ad2
            next_t = T_ad2(i) + step;
            break
        end
        
        if abs(angle-delta_angle_ad2) > abs(angle-angle_prev)
            step_size = step_size * 1.2;
        else
            step_size = step_size * 0.8;
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

R_ad2_est = IOD3V_V3(V_ad2,A,B,C,mu,'ordered',true);

disp(['Adaptive Simulation v2 complete! (' num2str(T_ad2(obsv_cap)-T_ad2(1)) ' days)'])

toc
disp(' ')


%% Simulation (Equispaced)

rng(rng_val);

Req =     zeros(obsv_cap,3);
Veq =     zeros(obsv_cap,3);
Teq =     linspace(start_day,T(obsv_cap),obsv_cap);

for i = 1:length(Teq)
    
    % propagate
    
    [Req(i,:),Veq(i,:)] = TimeProp_V4(r,v,mu,Teq(i));
    
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


%% Display

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
idx = obsv_cap;

% adaptive IOD

pos = Get_Orb_Points(R_est(idx,:),V(idx,:),mu,res,dur,start_day);
plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue')
scatter3(R(1:idx,1),R(1:idx,2),R(1:idx,3),'MarkerEdgeColor','Blue');

% adaptive IOD v2

pos = Get_Orb_Points(R_ad2_est(idx,:),V_ad2(idx,:),mu,res,dur,start_day);
plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Black')
scatter3(R_ad2(1:idx,1),R_ad2(1:idx,2),R_ad2(1:idx,3),'MarkerEdgeColor','Black');

% equispace IOD
pos = Get_Orb_Points(Req_est,Veq(1,:),mu,res,dur,start_day);
plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Green')
scatter3(Req(:,1),Req(:,2),Req(:,3),'MarkerEdgeColor','Green');

hold off

legend('True Orbit','Central Body',...
       'Adaptive IOD Orbit','Adaptive Sample Points',...
       'Adaptive IOD v2 Orbit','Adaptive v2 Sample Points',...
       'Equispace IOD Orbit','Equispace Sample Points')


%% Accuracy Check

time_pos = dur/2;
adap_sel = obsv_cap;

ap_true = TimeProp_V4(r,v,mu,time_pos);
ap_adap = TimeProp_V4(R_est(adap_sel,:),V(adap_sel,:),mu,time_pos-T(adap_sel));
ap_adv2 = TimeProp_V4(R_ad2_est(obsv_cap,:),V_ad2(obsv_cap,:),mu,time_pos-T_ad2(obsv_cap));
ap_eqsp = TimeProp_V4(Req_est,Veq(1,:),mu,time_pos-Teq(1));

ap_stat = zeros(1,3);
samp_size = 5;
for i = obsv_cap-samp_size+1:obsv_cap
    ap_stat = ap_stat + TimeProp_V4(R_est(i,:),V(i,:),mu,time_pos-T(i));
end
ap_stat = ap_stat / samp_size;

diff_adap = norm(ap_true - ap_adap);
diff_adv2 = norm(ap_true - ap_adv2);
diff_eqsp = norm(ap_true - ap_eqsp);
diff_stat = norm(ap_true - ap_stat);

disp(['   Adaptive IOD Error: ' num2str(diff_adap) ' km (' ...
        num2str(diff_adap/AU) ' AU)'])
disp(['Adaptive IOD v2 Error: ' num2str(diff_adv2) ' km (' ...
        num2str(diff_adv2/AU) ' AU)'])
disp(['  Equispace IOD Error: ' num2str(diff_eqsp) ' km (' ...
        num2str(diff_eqsp/AU) ' AU)'])
disp(['Statistical IOD Error: ' num2str(diff_stat) ' km (' ...
        num2str(diff_stat/AU) ' AU)'])
disp(' ')


%% Plot Accuracy Check Locations

hold on
scatter3(ap_true(:,1),ap_true(:,2),ap_true(:,3),'MarkerFaceColor','Red',    'DisplayName','True Position')
scatter3(ap_adap(:,1),ap_adap(:,2),ap_adap(:,3),'MarkerFaceColor','Blue',   'DisplayName','Adaptive Position')
scatter3(ap_adv2(:,1),ap_adv2(:,2),ap_adv2(:,3),'MarkerFaceColor','Black',  'DisplayName','Adaptive v2 Position')
scatter3(ap_eqsp(:,1),ap_eqsp(:,2),ap_eqsp(:,3),'MarkerFaceColor','Green',  'DisplayName','Equispace Position')
scatter3(ap_stat(:,1),ap_stat(:,2),ap_stat(:,3),'MarkerFaceColor','Magenta','DisplayName','Statistical Position')
hold off
legend('Location','Best')


%% Out of Plane Error

h = cross(r,v);

diff_adap_vec = (ap_true - ap_adap);
diff_adv2_vec = (ap_true - ap_adv2);
diff_eqsp_vec = (ap_true - ap_eqsp);
diff_stat_vec = (ap_true - ap_stat);

diff_adap_out = dot(diff_adap_vec,h)/norm(h);
diff_adv2_out = dot(diff_adv2_vec,h)/norm(h);
diff_eqsp_out = dot(diff_eqsp_vec,h)/norm(h);
diff_stat_out = dot(diff_stat_vec,h)/norm(h);

disp(['Out of Plane Adaptive IOD Error: ' newline '    ' ...
        num2str(diff_adap_out) ' km (' num2str(diff_adap_out/AU) ' AU)'])
disp(['Out of Plane Adaptive IOD v2 Error: ' newline '    ' ...
        num2str(diff_adv2_out) ' km (' num2str(diff_adv2_out/AU) ' AU)'])
disp(['Out of Plane Equispace IOD Error: ' newline '    ' ...
        num2str(diff_eqsp_out) ' km (' num2str(diff_eqsp_out/AU) ' AU)'])
disp(['Out of Plane Statistical IOD Error: ' newline '    ' ...
        num2str(diff_stat_out) ' km (' num2str(diff_stat_out/AU) ' AU)'])
disp(' ')