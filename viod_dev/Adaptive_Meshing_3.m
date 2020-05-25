%% Adaptive_Meshing_3.m
%
% Author:   Tiger Hou
%
% Compares two velocity orbit determination methods:
% --- measurements separated by equal time
% --- measurements separated by equal angular change in velocity


%% Initial Conditions

% clear workspace and load data
close all hidden
clear;clc
Parameters

% plot settings
ptSize = 26;

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

noise     = 1;             % 1 sigma sensor noise [m/s]
start_day = 400;         % measurement start day [days]

% observation conditions
obsv_cap = 12; % max number of observations allowed
delta_t = 5; % days, time between observations
delta_angle = deg2rad(10); % rad, change in inertial angle before next obsv

% seed random number generator
rng_val = 1;
rng(rng_val);

% plotting
ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit


%% Simulation (Equal Angular Change)
% takes velocity measurements at a fixed angular separation

tic

Rad = zeros(obsv_cap,3); % position (ground truth)
Vad = zeros(obsv_cap,3); % velocity (measured, w/ error)
Tad = zeros(obsv_cap,1); % time since periapsis

% iterate

next_t = start_day;

for i = 1:obsv_cap
    
    % propagate time, collect next measurement
    
    Tad(i) = next_t;
    [Rad(i,:),Vad(i,:)] = TimeProp_V4(r,v,mu,Tad(i));
    
    % store ground truth data to use for propagating
    
    rt = Rad(i,:);
    vt = Vad(i,:);
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    Vad(i,:) = Vad(i,:) + n_vec;
    vv = Vad(i,:);
    
    % propagate to find proper dt
    
    step = 0;
    step_size = 0.5;
    
    angle = 0;
    angle_prev = 0;
    
    while true
        
        step = step + step_size;
        [r_next,v_next] = TimeProp_V4(rt,vt,mu,step);
        
        angle = acos(dot(vv,v_next)/norm(vv)/norm(v_next));
        
        if angle > delta_angle
            next_t = Tad(i) + step;
            break
        end
        
        if abs(angle-delta_angle) > abs(angle-angle_prev)
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

Rad_est = VIOD(Vad,A,B,C,mu);

disp(['Equal time simulation complete! (' num2str(Tad(obsv_cap)-Tad(1)) ' days)'])

toc
disp(' ')


%% Simulation (Equal Time Change)

rng(rng_val);

Req =     zeros(obsv_cap,3);
Veq =     zeros(obsv_cap,3);
Teq =     linspace(start_day,Tad(obsv_cap),obsv_cap);

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

RReq = VIOD(Veq(1:i,:),A,B,C,mu);
Req_est = RReq(1,:);


%% Display

% true orbit

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

figure(1)
hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'--','LineWidth',1,'Color','Black')
scatter3(0,0,0,ptSize,'k','h','Filled')

view(cross(r,v))
pbaspect([1 1 1])
axis equal
grid on

% estimated orbits
idx = obsv_cap;

% adaptive velocity OD

pos = Get_Orb_Points(Rad_est(idx,:),Vad(idx,:),mu,res,dur,start_day);
plot3(pos(:,1), pos(:,2), pos(:,3),'-','LineWidth',1,'Color','Blue')
scatter3(Rad(:,1),Rad(:,2),Rad(:,3),ptSize,'MarkerEdgeColor','Blue');

% equispace velocity OD
pos = Get_Orb_Points(Req_est,Veq(1,:),mu,res,dur,start_day);
plot3(pos(:,1), pos(:,2), pos(:,3),':','LineWidth',1,'Color','Red')
scatter3(Req(:,1),Req(:,2),Req(:,3),ptSize,'MarkerEdgeColor','Red');

hold off

legend('True Orbit','Central Body',...
       'Equal Angle Orbit','Equal Angle Points',...
       'Equal Time Orbit','Equal Time Points')


%% Accuracy Check

time_pos = dur/2;

ap_true = TimeProp_V4(r,v,mu,time_pos);
ap_ad = TimeProp_V4(Rad_est(obsv_cap,:),Vad(obsv_cap,:),mu,time_pos-Tad(obsv_cap));
ap_eq = TimeProp_V4(Req_est,Veq(1,:),mu,time_pos-Teq(1));

diff_ad = norm(ap_true - ap_ad);
diff_eq = norm(ap_true - ap_eq);

disp(['Equal Angle VIOD Error: ' num2str(diff_ad) ' km (' ...
        num2str(diff_ad/AU) ' AU)'])
disp(['Equal Time VIOD Error:  ' num2str(diff_eq) ' km (' ...
        num2str(diff_eq/AU) ' AU)'])
disp(' ')


%% Plot Accuracy Check Locations

hold on
scatter3(ap_true(:,1),ap_true(:,2),ap_true(:,3),ptSize,'MarkerFaceColor','Black','DisplayName','True Position')
scatter3(ap_ad(:,1),ap_ad(:,2),ap_ad(:,3),ptSize,'MarkerFaceColor','Blue', 'DisplayName','Equal Angle Position')
scatter3(ap_eq(:,1),ap_eq(:,2),ap_eq(:,3),ptSize,'MarkerFaceColor','Red',  'DisplayName','Equal Time Position')
hold off
legend('Location','East')
% latexify(19,19)


%% Out of Plane Error

h = cross(r,v);

diff_ad_vec = (ap_true - ap_ad);
diff_eq_vec = (ap_true - ap_eq);

diff_ad_out = dot(diff_ad_vec,h)/norm(h);
diff_eq_out = dot(diff_eq_vec,h)/norm(h);

disp(['Out of Plane Equal Angle OD Error: ' newline '    ' ...
        num2str(diff_ad_out) ' km (' num2str(diff_ad_out/AU) ' AU)'])
disp(['Out of Plane Equal Time OD Error:  ' newline '    ' ...
        num2str(diff_eq_out) ' km (' num2str(diff_eq_out/AU) ' AU)'])
disp(' ')