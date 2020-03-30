close all hidden
clear;clc

Parameters


% -44 for r(2) is a test case
r  = [0.71,  0, 0.71] * AU; % position [km]
v  = [0   , -38, 0];      % velocity [km/s]
mu = sun.mu;                % gravitational parameter [km^3/s^2]

noise     = 0;              % sensor noise in each axis [m/s]
start_day = 200;            % measurement start day [days]

interval = 400;             % time interval scenario to plot
dur      = 6000;            % orbit propagation duration [days]
res      = 100;             % plotting time step

reps = 3;

%%
% ========================
% Plot True Orbit
% ========================

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Red')
plot3(0,0,0, 'bo')
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

pos_m = [];
vel_m = [];
for i = start_day:interval:start_day+2*interval
    [rm,vm] = TimeProp_Universal_V2(r, v, mu, i);
    pos_m = [pos_m; rm];
    vel_m = [vel_m; vm];
end
plot3(pos_m(:,1), pos_m(:,2), pos_m(:,3),'LineWidth',1.5,'Color','Green')

pbaspect([1 1 1])
axis equal
grid on
hold off

%%
% ========================
% Plot Guess Orbit
% ========================

[~,~,R,V] = Auto_Test_Ext(r,v,mu,interval,noise,start_day,reps);
r1 = R(1,:);
v1 = V(1,:);

pos = Get_Orb_Points(r1,v1,mu,res,dur,start_day);
e_vec = ((dot(v1,v1)-mu/norm(r1))*r1-dot(r1,v1)*v1)/mu;
e_vec = e_vec*AU;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue')

pbaspect([1 1 1])
axis equal
grid on
hold off