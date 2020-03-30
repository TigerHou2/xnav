close all
clear;clc

addpath('../')

Parameters


% ======================

r1 = earth.a; % departure planet
r2 = neptune.a; % arrival planet
mu = sun.mu; % gravitational parameter [km^3/s^2]

theta = pi/4;
phi   = pi/4;
r = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];

r = r * r1 * AU;
v = cross([0,0,1],r);
v = v / norm(v) * sqrt(mu*(2/norm(r)-2/(r1+r2)/AU));

ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit
start_day = 0.01;

step = res/dur;

T = linspace(step,dur,res);
T = T';

% ======================

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
pos = pos(:,1:3);
e = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;
e_norm = norm(e/AU);


vel_ang = zeros(size(T));
fpa = zeros(size(T));
true_anom = zeros(size(T));
f = zeros(size(T));

figure(100)
hold on
for i = 1:size(T,1)
    
    [R,V] = TimeProp_Universal_V2(r,v,mu,T(i));
    [~,~,~,~,~,ff] = Get_Orb_Params(R,V,mu);
    vel_ang(i) = acos(dot(V,e)/norm(V)/norm(e)) - pi/2;
    if dot(R,V) < 0
        vel_ang(i) = 2*pi - vel_ang(i);
    end
    quiver3(0,0,0,V(1),V(2),V(3))
    f(i) = ff;
end
hold off

figure(1)
plot(f,vel_ang)
figure(2)
plot(T,vel_ang)




