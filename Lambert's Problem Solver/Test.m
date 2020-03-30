% add frame translate to somewhere not at the origin

close all
clear;clc

Parameters

r = [1 0.2 0.05] * AU;
v = [8, 30, 7];

t1 = 10;
t2 = 50;
t3 = 150;

[r1, v1] = TimeProp(r, v, sun.mu, t1);
[r2, v2] = TimeProp(r, v, sun.mu, t2);
[r3, v3] = TimeProp(r, v, sun.mu, t3);

u1 = -r1 / norm(r1);
u2 = -r2 / norm(r2);

r1_test2V = IOD2V(v1,v2,u1,u2,sun.mu);
% r1_test3V = IOD3V(v1,v2,v3,sun.mu,'ordered',true);
r1_test3V = IOD3V(v1,v2,v3,sun.mu,'omega',[0,0,1],'prograde',true);

r1_test2V
r1_test3V
r1

disp("=======================")

% tf  = t1-t2; % days (not sure why it's t1-t2 and not the other way around)
% m   = 0; % complete nr. of orbits
% 
% [v1_test, v2_test, extremal_distances] = lambert(r1, r2, tf, m, sun.mu);
% 
% v1_test
% v1
% v2_test
% v2


pos = [];
vel = [];

for i = 0:5:500
    [r1,v1] = TimeProp(r, v, sun.mu, i);
    pos = [pos; r1];
    vel = [vel; v1];
end

e_vec = ((dot(v,v)-sun.mu/norm(r))*r - dot(r,v)*v)/sun.mu;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3))
plot3(0,0,0, 'bo')
r_sample = [r1;r2;r3];
v_sample = [v1;v2;v3];
plot3(r_sample(:,1), r_sample(:,2), r_sample(:,3), 'ro')
quiver3(r_sample(:,1), r_sample(:,2), r_sample(:,3), ...
        v_sample(:,1), v_sample(:,2), v_sample(:,3))

e_vec = e_vec*AU;
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3))
pbaspect([1 1 1])
axis equal
grid on
hold off