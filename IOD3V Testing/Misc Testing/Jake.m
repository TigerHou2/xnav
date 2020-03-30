close all hidden
clear;clc

Parameters;
mu = sun.mu;
res = 100;
dur = 400;
start_day = 20;
interval = 100;


%% Transfer 1
r = [-1.439e8, -4.397e7, 2.588e4];
v = [8.206, -30.898, -1.583];

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Green')
plot3(0,0,0,'bo')
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

view(cross(r,v))
pbaspect([1 1 1])
axis equal
grid on
hold off

%% Earth
r = [-1.439e8, -4.397e7, 2.588e4];
v = [8.248, -28.600, 0.0023];

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Blue')
plot3(0,0,0,'bo')
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

view(cross(r,v))
pbaspect([1 1 1])
axis equal
grid on
hold off

%% Mars
r = [1.931e8, -7.171e7, -6.225e6];
v = [9.302, 24.798, 0.292];

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Red')
plot3(0,0,0,'bo')
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

view(cross(r,v))
pbaspect([1 1 1])
axis equal
grid on
hold off