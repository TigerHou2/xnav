% % Propagator requirements
%
% 1. propagate a spacecraft trajectory for the restricted two-body problem
% given initial conditions: r, v, mu
%
% 2. output r, v for any of the following specified inputs: t (time since
% epoch), f (true anomaly)
%
% 3. verify the output of the solver
% 
%
% % Solver requirements
%
% 1. given mu and three inputs of r, output the orbital parameters of the
% spacecraft trajectory
%

% import physical parameters
Parameters

% pos = [];
% vel = [];
% r = [-4743, 4743, 0] * 1000;
% v = [-5.879, -4.223, 0] * 1000;
% 
% [r1, v1_ref] = TimeProp(r, v, earth.mu, 50);
% [r2, v2_ref] = TimeProp(r, v, earth.mu, 3100);
% dt = 3100-50;

r1 = [1,0,0];
r2 = [0.723*cosd(135), 0.723*sind(135),0];

[a, v1, v2] = Lambert_Solver(r1, r2, 1, 1.978, 0.01);
% [a, v1, v2] = Lambert_Solver(r1, r2, earth.mu, dt, 5e6, 100, 0.001);

v1
v1_ref

% for i = 0:50:5000
%     [r1,v1] = TimeProp(r, v, earth.mu, i);
%     pos = [pos; r1];
%     vel = [vel; v1];
% end
% 
% figure
% plot3(pos(:,1), pos(:,2), pos(:,3))
% figure
% origin = zeros(size(pos(:,1)));
% quiver3(origin, origin, origin ,pos(:,1), pos(:,2), pos(:,3))