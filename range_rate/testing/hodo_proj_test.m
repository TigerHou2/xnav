%% range-rate velocity projection
% Description: 
%   this section projects the hodograph onto a pulsar vector
%       and examines the shape of the projection w.r.t. true anomaly.
%   in theory, the projected results should be a sine wave regardless of
%       the orientation of the projection vector w.r.t. the hodogrpah.
%
% RESULTS:
%   confirmed, the projection vs. true anomaly plot is a sine wave.
%   furthermore, the sine wave's ...
%       amplitude = R*cos(beta), and
%       offset    = e*R*cos(beta)*cos(alpha-gamma), 
%   where
%       alpha = the right ascension of the pulsar vector
%       beta  = the declination of the pulsar vector
%       gamma = the direction of the hodograph offset
%   all w.r.t. the hodograph plane.

close all
clear;clc

R = 3;
c = [1,0,0];
angs = deg2rad(1:360)';
x = R * cos(angs) + c(1);
y = R * sin(angs) + c(2);
z = zeros(size(x));
points = [x,y,z];

test_vec = [1,5,5];

proj_length = nan(size(x));

for i = 1:length(x)
    proj_length(i) = dot(points(i,:),test_vec)/norm(test_vec);
end

alpha = atan2(test_vec(2),test_vec(1));
beta = atan2(test_vec(3),norm(test_vec(1:2)));
gamma = atan2(c(2),c(1));

e = norm(c) / R;

Lp = R*cos(beta)*(e*cos(alpha-gamma)+1);
Lm = R*cos(beta)*(e*cos(alpha-gamma)-1);
L_est = Lp - Lm;
L = max(proj_length)-min(proj_length);

plot(rad2deg(angs),proj_length)
hold on
plot(rad2deg(angs),vecnorm(points,2,2));
hold off
legend('Projected','Original')
grid on
grid(gca,'minor')

%%
close all
clear;clc

e_vect = [0,0.1,0.3,0.5,0.7,0.9,0.95,0.99];
f = deg2rad(181:540);
    
figure(1)
hold on

for i = 1:length(e_vect)
    e = e_vect(i);
    E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M = E - e*sin(E);
    plot(f,M);
end

hold off
xlabel('f')
ylabel('M')
