close all
clear;clc

R = 1;
c = [0,0,0];
angs = deg2rad(1:360)';
x = R * cos(angs) + c(1);
y = R * sin(angs) + c(2);
z = zeros(size(x));
points = [x,y,z];

test_vec = [1,2,2];

proj_length = nan(size(x));

for i = 1:length(x)
    proj_length(i) = dot(points(i,:),test_vec)/norm(test_vec);
end

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
