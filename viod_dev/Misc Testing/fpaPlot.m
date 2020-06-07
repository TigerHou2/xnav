close all
clear;clc

a = 1;
e = 0.7;
mu = 1;

f_vect = linspace(0,2*pi,1000000);
E_vect = 2 * atan(sqrt((1-e)/(1+e)).*tan(f_vect/2));
M_vect = E_vect - e.*sin(E_vect);

E_vect(E_vect<0) = E_vect(E_vect<0) + 2*pi;
M_vect(M_vect<0) = M_vect(M_vect<0) + 2*pi;

rp = a*(1-e^2)/(1+e*cos(0));
vp = sqrt(mu*(2/rp-1/a));

h = rp * vp;

r_vect = a*(1-e^2)./(1+e.*cos(f_vect));
v_vect = sqrt(mu.*(2./r_vect-1/a));

fpa_vect = acos(h./r_vect./v_vect);
fpa_vect(f_vect>pi) = -fpa_vect(f_vect>pi);

figure;
latexify(45,18)
subplot(1,3,1)
plot(f_vect,fpa_vect)
xlabel('True Anomaly')
axis equal
subplot(1,3,2)
plot(E_vect,fpa_vect)
xlabel('Eccentric Anomaly')
axis equal
subplot(1,3,3)
plot(M_vect,fpa_vect)
xlabel('Mean Anomaly')
axis equal
sgtitle('Flight Path Angle vs. Anomalies')

out = M_vect(fpa_vect==max(fpa_vect));
disp(['Mean anomaly at max FPA: ' num2str(out)])

intl_vect = f_vect-fpa_vect;

figure;
latexify
plot(M_vect,intl_vect)
xlabel('Mean Anomaly')
ylabel('Inertial Angle')
