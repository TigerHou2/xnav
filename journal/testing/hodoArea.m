close all
clear;clc

addpath('../fcns_orb')

e = 0.2;
R = 1;

dM = deg2rad(135);
f0 = deg2rad([-45, 135]);
E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
M0 = E0 - e*sin(E0);

Mf = M0 + dM;

Ef = kepler(Mf,e);
ff = 2 * atan(sqrt((1+e)/(1-e))*tan(Ef/2));
ff = mod(ff,2*pi);

f0 = f0 + pi/2;
ff = ff + pi/2;

angs = linspace(0,2*pi,1000);
ang1 = linspace(f0(1),ff(1),100);
ang2 = linspace(f0(2),ff(2),100);
c = [0;e*R];

scatter([0,c(1)],[0,c(2)],'k')
hold on
plot(cos(angs)+c(1),sin(angs)+c(2),'g','linewidth',0.7)
plot(cos(ang1)+c(1),sin(ang1)+c(2),'r','linewidth',1.8)
plot(cos(ang2)+c(1),sin(ang2)+c(2),'b','linewidth',1.8)
xline(0);
yline(0);
plot([0,cos(ang1(end))+c(1)],[0,sin(ang1(end))+c(2)],'r')
plot([0,cos(ang1(1))+c(1)],[0,sin(ang1(1))+c(2)],'r')
plot([0,cos(ang2(end))+c(1)],[0,sin(ang2(end))+c(2)],'b')
plot([0,cos(ang2(1))+c(1)],[0,sin(ang2(1))+c(2)],'b')
hold off
axis equal