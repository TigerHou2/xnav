close all
clear;clc

a = 1;
e = 0.7;
mu = 1;

pulsars = [3 1 2;...
           1 0 1;...
           1 2 0]';

ue  = [1 0 0]';
une = [0 1 0]';
ee  = e * ue;
en  = e * une;

res = 1000;
nump = size(pulsars,2);

rp = a*(1-e^2)/(1+e);
vp = sqrt(mu*(2/rp-1/a));
h = rp * vp;

unr_vect = nan(3,res);
v_vect = nan(3,res);
rr_vect = nan(nump,res);
f_vect = rad2deg(linspace(0,2*pi,res));

for i = 1:res
    unr_vect(:,i) = rotz(f_vect(i)) * une + en;
    v_vect(:,i) = unr_vect(:,i) / h * mu;
    for j = 1:nump
        rr_vect(j,i) = dot(v_vect(:,i),pulsars(:,j));
    end
end

figure;
latexify
hold on
for k = 1:nump
    plot(f_vect,rr_vect(k,:),'DisplayName',['Pulsar ' num2str(k)'])
end
xlabel('True Anomaly')
ylabel('Range-Rate')
legend