close all
clear;clc
addpath('..\fcns_misc')

mu = 1;
a = 1.56e5;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);

noise = 1e-6;
numObsv = 3;
numSims = 3000;

e0 = 0.5;
f0 = deg2rad(45);
t0 = 0.1 * 2*pi;

e = 0.4;
f = deg2rad(60);
t = 0.1 * 2*pi;

df0 = dm2df(f0,t0,e0);
df  = dm2df(f,t,e);

df0 = mod(df0,2*pi);
df  = mod(df ,2*pi);

% err0 = 1.42e-3;
err0 = viodSim(mu,a,e0,i,o,w,f0,t0,noise,numObsv,numSims);
err = err0 * df0^2/df^2 * (sqrt((1-e0)/(1+e0))+sqrt((1+e0)/(1-e0))) ...
                        / (sqrt((1-e )/(1+e ))+sqrt((1+e )/(1-e )));
errTrue = viodSim(mu,a,e,i,o,w,f,t,noise,numObsv,numSims);

disp(['Baseline error: ' num2str(100*err0) '%'])
disp(['Predicted error:' num2str(100*err ) '%'])
disp(['True error:     ' num2str(100*errTrue) '%'])