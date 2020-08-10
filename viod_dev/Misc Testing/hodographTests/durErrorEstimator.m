close all
clear;clc

addpath('..\functions')

mu = 1;
a = 1;
e = 0.5;
f0 = 0;

numObsv = 3;

durVect = [0.01 0.04 0.8];
durVect = durVect * 2*pi;

errVect = nan(size(durVect));

for i = 1:length(durVect)
    dur = durVect(i);
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
    M0 = E0 - e*sin(E0);
    M = M0 + dur;
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    f = fvect(end);
    df = f-f0;
    df = mod(df,2*pi);
    error = 1 / df^2;
    errVect(i) = error;
end

errVect(1:end-1) ./ errVect(2:end)