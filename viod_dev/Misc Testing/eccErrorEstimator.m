close all
clear;clc

addpath('functions')

eccs = [0.1 0.5 0.9];
a = 1;
mu = 1;
period = 0.1;
f0deg = 0;
f0 = deg2rad(f0deg);

period = period * 2*pi * sqrt(a^3/mu);
errVect = nan(size(eccs));

switch f0deg
    case 0
        ref = [8.417e-3 1.136e-3 1.420e-4];
    case 20
        ref = [8.635e-3 1.350e-3 1.517e-4];
    case 45
        ref = [9.280e-3 1.937e-3 1.835e-4];
    case 90
        ref = [1.230e-2 6.227e-3 4.149e-4];
    case 120
        ref = [1.573e-2 2.011e-2 1.567e-3];
    case 160
        ref = [2.038e-2 8.754e-2 1.224e-1];
    case 180
        ref = [2.079e-2 1.053e-1 6.672e-1];
    otherwise
        error('This f0 value is not catalogued!')
end

for i = 1:length(eccs)
    e = eccs(i);
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    E = kepler(M,e);
    Mvect = linspace(M0,M,6);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    f = fvect(end);
    df = f-f0;
    df = mod(df,2*pi);
    error = sqrt((1-e^2)/(1+2*e*cos(f0)+e^2)) ...
          * 1/df^2;
      
    % measure how even our f distribution
%     ro = 0;
%     for j = 1:length(fvect)
%         temp = fvect;
%         temp(j) = [];
%         ro = ro + min(abs(temp-fvect(j))) / length(fvect);
%     end
%     re = 0.5*sqrt((f-f0)/length(fvect));
%     R = ro/re;

    R = 1 - 0.5 * sum(abs((fvect-mean(fvect)))/mean(fvect)/length(fvect));

    disp(['R = ' num2str(R)])
    
    % calculate error
    error = error;
    errVect(i) = error;
    disp(['Error = ' num2str(error)])
end

errVect(1:end-1) ./ errVect(2:end)
ref(1:end-1) ./ ref(2:end)