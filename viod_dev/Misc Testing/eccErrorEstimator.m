close all
clear;clc

addpath('functions')

eccs = [0.1 0.5 0.9];
a = 1;
mu = 1;
period = 0.1;
f0deg = 160;
f0 = deg2rad(f0deg);

period = period * 2*pi * sqrt(a^3/mu);
errVect = nan(size(eccs));

switch f0deg
    case 0
        ref = [7.398e-3 1.042e-3 1.428e-4];
    case 20
        ref = [7.849e-3 1.287e-3 1.525e-4];
    case 45
        ref = [8.978e-3 1.935e-3 1.844e-4];
    case 90
        ref = [1.280e-2 6.463e-3 4.163e-4];
    case 120
        ref = [1.580e-2 2.024e-2 1.555e-3];
    case 160
        ref = [1.838e-2 8.070e-2 1.118e-1];
    case 180
        ref = [1.843e-2 9.488e-2 6.005e-1];
    otherwise
        error('This f0 value is not catalogued!')
end

refDat = [7.398e-3 1.042e-3 1.428e-4;
          7.849e-3 1.287e-3 1.525e-4;
          8.978e-3 1.935e-3 1.844e-4;
          1.280e-2 6.463e-3 4.163e-4;
          1.580e-2 2.024e-2 1.555e-3;
          1.838e-2 8.070e-2 1.118e-1;
          1.843e-2 9.488e-2 6.005e-1];
ref_f0 = [0 20 45 90 120 160 180]';

plot(ref_f0,refDat(:,1)./refDat(:,2));
plot(ref_f0,refDat(:,2)./refDat(:,3));

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
    error = sqrt((1-e^2)/(1+2*e*cos(f0)+e^2));
    disp(['df = ' num2str(rad2deg(df))])
%     error = error / df^2;

    figure;
    plot(Mvect,fvect)
    
    dfmin = min(fvect(2:end)-fvect(1:end-1));
    dfmax = max(fvect(2:end)-fvect(1:end-1));
    error = error / df^2;
      
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
%     disp(R)

    
    % calculate error
    error = error;
    errVect(i) = error;
    disp(['Error = ' num2str(error)])
    disp(' ')
end

errVect(1:end-1) ./ errVect(2:end)
ref(1:end-1) ./ ref(2:end)