% close all
clear;clc

addpath('..\functions')

mu = 1;
a = 1.56e5;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);
    
orbitParams = [a,e,i,o,w,f];

noise = 1e-6;

numObsv = 3;

numSims = 2000;
durVect = 0.05:0.05:0.9;
durVect = durVect * 2*pi;
selObsv = 1;

errVect = nan(numSims,length(durVect));
dfvect = nan(size(durVect));

for i = 1:length(durVect)
    dur = durVect(i);
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + dur;
    v = nan(numObsv,3);
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);
    
    if fend < f
        fend = fend + 2*pi;
    end
    dfvect(i) = fend-f;
%     fvect = linspace(f,fend,numObsv);
%     fvect = mod(fvect,2*pi);
    
    for s = 1:numSims
        nvect = randn(numObsv,3) * noise;
        for j = 1:numObsv
            orbitParams(6) = fvect(j);
            [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
        end
        orbitParams(6) = fvect(selObsv);
        rRef = Get_Orb_Vects(orbitParams,mu);
        r = hodo(v+nvect,mu);
        r = r(selObsv,:)';
        errVect(s,i) = norm(r-rRef) / norm(rRef) * 100;
    end
end

%%
figure(1)
hold on
% plot(rad2deg(durVect),mean(errVect),'LineWidth',1.5)
plot(rad2deg(dfvect),mean(errVect),'LineWidth',1.5)
xlabel('Mean Anomaly Span, degrees')
ylabel('Error Avg. \%')
set(gca,'FontSize',16)
grid on
hold off
