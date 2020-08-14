
clear;clc

addpath('..\functions')

f0 = deg2rad(0);

dfvect = deg2rad(linspace(10,350,20));

numObsv = 3;

mu = 1;
a = 1;
e = 0.9;
i = deg2rad(30);
omg = deg2rad(0);
w = deg2rad(0);
noise = 0.0002;
orbitParams = [a,e,i,omg,w,0];

rngSeed = 1;
selObsv = 1;
numSims = 2000;
v = nan(numObsv,3);
errVect = nan(numSims,length(dfvect));

for i = 1:length(dfvect)
    rng(rngSeed)
    fvect = linspace(f0,f0+dfvect(i),numObsv);
    ncube = randn(numObsv,3,numSims) * noise;
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    for s = 1:numSims
        nvect = ncube(:,:,s);
        for j = 1:numObsv
            orbitParams(6) = fvect(j);
            [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
        end
        r = hodo(v+nvect,mu);
        r = r(selObsv,:)';
        errVect(s,i) = norm(r-rRef) / norm(rRef);
    end
end

%% plotting
figure(1)
hold on
DF = rad2deg(dfvect);
power = 2 - 2 * DF / 360;
plot(DF,mean(errVect).*(DF.^power).*DF,'LineWidth',1.5)
% plot(DF,mean(errVect),'LineWidth',1.5)
xlabel('df')
ylabel('error metric')
ax = gca;
ax.FontSize = 18;
grid on
hold off