close all
clear;clc

addpath('..\functions')

mu = 1;
a = 1.56e5;
e = 0.5;
i = deg2rad(45);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(160);

orbitParams = [a,e,i,o,w,f];

noise = 1e-6;
dur = 0.009;
dur = dur * 2*pi;

disp(['df span: ' num2str(rad2deg(dm2df(f,dur,e)))])

E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M0 = E0 - e*sin(E0);
M = M0 + dur;

numSims = 2000;
obsVect = 3:3:100;
selObsv = 1;

errVect = nan(numSims,length(obsVect));

for i = 1:length(obsVect)
    numObsv = obsVect(i);
    v = nan(numObsv,3);
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);
    df = fend-f;
    df = mod(df,2*pi);
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    for j = 1:numObsv
        orbitParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    ncube = randn(numObsv,3,numSims) * noise;
    for s = 1:numSims
        nvect = ncube(:,:,s);
        r = hodo(v+nvect,mu);
        r = r(selObsv,:)';
        errVect(s,i) = norm(r-rRef) / norm(rRef) * 100;
%         if errVect(s,i) > 2000
%             disp([num2str(i) '.' num2str(s) ' large error warning: ' ...
%                   num2str(errVect(s,i)) '%'])
%         end
    end
end

figure(1)
hold on
plot(obsVect,mean(errVect),'LineWidth',1.5)
xlabel('Number of Observations')
ylabel('Error Avg. \%')
set(gca,'FontSize',16)
grid on
hold off

% figure(2)
% plot(Mvect,fvect)
% hold on
% plot([Mvect(1),Mvect(end)],[fvect(1),fvect(end)])
% hold off
% legend('Actual Line','Straight Reference','Location','Best')