% close all
clear;clc

addpath('..\functions')

mu = 1;
a = 1.56e5;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(180);

orbitParams = [a,e,i,o,w,f];

noise = 1e-6;
dur = 0.15;
dur = dur * 2*pi;

disp(['df span: ' num2str(rad2deg(dm2df(f,dur,e)))])

E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M0 = E0 - e*sin(E0);
M = M0 + dur;

numSims = 3000;
obsVect = 3:2:35;
selObsv = 1;

errVect = nan(numSims,length(obsVect));
svdErrVect = nan(numSims,length(obsVect));

for i = 1:length(obsVect)
    numObsv = obsVect(i);
    v = nan(numObsv,3);
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    
    % get ground truth velocity
    for j = 1:numObsv
        orbitParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    ncube = randn(numObsv,3,numSims) * noise;
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(:,:,s);
        r = hodoHyp(v+nvect,mu);
        r = r(selObsv,:)';
        errVect(s,i) = norm(r-rRef) / norm(rRef) * 100;
        
        % test svd
        [~,~,V] = svd(v+nvect,0);
        k = V(:,end);
        if k'*[0;0;1]<0
            k = -k;
        end
        svdErrVect(s,i) = norm(k-[0;0;1]);
    end
    
end

%% plotting

xVar = obsVect;
yRef = mean(errVect);
yVar = mean(svdErrVect);
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
yVar = yVar + offset;

figure(1)
hold on
% plot(xVar,yVar,'--','LineWidth',1.5,'HandleVisibility','off')
% plot(xVar,yRef,'LineWidth',1.5,'DisplayName',['Dur = ' num2str(dur*100,2) '\%'])
plot(xVar,yVar,'LineWidth',1.5)
plot(xVar,yRef,'LineWidth',1.5)
xlabel('Number of Observations')
ylabel('Error Avg. \%')
legend('Prediction','Simulation','Location','Best')
legend
set(gca,'FontSize',18)
grid on
hold off

% figure(2)
% plot(Mvect,fvect)
% hold on
% plot([Mvect(1),Mvect(end)],[fvect(1),fvect(end)])
% hold off
% legend('Actual Line','Straight Reference','Location','Best')