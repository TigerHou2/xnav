close all
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

numObsv = 10;

numSims = 3000;
durVect = 0.05:0.05:0.9;
durVect = durVect * 2*pi;
selObsv = 1;

errVect = nan(numSims,length(durVect));
dfVect = nan(size(durVect));
vMagVect = nan(1,length(durVect));
svdErrVect = nan(numSims,length(durVect));

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
    dfVect(i) = fend-f;

    for j = 1:numObsv
        orbitParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    
    for s = 1:numSims
        nvect = randn(numObsv,3) * noise;
        orbitParams(6) = fvect(selObsv);
        rRef = Get_Orb_Vects(orbitParams,mu);
        r = hodo(v+nvect,mu);
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
    
    vMagVect(i) = min(vecnorm(v,2,2));
    
end

%%
figure(1)
hold on

xVar = rad2deg(durVect);
yRef = mean(errVect);

yVar = mean(svdErrVect)./dfVect.^2;
% yVar = 1 ./ dfVect .^2;
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

plot(xVar,yVar,'LineWidth',1.5)
plot(xVar,yRef,'LineWidth',1.5)

xlabel('Mean Anomaly Span, degrees')
ylabel('Error Avg. \%')
legend('Prediction','Simulation','Location','Best')
set(gca,'FontSize',18)
grid on
hold off
