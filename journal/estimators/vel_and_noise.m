%% vel_and_noise.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Predicts circle fitting error for different velocity-to-noise ratios.
%

%% initialization
close all
clear;clc
addpath('..\fcns_misc')
savePath = 'plots\';
latexify

%% setup
R = 1;
% e = 0.5;
c = [0;1];
df = 0.1*2*pi;
% f0 = deg2rad(90);
numObsv = 10;
numSims = 3000;

nvrVect = linspace(1e-4,1e-2,100);


eVect = [0.3 0.8];
fVect = deg2rad([0 90 180]);

for ee = 1:length(eVect)
for ff = 1:length(fVect)
    
if ee == 2 && ff == 3
    continue
end
    
e = eVect(ee);
f0 = fVect(ff);
    
c = c / norm(c) * e * R;
f0 = f0 + atan2(c(2),c(1));
f1 = f0 + df;
fvect = linspace(f0,f1,numObsv);
pos = [cos(fvect);sin(fvect)] * R + c;
posRef = [cos(f0);sin(f0)] * R + c;

ncube = randn(2,numObsv,numSims);
ncube = ncube ./ vecnorm(ncube,2,1);

errDat = nan(numSims,length(nvrVect));

for i = 1:length(nvrVect)
    nGauss = normrnd(0,nvrVect(i),1,numObsv,numSims);
    nGauss = repmat(nGauss,2,1,1);
    ncube_loc = ncube .* nGauss;
    
    for s = 1:numSims
        noise = ncube_loc(:,:,s);
        [x,y,g_R] = hyperfit(pos+noise);

        g_c = [x;y];
        g_e = norm(g_c)/g_R;
        g_c = g_c / norm(g_c) * g_e * g_R;

        g_pos = [cos(f0);sin(f0)] * g_R + g_c;
        errDat(s,i) = norm(posRef-g_pos) / norm(posRef) * 100;
    end
end

xVar = nvrVect;
yVar = sqrt(mean(errDat.^2));
figure(1)
hold on
plot(xVar,yVar,'LineWidth',1.25,...
       'DisplayName',['e=' num2str(e) ', f0=' num2str(rad2deg(f0)) '$^o$'])
hold off

end
end

xlabel('Noise-to-Velocity Ratio')
ylabel('Mean Squared Velocity Error \%')
legend('Location','Best')
latexify(13,10,14)
setgrid
expand
svnm = [savePath 'NVR'];
print(svnm,'-dpdf','-bestfit')