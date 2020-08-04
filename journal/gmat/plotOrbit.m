%% plotOrbit.m
function plotOrbit(numSims,rngSeed,perturbed)
%
% Author:
%   Tiger Hou
%
% Description:
%   Generates plots for velocity initial orbit determination (VIOD) error
%   with respect to variations in orbital parameters, including semi-major
%   axis, eccentricity, and true anomaly. 
%   This version accounts for perturbations from other solar system objects
%   as well as solar radiation pressure and Earth oblateness. Propagation
%   is performed in GMAT.
%
% Notes:
%   In the case of varying orbital parameters, we choose to hold spacecraft
%   parameters constant. This means all simulations will be conducted with
%   equal noise, number of measurements, and duration of observation.
%   Specifically, we will use the following:
%       - Noise = 1e-6 DU/TU
%       - Number of measurements = 3
%       - Duration of observation = 0.1 orbital period

%% initialization

% load orbital parameter config case studies
[smaVect,eccVect,taVect,~,~,~,~,~,~,noise,numObsv] = load_orbit_cases(1);
[~,~,~,~,names,~,~,~,~,~] = load_orbit_cases(0);

%% collect ground truth data from .mat file and generate perturbations

% load cell matrix of position, velocity, and mu data
if perturbed == 0
    load('temp\datOrbit.mat','rArray','vArray','muArray');
else
    load('temp\datOrbitPerturbed.mat','rArray','vArray','muArray');
end

% initialize cell array of noise values
nArray = cell(1,numSims);
% populate noise cell array
rng(rngSeed);
for n = 1:numSims
    noiseVect = randn(numObsv,3);
    noiseVect = noiseVect ./ vecnorm(noiseVect,2,2) * noise;
    nArray{n} = noiseVect;
end %numSims

% data logging
if perturbed == 0
    filepath = 'data\orbitParams.mat';
else
    filepath = 'data\orbitParamsPerturbed.mat';
end
meanOrig = nan(length(eccVect),length(taVect),length(smaVect));
meanHodo = nan(length(eccVect),length(taVect),length(smaVect));
 stdOrig = nan(length(eccVect),length(taVect),length(smaVect));
 stdHodo = nan(length(eccVect),length(taVect),length(smaVect));

%% apply perturbations and perform VIOD

% line formats
origWidth = 1.5;
hodoWidth = 0.3;
origFormat = '-x';
hodoFormat = '--k';

% color formats
cmap = hsv(length(smaVect));

% initialize cell array for custom labeling
labelArray = cell(length(smaVect),length(eccVect),2);

for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
    
% initialize position error matrix
errOrig = nan(numSims,length(taVect));
errHodo = nan(numSims,length(taVect));

% perform Monte Carlo sim
for n = 1:numSims
for k = 1:length(taVect)
    rRef = rArray{i,j,k};
    v = vArray{i,j,k};
    mu = muArray{i,j,k};
    noiseVect = nArray{n};
    % do orbit determination
    % --- note the scaling for the original method: this is because
    % --- precision issues arise when using canonical units. 
    rOrig = viod((v+noiseVect)*1e4,mu*1e12)/1e4;
    rHodo = hodo(v+noiseVect,mu);
    % we choose to compare the position estimate at the first measurement
    rOrig = rOrig(1,:);
    rHodo = rHodo(1,:);
    % store error data
    errOrig(n,k) = norm(rOrig-rRef) / norm(rRef) * 100;
    errHodo(n,k) = norm(rHodo-rRef) / norm(rRef) * 100;
end %taVect
end %numSims

% plot results for each case
color = cmap(i,:);
taDeg = rad2deg(taVect);

% --- mean error
figure(1)
latexify('plotSize',[30 18])
hold on
plot(taDeg,mean(errOrig),origFormat,'Color',color,'LineWidth',origWidth)
plot(taDeg,mean(errHodo),hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Energy Method','Hodograph Method','Location','SouthEast')
xlabel('True Anomaly, deg')
ylabel('Position Error Avg. \%')
% store annotation data point
labelArray{i,j,1} = [taDeg(end),mean(errOrig(:,end))];
latexify('fontSize',18)

% --- standard deviation
figure(2)
latexify('plotSize',[30 18])
hold on
plot(taDeg,std(errOrig),origFormat,'Color',color,'LineWidth',origWidth)
plot(taDeg,std(errHodo),hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Energy Method','Hodograph Method','Location','SouthEast')
xlabel('True Anomaly, deg')
ylabel('Position Error StDev \%')
% store annotation data point
labelArray{i,j,2} = [taDeg(end),std(errOrig(:,end))];
latexify('fontSize',18)

% --- log data
meanOrig(j,:,i) = mean(errOrig)' / 100;
meanHodo(j,:,i) = mean(errHodo)' / 100;
 stdOrig(j,:,i) =  std(errOrig)' / 100;
 stdHodo(j,:,i) =  std(errHodo)' / 100;

end %eccVect
end %smaVect

% save data to file
save(filepath,'meanOrig','meanHodo','stdOrig','stdHodo',...
              'smaVect','eccVect','taVect');

% print data in latex table formatting
regexprep(...
regexprep(...
regexprep(latex(vpa(sym(meanOrig(:,:,1)),10)),...
    '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
    '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
    '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')
regexprep(...
regexprep(...
regexprep(latex(vpa(sym(meanOrig(:,:,2)),10)),...
    '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
    '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
    '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')
regexprep(...
regexprep(...
regexprep(latex(vpa(sym(meanOrig(:,:,3)),10)),...
    '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
    '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
    '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')

%% add annotations

figure(1)
set(gca, 'YScale', 'log')
setgrid
expand(0.07,0.12,0.07,0.05)
axis tight
axlim = axis(gca) .* [1 1 0.3 3];
xlim(axlim(1:2))
ylim(axlim(3:4))
for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
    p = labelArray{i,j,1};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['e = ' num2str(ECC,4)],...
               'EdgeColor','none',...
               'Interpreter','latex',...
               'FontSize',14)
end
end

figure(2)
set(gca, 'YScale', 'log')
setgrid
expand(0.07,0.12,0.07,0.05)
axis tight
axlim = axis(gca) .* [1 1 0.3 3];
xlim(axlim(1:2))
ylim(axlim(3:4))
for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
    p = labelArray{i,j,2};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['e = ' num2str(ECC,4)],...
               'EdgeColor','none',...
               'Interpreter','latex',...
               'FontSize',14)
end
end

%% add orbit SMA legend

figure(1)
ax = gca;
ax2 = copyobj(ax,gcf);
delete(get(ax2,'Children')) % delete its children
hold(ax2,'on')
for i = 1:length(smaVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',names{i});
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(smaVect)
    set(h(i),'XData',[],'YData',[]);
end
set(ax2,'Color','none',...
    'XTick',[],'YTick',[],...
    'Box','off') % make legend axes transparent
ax2.XLabel.String = '';
ax2.YLabel.String = '';

figure(2)
ax = gca;
ax2 = copyobj(ax,gcf);
delete(get(ax2,'Children')) % delete its children
hold(ax2,'on')
for i = 1:length(smaVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',names{i});
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(smaVect)
    set(h(i),'XData',[],'YData',[]);
end
set(ax2,'Color','none',...
    'XTick',[],'YTick',[],...
    'Box','off') % make legend axes transparent
ax2.XLabel.String = '';
ax2.YLabel.String = '';
end %plotOrbit.m