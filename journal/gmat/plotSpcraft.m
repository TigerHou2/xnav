%% plotSpcraft.m
function plotSpcraft(numSims,rngSeed,perturbed)
%
% Author:
%   Tiger Hou
%
% Description:
%   Generates plots for velocity initial orbit determination (VIOD) error
%   with respect to variations in spacecraft parameters, including duration
%   of measurement, number of measurements, and measurement scheduling.
%   This version accounts for perturbations from other solar system objects
%   as well as solar radiation pressure and Earth oblateness. Propagation
%   is performed in GMAT.
%
% Notes:
%   In the case of varying spacecraft parameters, we choose to hold orbital
%   parameters constant. This means all simulations will be conducted with
%   the same semi-major axis, eccentricity, and starting true anomaly.
%   Specifically, we will use the following:
%       - SMA = 1.56e5 DU
%       - ECC = 0.5
%       - TA  = 0

%% initialization

% load spacecraft parameter config case studies
[noiseVect,durVect,obsVect,SMA] = load_spcraft_cases(1);

%% collect ground truth data from .mat file and generate perturbations

% load cell matrix of position, velocity, and mu data
if perturbed == 0
    load('temp\datSpcraft.mat','rArray','vArray','muArray');
else
    load('temp\datSpcraftPerturbed.mat','rArray','vArray','muArray');
end

rng(rngSeed);

% data logging
if perturbed == 0
    filepath = 'data\spcraftParams.mat';
else
    filepath = 'data\spcraftParamsPerturbed.mat';
end
meanOrig = nan(length(obsVect),length(noiseVect),length(durVect));
meanHodo = nan(length(obsVect),length(noiseVect),length(durVect));
 stdOrig = nan(length(obsVect),length(noiseVect),length(durVect));
 stdHodo = nan(length(obsVect),length(noiseVect),length(durVect));

%% apply perturbations and perform VIOD

% line formats
origWidth = 1.5;
hodoWidth = 0.3;
origFormat = '-x';
hodoFormat = '--k';

% color formats
cmap = hsv(length(durVect));

% initialize cell array for custom labeling
labelArray = cell(length(durVect),length(obsVect),2);

for i = 1:length(durVect)
    dur = durVect(i);
for j = 1:length(obsVect)
    numObsv = obsVect(j);
    
% initialize position error matrix
errOrig = nan(numSims,length(noiseVect));
errHodo = nan(numSims,length(noiseVect));

% initialize cell array of noise values
nArray = cell(1,numSims);
% populate noise cell array
rng(rngSeed);
for n = 1:numSims
    noise_unit = randn(numObsv,3);
    noise_unit = noise_unit ./ vecnorm(noise_unit,2,2);
    nArray{n} = noise_unit;
end %numSims

% perform Monte Carlo sim
for n = 1:numSims
    noise_vec = nArray{n};
for k = 1:length(noiseVect)
    rRef = rArray{i,j,1};
    v = vArray{i,j,1};
    mu = muArray{i,j,1};
    noise = noise_vec * noiseVect(k);
    % do orbit determination
    % --- note the scaling for the original method: this is because
    % --- precision issues arise when using canonical units. 
    rOrig = viod((v+noise)*1e4,mu*1e12)/1e4;
    rHodo = hodo(v+noise,mu);
    % we choose to compare the position estimate at the first measurement
    rOrig = rOrig(1,:);
    rHodo = rHodo(1,:);
    % store error data
    errOrig(n,k) = norm(rOrig-rRef);
    errHodo(n,k) = norm(rHodo-rRef);
end %noiseVect
end %numSims

% plot results for each case
color = cmap(i,:);

% --- mean error
figure(1)
latexify('plotSize',[30 18])
hold on
plot(noiseVect,mean(errOrig)/SMA,origFormat,'Color',color,'LineWidth',origWidth)
plot(noiseVect,mean(errHodo)/SMA,hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Energy Method','Hodograph Method','Location','SouthEast')
xlabel('Noise, DU/TU')
ylabel('Position Error Avg., fraction of SMA')
% store annotation data point
labelArray{i,j,1} = [noiseVect(end),mean(errOrig(:,end))/SMA];
latexify('fontSize',18)

% --- standard deviation
figure(2)
latexify('plotSize',[30 18])
hold on
plot(noiseVect,std(errOrig)/SMA,origFormat,'Color',color,'LineWidth',origWidth)
plot(noiseVect,std(errHodo)/SMA,hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Energy Method','Hodograph Method','Location','SouthEast')
xlabel('Noise, DU/TU')
ylabel('Position Error StDev, fraction of SMA')
% store annotation data point
labelArray{i,j,2} = [noiseVect(end),std(errOrig(:,end))/SMA];
latexify('fontSize',18)

% --- log data
meanOrig(j,:,i) = mean(errOrig)'/SMA;
meanHodo(j,:,i) = mean(errHodo)'/SMA;
 stdOrig(j,:,i) =  std(errOrig)'/SMA;
 stdHodo(j,:,i) =  std(errHodo)'/SMA;

end %obsVect
end %durVect

% save data to file
save(filepath,'meanOrig','meanHodo','stdOrig','stdHodo',...
              'durVect','obsVect','noiseVect');

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
for i = 1:length(durVect)
    dur = durVect(i);
for j = 1:length(obsVect)
    numObsv = obsVect(j);
    p = labelArray{i,j,1};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['\# = ' num2str(numObsv,4)],...
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
for i = 1:length(durVect)
    dur = durVect(i);
for j = 1:length(obsVect)
    numObsv = obsVect(j);
    p = labelArray{i,j,2};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['\# = ' num2str(numObsv,4)],...
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
for i = 1:length(durVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',['Dur = ' num2str(durVect(i))]);
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(durVect)
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
for i = 1:length(durVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',['Dur = ' num2str(durVect(i))]);
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(durVect)
    set(h(i),'XData',[],'YData',[]);
end
set(ax2,'Color','none',...
    'XTick',[],'YTick',[],...
    'Box','off') % make legend axes transparent
ax2.XLabel.String = 5h'';
ax2.YLabel.String = '';
end %plotSpcraft.m