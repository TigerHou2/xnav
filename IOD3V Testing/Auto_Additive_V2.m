%% Notes
%
% Created:  01/30/2020
% Modified: 02/03/2020
% Author:   Tiger Hou
%
% Changes:
%   Previous Version: Auto_Additive_Obsv.m
%       - Runs Auto_Trial.m instead of Auto_Test_Ext.m, which brings the
%           orbit measurement point calculations outside of the loop.
%           This drastically improves speed for large sample sizes.
%
% Dependencies (in call order):
%       Parameters.m
%       PlanetScan.m
%       TimeProp_Universal_V2.m
%       Auto_Trial.m
%       Data_Consolidate_Additive.m
%       Plot_Data.m
%       Sort_Orbit.m
%       Get_Orb_Params.m
%       Cut_Noodle.m
%
% This script calculates the error magnitude of six orbital elements
%   for a range of velocity observation counts at fixed intervals
%   for a fixed orbit and fixed starting conditions.
% This script visualizes results both in standard deviation plots
%   as well as orbit simulations. 
%
% Sources:
%	- Initial Orbit Determination from Three Velocity Vectors


%% Initial Conditions
%       Set up the Monte Carlo environment and load data.

close all hidden
clear;clc

Parameters

% ======================================
% Dependencies
% ======================================

use_progbar = false;

% ======================================
% Orbital Parameters
% ======================================

r1 = earth.a; % departure planet
r2 = neptune.a; % arrival planet
mu = sun.mu; % gravitational parameter [km^3/s^2]

theta = pi/4;
phi   = pi/4;
r = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];

r = r * r1 * AU;
v = cross([0,0,1],r);
v = v / norm(v) * sqrt(mu*(2/norm(r)-2/(r1+r2)/AU));

% additional remarks
planet = PlanetScan(r1,r2,planets);
if not(isempty(planet))
    notes = [planet.start '-' planet.end ' Transfer Orbit'];
else
    notes = input('Simulation Description', 's');
end


% ======================================
% Orbit Visualization
% ======================================

ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
obsvPlot = 070;                         % observation count to plot
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit


% ======================================
% Testing Conditions
% ======================================

noise     = 10;             % 1 sigma sensor noise [m/s]
start_day = 0.01;            % measurement start day [days]

% ======================================
% Monte Carlo Parameters
% ======================================

% 83563
rng_val  = 01;              % rng seed
obsv_int = 60;             % observe every n-day interval
obsv_min = 010;             % minimum number of observations
obsv_max = 090;             % maximum number of observations
trials   = 050;          	% trials per interval (~0.07s / trial)

% ======================================
% Data Settings
% ======================================

use_f     = true;       % use true anomaly instead of days on x-axis

data_files = dir('Data/*.mat');
fName = [   planet.start(1) planet.end(1) ...
            '_St'   num2str(start_day) ...
            '_Ns'   num2str(noise) ...
            '_Ob'   num2str(obsv_min) '-' ...
                    num2str(obsv_max) ...
            '_Int'  num2str(obsv_int) ...
            '_Samp' num2str(trials) ...
            '_Rng'  num2str(rng_val) ...
            '.mat'];
for i = 1:size(data_files,1)
    if strcmp(fName,data_files(i).name)
        f = warndlg('An identical file already exists.', 'Warning');
        break
    end
end

disp(['Output File Name: ' fName])


%% Monte Carlo Simulation
%       Collect data through Monte Carlo and perform rudimentary
%       statistical analysis without sorting data into elliptical
%       and hyperbolic cases.

profile on
tic

interval_count = length(obsv_min:obsv_max);

% initialize a progress bar window to show time elapsed and % completion
% if use_progbar
%     prog_bar = parfor_wait(50,'Waitbar',true);
% end
tot_iter = interval_count*trials;

% initializa the dataCell
% the dataCell stores sim data for all trials for one particular interval
dataCell(trials,:) = {0, 0, zeros(1,6), 0};
for i = 1:trials
    dataCell(i,:) = {0, 0, zeros(1,6), 0};
end

% initialize the cellofCells
% the cellofCells stores dataCells for all intervals
cellofCells(interval_count) = {dataCell};
for i = 1:interval_count
    cellofCells(i) = {dataCell};
end

% we need to leftshift the interval because there's a bug in MATLAB
% which causes indexing to not work with two operations
% e.g. (i-2) will work, but (i-3+1) will not work.
obsv_min_LS = obsv_min - 1;

% find position and velocity vectors on at observation points
t = zeros(1,obsv_max);  % initialize timestamp vector for observations
for i = 1:obsv_max
    % calculate the timestamps as days since periapse
    t(i) = start_day + (i-1)*obsv_int;
end
R = zeros(obsv_max,3);     % initialize matrix of position vectors
V = zeros(obsv_max,3);     % initialize matrix of velocity vectors
for i = 1:obsv_max
    % calculate position and velocity vectors from orbit propagator
    [R(i,:),V(i,:)] = TimeProp_Universal_V2(r,v,mu,t(i));
end

% simulation
% iterate over all intervals
rng(rng_val);   % seed the rng. no need to seed inside loop because
                %   each loop the observations change size
                %   and also the loop is parallelized.
parfor i = obsv_min:obsv_max
    DT = obsv_int * (i);    % the total interval of all observations
    disp(['Current Interval: ' num2str(i)]) % print to console
    
    C = combnk(1:i,2);           % permutations of two velocity vectors i,j
    comb_n = size(C,1);
    % using sparse matrix to reduce memory consumption
    A = spalloc(comb_n*4,i*2,7*i*(i-1)); % LHS matrix of eq. 36, see line 6
    B = zeros(comb_n*4,1);               % RHS matrix of eq. 36, see line 6
    
    % iterate over all trials
    for j = 1:trials
        % get error data for each trial
        % note that only the first i vectors of R and V are used
        %   because we are only making i observations
        [err,df] = Auto_Trial(R(1:i,:),V(1:i,:),A,B,C,mu,noise,i); % trial
        currCell = cellofCells{i-obsv_min_LS};    % fetch data cell
        currCell(j,:) = {DT, df, err(1:6), boolean(err(7))}; % push data
        cellofCells(i-obsv_min_LS) = {currCell};  % commit data into cell
    end
end

% restructure data from cell format into struct
dataCells = cat(1,cellofCells{:});
colHeadings = {'int', 'df', 'orb', 'is_elliptic'};
data = cell2struct(dataCells, colHeadings, 2);

% sound notification when complete
load handel
sound(y*0.2,Fs)

% Consolidate Data
%       Sort simulation data by elliptical and hyperbolic cases
%       then perform statistical analysis of mean and standard deviation.

init_cond.rng = rng_val;        % rng seed for sim
init_cond.day = start_day;      % measurements start x days after periapse
init_cond.min = obsv_min;       % minimum number of measurement days
init_cond.max = obsv_max;       % maximum number of measurement days
init_cond.obi = obsv_int;       % observe every x days
init_cond.num = trials;         % trials to run for each measurement bin
init_cond.rms = noise;          % noise applied
init_cond.mu  = mu;             % gravitational parameter of orbiting body
init_cond.pos = r;              % initial position vector
init_cond.vel = v;              % initial velocity vector
init_cond.notes = notes;        % notes on simulation setup
init_cond.use_f = use_f;        % use true anomaly as bins instead of days
init_cond.dur   = dur;          % orbital period of true orbit
init_cond.res   = res;          % number of points to plot on IOD orbits
init_cond.interval = obsvPlot;  % observation count used for display
init_cond.fName = fName;        % file name, which contains setup info

[data_ell,data_hyp,mean_ell,std_ell,mean_hyp,std_hyp] = ...
    Data_Consolidate_Additive(data,init_cond);

save(['Data/' fName], 'init_cond', 'data');

toc

disp('Done!')


%% Plot Error Values

close all
Plot_Data(mean_ell, std_ell, 'Elliptic', use_f); % elliptic orbit data

% save fullscreen image of elliptic orbit data
set(gcf, 'Position', get(0, 'Screensize'));
F    = getframe(gcf);
imwrite(F.cdata, ['Plots/' fName(1:end-4) '.png'], 'png')

Plot_Data(mean_hyp, std_hyp, 'Hyperbolic', use_f); % hyperbolic data
Sort_Orbit(data_ell, data_hyp); % stacked histogram of ell vs. hyp orbits


%% Plot True Orbit

pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

figure(4)
hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',2.5,'Color','Red')
plot3(0,0,0,'bo')
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

pos_m = [];
vel_m = [];
for i = start_day:obsvPlot*obsv_int/40:start_day+obsvPlot*obsv_int
    [rm,vm] = TimeProp_Universal_V2(r, v, mu, i);
    pos_m = [pos_m; rm];
    vel_m = [vel_m; vm];
end
plot3(pos_m(:,1), pos_m(:,2), pos_m(:,3),'LineWidth',4.0,'Color','Green')

view(cross(r,v))
pbaspect([1 1 1])
axis equal
grid on
hold off


%% Plot Sample Guess Orbit

tic

reps = 12;

o = obsvPlot;

% find position and velocity vectors on at observation points
t = zeros(1,o);  % initialize timestamp vector for observations
for i = 1:o
    % calculate the timestamps as days since periapse
    t(i) = start_day + (i-1)*obsv_int;
end
Rsamp = zeros(o,3);     % initialize matrix of position vectors
Vsamp = zeros(o,3);     % initialize matrix of velocity vectors
for i = 1:o
    % calculate position and velocity vectors from orbit propagator
    [Rsamp(i,:),Vsamp(i,:)] = TimeProp_Universal_V2(r,v,mu,t(i));
end

C = combnk(1:o,2);   % permutations of two velocity vectors i,j
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
% LHS matrix of eq. 36, see line 6
A = spalloc(comb_n*4,o*2,7*o*(o-1));
B = zeros(comb_n*4,1);      % RHS matrix of eq. 36, see line 6

for k = 1:reps

    dt = obsvPlot * obsv_int;
    [~,~,R0,V0] = Auto_Trial(Rsamp(1:o,:),Vsamp(1:o,:),A,B,C,mu,noise,o);
    r1 = R0(1,:);
    v1 = V0(1,:);

    pos = Get_Orb_Points(r1,v1,mu,res,dur,start_day);
    e_vec = ((dot(v1,v1)-mu/norm(r1))*r1-dot(r1,v1)*v1)/mu*AU;

    hold on

    plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue')
    quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Blue')

    pbaspect([1 1 1])
    axis equal
    grid on
    hold off
    
end

toc
disp('Orbit sampling done!')


%% Cut Space Noodle

hold on
point_cloud = ...
    Cut_Noodle(r,v,mu,interval,noise,start_day,v_count,100);
hold off