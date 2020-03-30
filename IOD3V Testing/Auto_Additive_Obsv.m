%% Fixes
%   - fix rng for sample orbit plotting and noodle plotting
%   - run Neptune at the same starting day
%
% From research meeting:
%   - fix df, change f
%   - look into more than 3 velocity measurements more thoroughly
%   - how to estimate standard deviation given true anomaly estimate?
%
%   - replaced spherical coordinate noise generator with cartesian
%       so the distribution is uniform on a spherical surface
%   - 

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
obsvPlot = 150;                         % observation count to plot
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit


% ======================================
% Testing Conditions
% ======================================

noise     = 10;             % 1 sigma sensor noise [m/s]
start_day = 020;            % measurement start day [days]

% ======================================
% Monte Carlo Parameters
% ======================================

% 83563
rng_val  = 01;              % rng seed
obsv_int = 02;              % observe every n-day interval
obsv_min = 030;             % minimum number of observations
obsv_max = 120;             % maximum number of observations
trials   = 010;          	% trials per interval (~0.07s / trial)

% ======================================
% Data Settings
% ======================================

save_data = false;      % save data?
use_f     = true;       % use true anomaly instead of days on x-axis


%% Monte Carlo Simulation
%       Collect data through Monte Carlo and perform rudimentary
%       statistical analysis without sorting data into elliptical
%       and hyperbolic cases.
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

% simulation
% iterate over all intervals
parfor i = obsv_min:obsv_max
    DT = obsv_int * (i);  % the total interval
    rng(rng_val);
%     if use_progbar
%         prog_bar.Send;
%     end
    disp(['Current Interval: ' num2str(i)])
    % iterate over all trials
    for j = 1:trials
        [err,df] = Auto_Test_Ext(r,v,mu,DT,noise,start_day,i);
        currCell = cellofCells{i-obsv_min_LS};
        currCell(j,:) = {DT, df, err(1:6), boolean(err(7))};
        cellofCells(i-obsv_min_LS) = {currCell};
    end
end

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

[data_ell,data_hyp,mean_ell,std_ell,...
    mean_hyp,std_hyp] = Data_Consolidate_Additive(data,init_cond);

if save_data
    data_files = dir('Data/Data-*.mat');
    if size(data_files,1) == 0
        label = 1;
    else
        data_names = [data_files.name];
        label = max(str2num(cell2mat(regexp(data_names,'\d+','match')')))+1;
    end
    save(['Data/Data-' num2str(label,'%03.f')],...
          'init_cond', 'data');
end

toc

disp('Done!')


%% Plot Error Values

close all
Plot_Data(mean_ell, std_ell, 'Elliptic', use_f);
Plot_Data(mean_hyp, std_hyp, 'Hyperbolic', use_f);
Sort_Orbit(data_ell, data_hyp);


%% Plot True Orbit


pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

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

for k = 1:reps

    dt = obsvPlot * obsv_int;
    [~,~,R,V] = Auto_Test_Ext(r,v,mu,dt,noise,start_day,obsvPlot);
    r1 = R(1,:);
    v1 = V(1,:);

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