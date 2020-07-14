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

use_progbar = true;

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
interval = 030;                         % time interval scenario to plot
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit


% ======================================
% Testing Conditions
% ======================================

noise     = 10;             % 3 sigma sensor noise [m/s]
start_day = 020;            % measurement start day [days]

% ======================================
% Monte Carlo Parameters
% ======================================

% 83563
rng_val = 88127;                % rng seed
int_min = 030;                  % minimum interval [days]
int_max = 200;                  % maximum interval [days]
int_res = (int_max-int_min)/50; % do trials every n-day interval
v_count = 03;                   % number of velocity measurements per trial
trials  = 010;                  % trials per interval (~0.07s / trial)

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

% initialize a progress bar window to show time elapsed and % completion
if use_progbar
    prog_bar = parfor_wait(50,'Waitbar',true);
end
tot_iter = (length(int_min:int_res:int_max))*trials;
dataCell(trials,:) = {0, 0, zeros(1,6), 0};
for i = 1:trials
    dataCell(i,:) = {0, 0, zeros(1,6), 0};
end
cellofCells(51) = {dataCell};
for i = 1:51
    cellofCells(i) = {dataCell};
end

parfor i = 1:51
    itv = int_min + (i-1)*int_res;
    rng(rng_val);
    if use_progbar
        prog_bar.Send;
    end
    for j = 1:trials
        [err,df] = Auto_Test_Ext(r,v,mu,itv,noise,start_day,v_count);
        currCell = cellofCells{i};
        currCell(j,:) = {itv, df, err(1:6), boolean(err(7))};
        cellofCells(i) = {currCell};
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
init_cond.min = int_min;        % minimum number of measurement days
init_cond.max = int_max;        % maximum number of measurement days
init_cond.irs = int_res;        % increment in measurement days per bin
init_cond.vct = v_count;        % number of observations per measurement, at least 3
init_cond.num = trials;         % trials to run for each measurement bin
init_cond.rms = noise;          % noise applied
init_cond.mu  = mu;             % gravitational parameter of orbiting body
init_cond.pos = r;              % initial position vector
init_cond.vel = v;              % initial velocity vector
init_cond.notes = notes;        % notes on simulation setup
init_cond.use_f = use_f;        % use true anomaly as bins instead of days
init_cond.dur   = dur;          % orbital period of true orbit
init_cond.res   = res;          % number of points to plot on IOD orbits
init_cond.interval = interval;  % measurement day used for display

[data_ell,data_hyp,mean_ell,std_ell,...
    mean_hyp,std_hyp] = Data_Consolidate_V2(data,init_cond);

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
for i = start_day:(interval/(v_count-1)/5):start_day+interval
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

reps = 12;

for k = 1:reps

    [~,~,R,V] = Auto_Test_Ext(r,v,mu,interval,noise,start_day,v_count);
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


%% Cut Space Noodle

hold on
point_cloud = ...
    Cut_Noodle(r,v,mu,interval,noise,start_day,v_count,100);
hold off