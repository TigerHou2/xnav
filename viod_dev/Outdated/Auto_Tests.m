%% Fixes
%
%   - plot each of the 100 cases across all time intervals to get 100
%   strands. the goal is to see trends for each rng sample.
%   - check zero noise gives true orbit
%   - fix orbit direction error in IOD3V_V2
%
%   - find cause for spikes
%   - find max eccentricity / farthest planet where no hyperbolic cases
%   start to happen


%% Initial Conditions
%       Set up the Monte Carlo environment and load data.

close all hidden
clear;clc

Parameters

% ======================================
% Orbital Parameters
% ======================================

theta = pi/4;
phi   = pi/4;
r = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];

r1 = earth.a; % departure planet
r2 = jupiter.a; % arrival planet

mu = sun.mu; % gravitational parameter [km^3/s^2]

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
% Testing Conditions
% ======================================

noise     = 10;              % sensor noise in each axis [m/s]
start_day = 000;            % measurement start day [days]

% ======================================
% Monte Carlo Parameters
% ======================================

rng_val = 38322;            % rng seed
int_min = 010;              % minimum interval [days]
int_max = 100;              % maximum interval [days]
int_res = 2;                % do trials every nth interval
v_count = 3;                % number of velocity measurements per trial
trials  = 008;              % trials per interval (approx. 0.2 seconds per trial)

% ======================================
% Data Settings
% ======================================

save_data = false;           % save data?
use_f     = true;           % use true anomaly instead of days on x-axis

% ======================================
% Orbit Visualization
% ======================================

ah = norm(r)/(2 - norm(r)*dot(v,v)/mu);
interval = 030;                         % time interval scenario to plot
dur      = 2*pi*sqrt(ah^3/mu)/3600/24;  % orbit propagation duration [days]
res      = 100;                         % # points to plot for orbit


%% Monte Carlo Simulation
%       Collect data through Monte Carlo and perform rudimentary
%       statistical analysis without sorting data into elliptical
%       and hyperbolic cases.

tot_iter = (length(int_min:int_res:int_max))*trials;
prog_bar = waitbar(0,'Monte Carlo Progress: 0%');
data(tot_iter).int = 0;
idx = 0;

for i = int_min:int_res:int_max
    rng(rng_val);
    for j = 1:trials
        idx = idx + 1;
        waitbar(idx/tot_iter,prog_bar,...
            ['Monte Carlo Progress:' num2str(idx/tot_iter*100) '%']);
        [err,df] = Auto_Test_Ext(r,v,mu,i,noise,start_day,v_count);
        data(idx).int = i;
        data(idx).df  = df;
        data(idx).orb = err(1:6);
        data(idx).is_elliptic = boolean(err(7));
    end
end

load handel
sound(y,Fs)

% Consolidate Data
%       Sort simulation data by elliptical and hyperbolic cases
%       then perform statistical analysis of mean and standard deviation.

init_cond.rng = rng_val;
init_cond.day = start_day;
init_cond.min = int_min;
init_cond.max = int_max;
init_cond.irs = int_res;
init_cond.vct = v_count;
init_cond.num = trials;
init_cond.rms = noise;
init_cond.mu  = mu;
init_cond.pos = r;
init_cond.vel = v;
init_cond.notes = notes;
init_cond.use_f = use_f;
init_cond.dur   = dur;
init_cond.res   = res;
init_cond.interval = interval;

[data_ell,data_hyp,mean_ell,std_ell,...
    mean_hyp,std_hyp] = Data_Consolidate(data,init_cond);

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


%% Plot Error Values

close all
Plot_Data(mean_ell, std_ell, 'Elliptic', use_f);
Plot_Data(mean_hyp, std_hyp, 'Hyperbolic', use_f);


%% Plot True Orbit


pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Red')
plot3(0,0,0,'bo')
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

pos_m = [];
vel_m = [];
for i = start_day:(interval/(v_count-1)/5):start_day+interval
    [rm,vm] = TimeProp_Universal_V2(r, v, mu, i);
    pos_m = [pos_m; rm];
    vel_m = [vel_m; vm];
end
plot3(pos_m(:,1), pos_m(:,2), pos_m(:,3),'LineWidth',2.5,'Color','Green')

pbaspect([1 1 1])
axis equal
grid on
hold off


%% Plot Sample Guess Orbit

reps = 5;

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