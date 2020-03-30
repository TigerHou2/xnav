%% Fixes

% currently we think it's producing a hyperbolic orbit with negative a and
% therefore imaginary E

% plot velocity vectors to check if collinear using rng(1) on day 22
%   > the error was not due to collinearity. velocity vectors were already
%     imaginary prior to input to IOD3V.m
%   > error still comes from TimeProp.m, therefore it's an error with
%     the hyperbolic orbit
%   > solved by implementing universal solver TimeProp_Universal.m
%     10/11/2019

% check days 12, 15, 16 for eccentricity and arg. of pe. error
% check results of changing seed for rng
%   > it appears that all anomalies still occur at the same days
%   > removing absolute value confirms the mean error trends to zero
%   > FALSE - it was because of a bug in the code causing rng to be
%   re-seeded every time

% check semi-major axis error spike at 55 days
% check Neptune 2000
%   especially true anomaly, arg. of pe, ecc


%% Initial Conditions
%       Set up the Monte Carlo environment and load data.

close all
clear;clc

Parameters

% ======================================
% Orbital Parameters
% ======================================

r  = [0.71,  0, 0.71] * AU; % position [km]
v  = [0   , 41.34, 0];      % velocity [km/s]
mu = sun.mu;                % gravitational parameter [km^3/s^2]

% additional remarks
notes = 'Earth-Neptune Transfer Orbit';


% ======================================
% Testing Conditions
% ======================================

noise     = 6;              % sensor noise in each axis [m/s]
start_day = 200;            % measurement start day [days]


% ======================================
% Monte Carlo Parameters
% ======================================

rng_val = 86583;            % rng seed
int_min = 10;               % minimum interval [days]
int_max = 60;               % maximum interval [days]
trials  = 5;              % trials per interval (approx. 0.2 seconds per trial)


% ======================================
% Data Settings
% ======================================

save_data = false;          % save data?
use_f     = true;           % use true anomaly instead of days on x-axis


% ======================================
% Orbit Visualization
% ======================================

interval = 100;             % time interval scenario to plot
dur      = 6000;            % orbit propagation duration [days]


% ======================================
% Initial Condition Compilation
% ======================================
init_cond = {rng_val, dur, start_day, int_min, int_max, ...
             trials, noise, mu, r, v};


%% Monte Carlo Simulation
%       Collect data through Monte Carlo and perform rudimentary
%       statistical analysis without sorting data into elliptical
%       and hyperbolic cases.

mean_tot = [];
std_tot  = [];
all_data = [];

for i = int_min:int_max
    data = [];
    rng(rng_val);
    for j = 1:trials
        disp([i j])
        [err,df] = Auto_Test(r,v,mu,i,noise,start_day);
        if use_f
            x = df;
            z = i;
        else
            x = i;
            z = df;
        end
        all_data = [all_data; [x err z]];
        data = [data; [x err z]];
    end
    mean_tot = [mean_tot; [x,mean(data(:,2:end))]];
    std_tot  = [ std_tot; [x, std(data(:,2:end))]];
end


%% Consolidate Data
%       Sort simulation data by elliptical and hyperbolic cases
%       then perform statistical analysis of mean and standard deviation.

data_ell = all_data( boolean(all_data(:,end-1)),:);
data_hyp = all_data(~boolean(all_data(:,end-1)),:);

% if true anomaly is used for the x-axis, the column containing information
% on the time interval changes. 
if use_f
    day_info = size(data,2);
else
    day_info = 1;
end

curr_day = int_min;
temp = [];
mean_ell = [];
std_ell = [];

for i = 1:length(data_ell)
    if not (data_ell(i,day_info) == curr_day)
        if size(temp,1) <= 1
            mean_ell = [mean_ell; temp];
            std_ell  = [std_ell ; temp];
        else
            mean_ell = [mean_ell; mean(temp)];
            std_ell  = [std_ell ; std(temp)];
        end
        temp = [];
        curr_day = curr_day + 1;
    end
    temp = [temp; data_ell(i,:)];
end

curr_day = int_min;
temp = [];
mean_hyp = [];
std_hyp = [];

for i = 1:length(data_hyp)
    if not (data_hyp(i,day_info) == curr_day)
        if size(temp,1) <= 1
            mean_hyp = [mean_hyp; temp];
            std_hyp  = [std_hyp ; temp];
        else
            mean_hyp = [mean_hyp; mean(temp)];
            std_hyp  = [std_hyp ; std(temp)];
        end
        temp = [];
        curr_day = curr_day + 1;
    end
    temp = [temp; data_hyp(i,:)];
end

if save_data
    data_files = dir('Data/Data-*.mat');
    if size(data_files,1) == 0
        label = 1;
    else
        data_names = [data_files.name];
        label = max(str2num(cell2mat(regexp(data_names,'\d+','match')')))+1;
    end
    save(['Data/Data-' num2str(label,'%03.f')]);
end

%% Plot Error Values

close all
Plot_Data(mean_ell, std_ell, 'Elliptic', use_f);
Plot_Data(mean_hyp, std_hyp, 'Hyperbolic', use_f);
Plot_Data(mean_tot, std_tot, 'All Data', use_f);

%% Plot True Orbit

pos = [];
vel = [];
for i = 0:res:dur
    [r1,v1] = TimeProp_Universal_V2(r, v, mu, i);
    pos = [pos; r1];
    vel = [vel; v1];
end
e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

hold on

plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Red')
plot3(0,0,0, 'bo')
quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

pos_m = [];
vel_m = [];
for i = start_day:interval:start_day+2*interval
    [rm,vm] = TimeProp_Universal_V2(r, v, mu, i);
    pos_m = [pos_m; rm];
    vel_m = [vel_m; vm];
end
plot3(pos_m(:,1), pos_m(:,2), pos_m(:,3),'LineWidth',1.5,'Color','Green')

pbaspect([1 1 1])
axis equal
grid on
hold off

%% Plot Sample Guess Orbit

reps = 3;

for k = 1:reps

    pos = [];
    vel = [];

    [~,~,r_IOD,v_IOD,v2,v3,r2,r3] = Auto_Test(r,v,mu,interval,noise,start_day);
    v1 = v_IOD;
    r1 = r_IOD;

    for i = 0:res:dur
        [r_,v_] = TimeProp_Universal_V2(r_IOD, v_IOD, mu, i);
        pos = [pos; r_];
        vel = [vel; v_];
    end
    e_vec = ((dot(v_IOD,v_IOD)-mu/norm(r_IOD))*r_IOD-dot(r_IOD,v_IOD)*v_IOD)/mu;
    e_vec = e_vec*AU;

    hold on

    plot3(pos(:,1), pos(:,2), pos(:,3),'Color','Blue')

    pbaspect([1 1 1])
    axis equal
    grid on
    hold off
    
end