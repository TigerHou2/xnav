clear;clc
Parameters

[filename,filepath] = uigetfile;
if not(filepath)
    disp('No file selected.')
    return
end
load(fullfile(filepath,filename))

options = {'Plot Statistics', 'Plot True Orbit', 'Plot Swarm'};
indx = listdlg('ListString',options,...
               'PromptString','Display Options',...
               'CancelString','No Plotting');

plot_stats = ismember(1,indx);
plot_orb   = ismember(2,indx);
plot_swarm = ismember(3,indx);

rng_val   = init_cond.rng;
start_day = init_cond.day;
int_min   = init_cond.min;
int_max   = init_cond.max;
int_res   = init_cond.irs;
v_count   = init_cond.vct;
trials    = init_cond.num;
noise     = init_cond.rms;
mu        = init_cond.mu;
r         = init_cond.pos;
v         = init_cond.vel;
notes     = init_cond.notes;
use_f     = init_cond.use_f;
dur       = init_cond.dur;
res       = init_cond.res;
interval  = init_cond.interval;

[r1,v1] = TimeProp_Universal_V2(r,v,mu,start_day);
[~,~,~,~,~,f] = Get_Orb_Params(r1,v1,mu);

f = msgbox({['File Name: ' filename] ...
             '' ...
             '========================' ...
             '- Setup -' ...
            ['Position: ' mat2str(r) ' km'] ...
            ['Velocity: ' mat2str(v) ' km/s'] ...
            ['G. Param: ' num2str(mu) ' km^3/s^2'] ...
            ['Description: ' notes] ...
             '' ...
             '========================' ...
             '- Test Conditions -' ...
            ['Noise: ' num2str(noise) ' m/s'] ...
            ['Obsv. Start Day: ' num2str(start_day)] ...
            ['Obsv. Start True Anomaly: ' num2str(rad2deg(f)) ' deg'] ...
             '' ...
             '========================' ...
             '- Monte Carlo Setup -' ...
            ['RNG Seed: ' num2str(rng_val)] ...
            ['Min. Interval: ' num2str(int_min) ' days'] ...
            ['Max. Interval: ' num2str(int_max) ' days'] ...
            ['Interval Resolution: ' num2str(int_res) ' days'] ...
            ['Measurements per Trial: ' num2str(v_count)] ...
            ['Trials per Interval: ' num2str(trials)] ...
            ''});
tx = findall(f, 'Type', 'Text');
tx.FontSize = 11;
tx.FontName = 'FixedWidth';
set(f, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.04, 0.5, 0.48]);

if plot_stats
    disp('Plotting statistics...')
    [~,~,mean_ell,std_ell,mean_hyp,std_hyp] = ...
        Data_Consolidate_V2(data,init_cond);
    Plot_Data(mean_ell, std_ell, 'Elliptic', use_f, filename);
    set(gcf,'Units','Normalized','OuterPosition',[0, 0.52, 0.5, 0.48]);
    Plot_Data(mean_hyp, std_hyp, 'Hyperbolic', use_f, filename);
    set(gcf,'Units','Normalized','OuterPosition',[0, 0.04, 0.5, 0.48]);
end

if plot_orb
    disp('Plotting true orbit...')
    pos = Get_Orb_Points(r,v,mu,res,dur,start_day);
    e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu * AU;

    figure('Name',['True Orbital Trajectory ' filename])
    
    hold on

    plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Red')
    plot3(0,0,0, 'bo')
    quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')

    pos_m = [];
    vel_m = [];
    for i = start_day:interval/5:start_day+2*interval
        [rm,vm] = TimeProp_Universal_V2(r, v, mu, i);
        pos_m = [pos_m; rm];
        vel_m = [vel_m; vm];
    end
    plot3(pos_m(:,1), pos_m(:,2), pos_m(:,3),'LineWidth',2,'Color','Green')

    view(cross(r,v))
    pbaspect([1 1 1])
    axis equal
    grid on
    hold off
    
    set(gcf,'Units','Normalized','OuterPosition',[0.5, 0.52, 0.5, 0.48]);
end

if plot_swarm
    swarm(trials).a   = [];
    swarm(trials).e   = [];
    swarm(trials).i   = [];
    swarm(trials).omg = [];
    swarm(trials).w   = [];
    swarm(trials).f   = [];
    data_size = size(data,2);
    for i = 1:trials
        for j = i:trials:data_size
            swarm(i).a   = [swarm(i).a   data(j).orb(1)];
            swarm(i).e   = [swarm(i).e   data(j).orb(2)];
            swarm(i).i   = [swarm(i).i   data(j).orb(3)];
            swarm(i).omg = [swarm(i).omg data(j).orb(4)];
            swarm(i).w   = [swarm(i).w   data(j).orb(5)];
            swarm(i).f   = [swarm(i).f   data(j).orb(6)];
        end
    end
    
    if use_f
        msg_x = 'True Anomaly Change (rad)';
        ser = [data(1:trials:end).df];
    else
        msg_x = 'Measurement Interval (days)';
        ser = [int_min:int_max];
    end
    
    figure('Name',['Population Data ' filename])
    for i = 1:trials
        subplot(2,3,1)
        hold on
        plot(ser, swarm(i).a)
        title('Semi-Major Axis')
        xlabel(msg_x)
        ylabel('Error (km)')
        grid on
        hold off
        
        subplot(2,3,2)
        hold on
        plot(ser, swarm(i).e)
        title('Eccentricity')
        xlabel(msg_x)
        ylabel('Error')
        grid on
        hold off

        subplot(2,3,3)
        hold on
        plot(ser, swarm(i).i)
        title('Inclination')
        xlabel(msg_x)
        ylabel('Error (rad)')
        grid on
        hold off

        subplot(2,3,4)
        hold on
        plot(ser, swarm(i).omg)
        title('Longitude of Ascending Node')
        xlabel(msg_x)
        ylabel('Error (rad)')
        grid on
        hold off

        subplot(2,3,5)
        hold on
        plot(ser, swarm(i).w)
        title('Argument of Periapse')
        xlabel(msg_x)
        ylabel('Error (rad)')
        grid on
        hold off

        subplot(2,3,6)
        hold on
        plot(ser, swarm(i).f)
        title('True Anomaly')
        xlabel(msg_x)
        ylabel('Error (rad)')
        grid on
        hold off
    end
end

disp(['Done!' newline])