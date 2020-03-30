function Data_Disp(filepath,plot)
%DATA_DISP Summary of this function goes here
%   Detailed explanation goes here

clc;
load(filepath)

[rng_val, dur, start_day, int_min, int_max, ...
 trials, noise, mu, r, v] = deal(init_cond{:});

[r1,v1] = TimeProp_Universal_V2(r,v,mu,start_day);
[~,~,~,~,~,f] = Get_Orb_Params(r1,v1,mu);

fprintf([ 'Initial Conditions\n' ...
          '========================\n' ...
          '- Setup -\n' ...
          'Position: ' mat2str(r) ' km\n' ...
          'Velocity: ' mat2str(v) ' km/s\n' ...
          'G. Param: ' num2str(mu) ' km^3/s^2\n' ...
          'Description: ' notes '\n\n' ...
          '========================\n' ...
          '- Test Conditions -\n' ...
          'Noise: ' num2str(noise) ' m/s\n' ...
          'Obsv. Start Day: ' num2str(start_day) '\n' ...
          'Obsv. Start True Anomaly: ' num2str(rad2deg(f)) ' deg\n\n' ...
          '========================\n' ...
          '- Monte Carlo Setup -\n' ...
          'RNG Seed: ' num2str(rng_val) '\n' ...
          'Min. Interval: ' num2str(int_min) ' days\n' ...
          'Max. Interval: ' num2str(int_max) ' days\n' ...
          'Trials per Interval: ' num2str(trials) '\n\n']);

if plot
    Plot_Data(mean_ell, std_ell, 'Elliptic', use_f);
    Plot_Data(mean_hyp, std_hyp, 'Hyperbolic', use_f);
    Plot_Data(mean_tot, std_tot, 'All Data', use_f);
    
    pos = [];
    vel = [];
    for i = 0:res:dur
        [r1,v1] = TimeProp_Universal_V2(r, v, mu, i);
        pos = [pos; r1];
        vel = [vel; v1];
    end
    e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu;
    figure('Name','True Orbit Trajectory')
    hold on
    plot3(pos(:,1), pos(:,2), pos(:,3),'LineWidth',1.5,'Color','Red')
    plot3(0,0,0, 'bo')
    e_vec = e_vec*AU;
    quiver3(0,0,0,e_vec(1),e_vec(2),e_vec(3),'Color','Red')
    pbaspect([1 1 1])
    axis equal
    grid on
    hold off
end
  
end

