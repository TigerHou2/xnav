function Plot_Data(mean,std,name,use_f)
%PLOT_DATA Summary of this function goes here
%   Detailed explanation goes here

figure('Name',name)

if isempty(mean) || isempty(std)
    disp(['No data available for "' name '"' newline])
    return
end

stat_mean = mean;
stat_std  = std;

if use_f
    msg_x = 'True Anomaly Change Between Measurements (rad)';
else
    msg_x = 'Time Between Measurements (days)';
end

subplot(2,3,1)
errorbar(stat_mean(:,1),stat_mean(:,2),stat_std(:,2))
title('Semi-Major Axis')
xlabel(msg_x)
ylabel('Error (km)')
grid on

subplot(2,3,2)
errorbar(stat_mean(:,1),stat_mean(:,3),stat_std(:,3))
title('Eccentricity')
xlabel(msg_x)
ylabel('Error')
grid on

subplot(2,3,3)
errorbar(stat_mean(:,1),stat_mean(:,4),stat_std(:,4))
title('Inclination')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

subplot(2,3,4)
errorbar(stat_mean(:,1),stat_mean(:,5),stat_std(:,5))
title('Longitude of Ascending Node')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

subplot(2,3,5)
errorbar(stat_mean(:,1),stat_mean(:,6),stat_std(:,6))
title('Argument of Periapse')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

subplot(2,3,6)
errorbar(stat_mean(:,1),stat_mean(:,7),stat_std(:,7))
title('True Anomaly')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

end