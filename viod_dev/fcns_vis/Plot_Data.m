function Plot_Data(mean,std,name,use_f,varargin)
%PLOT_DATA Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addOptional(p,'filename','',@ischar);
parse(p,varargin{:});

figure('Name',[name ' ' p.Results.filename])

if isempty(mean) || isempty(std)
    disp(['No data available for "' name '"' newline])
    return
end

if use_f
    msg_x = 'True Anomaly Change (rad)';
    idx = 2;
else
    msg_x = 'Measurement Interval (days)';
    idx = 1;
end

subplot(2,3,1)
errorbar(mean(:,idx),mean(:,3),std(:,3))
title('Semi-Major Axis')
xlabel(msg_x)
ylabel('Error (km)')
grid on

subplot(2,3,2)
errorbar(mean(:,idx),mean(:,4),std(:,4))
title('Eccentricity')
xlabel(msg_x)
ylabel('Error')
grid on

subplot(2,3,3)
errorbar(mean(:,idx),mean(:,5),std(:,5))
title('Inclination')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

subplot(2,3,4)
errorbar(mean(:,idx),mean(:,6),std(:,6))
title('Longitude of Ascending Node')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

subplot(2,3,5)
errorbar(mean(:,idx),mean(:,7),std(:,7))
title('Argument of Periapse')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

subplot(2,3,6)
errorbar(mean(:,idx),mean(:,8),std(:,8))
title('True Anomaly')
xlabel(msg_x)
ylabel('Error (rad)')
grid on

end