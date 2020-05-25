clear;clc

[filename1,filepath] = uigetfile;
if not(filepath)
    disp('No file selected.')
    return
end
load(fullfile(filepath,filename1))

[data_ell_one,data_hyp_one,...
 mean_ell_one,std_ell_one,...
 mean_hyp_one,std_hyp_one] = Data_Consolidate_Additive(data,init_cond);

[filename2,filepath] = uigetfile;
if not(filepath)
    disp('No file selected.')
    return
end
load(fullfile(filepath,filename2))

[data_ell_two,data_hyp_two,...
 mean_ell_two,std_ell_two,...
 mean_hyp_two,std_hyp_two] = Data_Consolidate_Additive(data,init_cond);

idx = 2; % we will use true anomaly interval measurements instead of days
msg_x = 'True Anomaly Change (rad)';
ylabels = { 'Error (km)', 'Error', ...
            'Error (rad)', 'Error (rad)', ...
            'Error (rad)', 'Error (rad)'};
titles  = { 'Semi-Major Axis', 'Eccentricity', ...
            'Inclination', 'Longitude of Ascending Node', ...
            'Argument of Periapse', 'True Anomaly'};

figure(300)
for i = 1:6
    if isempty(mean_ell_one) || isempty(std_ell_one) || ...
       isempty(mean_ell_two) || isempty(std_ell_two)
        disp(['No data available for elliptic orbits' newline])
        return
    end
    subplot(2,3,i)
    hold on
    errorbar(mean_ell_one(:,idx),mean_ell_one(:,i+2),...
                std_ell_one(:,i+2),'Color','r')
    errorbar(mean_ell_two(:,idx),mean_ell_two(:,i+2),...
                std_ell_two(:,i+2))
    hold off
    title(titles{i})
    xlabel(msg_x)
    ylabel(ylabels{i})
    grid on
    legend('Dataset 1','Dataset 2')
end

figure(301)
for i = 1:6
    if isempty(mean_hyp_one) || isempty(std_hyp_one) || ...
       isempty(mean_hyp_two) || isempty(std_hyp_two)
        disp(['No data available for hyperbolic orbits' newline])
        return
    end
    subplot(2,3,i)
    hold on
    errorbar(mean_hyp_one(:,idx),mean_hyp_one(:,i+2),...
                std_hyp_one(:,i+2),'Color','r')
    errorbar(mean_hyp_two(:,idx),mean_hyp_two(:,i+2),...
                std_hyp_two(:,i+2))
    hold off
    title(titles{i})
    xlabel(msg_x)
    ylabel(ylabels{i})
    grid on
    legend('Dataset 1','Dataset 2')
end
    