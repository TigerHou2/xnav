%% nbody_err_test.m
%
% Author: Tiger Hou
% 
% Description:
%   Runs n-body simulations to determine the error reduction rate for
%   local interpolation numerical methods such as RK4.
%
%% Setup

% clear console and workspace
close all hidden
clear;clc

% load n-body case
case_select = 15;
sub_case = 11;
nbody_test_cases;

% order reduction for gravity model
rv0 = [r0;v0];

% define various timestep resolutions
dtvect = [dt*2^7, dt*2^6, dt*2^5, dt*2^4, dt*2^3, dt*2^2, dt*2, dt];


%% Simulation

rv_end = zeros(6,size(rv0,2),length(dtvect));

for j = 1:length(dtvect)
    
    rv_end(:,:,j) = rk4(@gravity,rv0,T,dtvect(j),m,G);
    
end

filename = strcat('Err_Data\Case_',num2str(case_select));
save(filename, 'rv_end', 'dtvect', 'case_select', 'rv0');


%% Error Plots

errs = zeros(size(rv0,2),length(dtvect)-1);

subplot_rows = floor(sqrt(size(rv0,2)));
subplot_cols = ceil((size(rv0,2))/subplot_rows);
length_baseline = mean(vecnorm(r0,2,1));

for j = 1:length(dtvect)-1
    
    for k = 1:size(rv0,2)
        
        errs(k,j) = norm(rv_end(1:3,k,j)-rv_end(1:3,k,end)) / ...
                    length_baseline;
    
    end
    
end

f = figure(1);

for k = 1:size(rv0,2)
    
    loglog(dtvect(1:end-1),errs(k,:), ...
           '.-', 'markersize', 32, 'linewidth', 2)
    hold on
    
end

hold off

% title(['Error Convergence of the ' num2str(size(rv0,2)) '-Body System'], ...
%       'fontsize', 40, 'interpreter' , 'latex')
xlabel('$\Delta t$', 'fontsize', 50, 'interpreter' , 'latex')
ylabel('Normalized Error', 'fontsize', 50, 'interpreter' , 'latex')
legend(lgd, 'fontsize', 36, 'interpreter' , 'latex', 'Location', 'Southeast')

grid(gca,'minor')
grid on
set( gca, 'Color', [1 1 1] )
set( gca, 'fontsize', 50, 'ticklabelinterpreter', 'latex' )
pbaspect([1,1,1])

set(f, 'PaperPositionMode', 'manual')
set(f, 'Color', [1 1 1])
set(f, 'PaperUnits', 'centimeters')
set(f, 'PaperSize', [30 30])
set(f, 'Units', 'centimeters' )
set(f, 'Position', [0 0 30 30])
set(f, 'PaperPosition', [0 0 30 30])

lim_y = ylim;
lb = floor(log10(lim_y(1)));
ub =  ceil(log10(lim_y(2)));
set( gca, 'YTick', ...
     logspace(lb,ub,ub-lb+1));
% set( gca, 'YTick', ...
%      logspace(lb,ub-1,(ub-lb+1)/2));
 
lim_x = xlim;
lb = floor(log10(lim_x(1)));
ub =  ceil(log10(lim_x(2)));
set( gca, 'XTick', ...
     logspace(lb,ub,ub-lb+1));

svnm = strcat('.\Figures\case', num2str(case_select), '_error');
print( '-dpng', svnm, '-r200' )

