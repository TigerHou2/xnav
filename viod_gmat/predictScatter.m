function [diff_fp] = predictScatter(num,V,mu,noise,tspan,dt,rf)
%PREDICTSCATTER Compares methods of estimating future position after viod
%
% Author:
%   Tiger Hou
%

% clean up unused variable space
V = vclean(V);

% convert seconds to days for propagation
tspan = tspan/3600/24;
dt = dt/3600/24;

diff_fp = nan(num,3); % first measurement point propagation
diff_cm = nan(num,3); % mean value of propagation from all points
diff_mr = nan(num,3); % mirror of fp about cm

for i = 1:num
    
    % add noise
    Vn = addnoise(V,noise);
    Re = viod(Vn,mu);

    % find each position estimate
    r_vect = nan(size(Re));
    for j = 1:size(Re,1)
        r_vect(j,:) = propagate(Re(j,:),Vn(j,:),mu,tspan-(j-1)*dt);
    end
    r_vect = rf - r_vect;

    % use the three methods
    diff_fp(i,:) = r_vect(1,:);
    diff_cm(i,:) = mean(r_vect);
    
end

diff_fp = diff_fp;
diff_cm = diff_cm;
diff_mr = 2 * diff_cm - diff_fp;

figure;
latexify
hold on
scatter3(diff_fp(:,1),diff_fp(:,2),diff_fp(:,3),'kx')
scatter3(diff_cm(:,1),diff_cm(:,2),diff_cm(:,3),'r+')
scatter3(diff_mr(:,1),diff_mr(:,2),diff_mr(:,3),'bo')
hold off
grid(gca,'minor')
grid on
xlabel('x axis estimate error, $km$')
ylabel('y axis estimate error, $km$')
zlabel('z axis estimate error, $km$')
title('Future Position Estimate Error')
legend('First Point','Mean Value','First Point Mirrored')
latexify

end

