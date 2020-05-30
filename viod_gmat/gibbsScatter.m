function [r_vect] = gibbsScatter(R,mu,dt,rf,v_noisy,dt_meas)
%GIBBSSCATTER Creates scatter of positions from solving Gibb's problem
%
% Author: 
%   Tiger Hou
%

R = vclean(R);
v_noisy = vclean(v_noisy);
r_vect = nan(size(R,1)-3,3);
r_vect_fromV = nan(size(R));

for i = 1:size(R,1)-3
    rr = R(1:end-i,:);
    v1 = gibbs(rr,mu);
    r_vect(i,:) = propagate(rr(1,:),v1,mu,dt);
end
r_vect = rf - r_vect;

for i = 1:size(R,1)
    r_vect_fromV(i,:) = propagate(R(i,:),v_noisy(i,:),mu,dt-(i-1)*dt_meas);
end
r_vect_fromV = rf - r_vect_fromV;

figure;
latexify
hold on
plot3(r_vect(:,1),r_vect(:,2),r_vect(:,3),'r-o')
plot3(r_vect_fromV(:,1),r_vect_fromV(:,2),r_vect_fromV(:,3),'k--o')
scatter3(r_vect(1,1),r_vect(1,2),r_vect(1,3),'ro','Filled')
scatter3(r_vect_fromV(1,1),r_vect_fromV(1,2),r_vect_fromV(1,3),'ko','Filled')
hold off
grid(gca,'minor')
grid on
xlabel('x axis estimate error, $km$')
ylabel('y axis estimate error, $km$')
zlabel('z axis estimate error, $km$')
title('Future Position Estimate Error')
legend('Gibbs','Velocity','Location','Best')
latexify

end

