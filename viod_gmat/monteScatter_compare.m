function monteScatter_compare(v,mu,noise,rRef,itr,idx,bias)
%MONTESCATTER_COMPARE Generates a scatter plot od OD errors from Monte Carlo sims
%
% Author: 
%   Tiger Hou
%
% Description:
%   A modified version of monteScatter which compares two methods of viod.
%
% Arguments:
%   v     - km/s    , N-by-3 array of velocity measurements
%   mu    - km^3/s^2, scalar value of gravitational parameter
%   rRef  - km/s    , N-by-3 array of ground truth position measurements
%   itr   - nd      , number of Monte Carlo sims to run
%   idx   - nd      , index of measurement to calculate error
%
% Optional Arguments:
%   bias  - m/s     , 1-by-3 array of velocity measurement biases
%
% Notes: 
%   For arguments v and rRef, each row corresponds to the velocity and
%   position at a measurement location in orbit. 
%
% Outputs:
%   Scatter plot of OD errors.

figure;
latexify(45,18)

r2_vect_hodo = nan(itr,3);
r2_vect_orig = nan(itr,3);

if ~exist('idx','var')
    idx = 2;
end
if ~exist('bias','var')
    bias = zeros(1,3);
end

for i = 1:itr
    
    v_noisy = addnoise(v,noise,bias);
    
    rEst_hodo = hodo_od(v_noisy, mu);
    rEst_orig = viod(v_noisy,mu);
    
    r2_hodo = rRef(idx,:) - rEst_hodo(idx,:);
    r2_orig = rRef(idx,:) - rEst_orig(idx,:);
    
    r2_vect_hodo(i,:) = r2_hodo;
    r2_vect_orig(i,:) = r2_orig;
    
end

subplot(1,3,1)
hold on
scatter(r2_vect_hodo(:,1),r2_vect_hodo(:,2),'r+');
scatter(r2_vect_orig(:,1),r2_vect_orig(:,2),'kx');
hold off
axis equal
xlabel('$\delta_x, km$')
ylabel('$\delta_y, km$')
legend('Hodograph VIOD','Original VIOD')

subplot(1,3,2)
hold on
scatter(r2_vect_hodo(:,1),r2_vect_hodo(:,3),'r+');
scatter(r2_vect_orig(:,1),r2_vect_orig(:,3),'kx');
hold off
axis equal
xlabel('$\delta_x, km$')
ylabel('$\delta_z, km$')
legend('Hodograph VIOD','Original VIOD')

subplot(1,3,3)
hold on
scatter(r2_vect_hodo(:,2),r2_vect_hodo(:,2),'r+');
scatter(r2_vect_orig(:,2),r2_vect_orig(:,2),'kx');
hold off
axis equal
xlabel('$\delta_y, km$')
ylabel('$\delta_z, km$')
legend('Hodograph VIOD','Original VIOD')

end

