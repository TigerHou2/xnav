function [r,runTime] = monteScatter(v,mu,noise,rRef,itr,idx,bias)
%MONTESCATTER Generates a scatter plot od OD errors from Monte Carlo sims
%
% Author: 
%   Tiger Hou
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


tic

latexify
r2_vect = nan(itr,3);

if ~exist('idx','var')
    idx = 2;
end
if ~exist('bias','var')
    bias = zeros(1,3);
end

for i = 1:itr
    
    v_noisy = addnoise(v,noise,bias);
    
    rEst = hodo_od(v_noisy, mu);
    
    if i == 1
        r = rEst; % grab sample data to return
    end
    
    r2 = rRef(idx,:) - rEst(idx,:);
    
    r2_vect(i,:) = r2;
    
end

scatter3(r2_vect(:,1),r2_vect(:,2),r2_vect(:,3),'kx');
axis equal
xlabel('$\delta_x, km$')
ylabel('$\delta_y, km$')
zlabel('$\delta_z, km$')
latexify

runTime = toc;

end

