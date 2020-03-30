function [R_est,R,V,T] = viod_adap_pred(r,v,mu,start_day,obsv_cap,da,noise,x)
%VIOD_ADAP_PRED Summary of this function goes here
%   Detailed explanation goes here

% perform first three observations

R =     zeros(obsv_cap,3); % position (ground truth)
R_est = zeros(obsv_cap,3); % position (IOD estimate)
V =     zeros(obsv_cap,3); % velocity (measured, w/ error)
T =     zeros(obsv_cap,1); % time since periapsis


C = combnk(1:3,2);
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
A = spalloc(comb_n*4,6,42);
B = zeros(comb_n*4,1);


% iterate for first three measurements
% assumes true orbit is not hyperbolic, so if IOD gives hyperbolic
% then the program increases the initial delta-t
% and reworks the first three measurements.
%
% we're kind of cheating the system by going back in time, 
% but on an actual mission we would just take a lot of measurements
% and use the start, midpoint, and end.

delta_t = 5;

while true
    
    for i = 1:3

        % propagate forward in orbit by delta-t, get position and velocity

        T(i) = start_day + (i-1) * delta_t;
        [R(i,:),V(i,:)] = TimeProp_V3(r,v,mu,T(i),x);

        % add noise to simulate velocity measurements w/ error

        n_vec = randn(1,3);
        n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
        V(i,:) = V(i,:) + n_vec;

    end

    % perform IOD using velocity measurements to get IOD positions
    
    RR = IOD3V_V3(V(1:3,:),A,B,C,mu,'ordered',true);
    
    % get geometry of IOD orbit from IOD position and velocity measurements
    
    rr = RR(1,:);
    vv = V(1,:);
    ah_est = norm(rr)/(2 - norm(rr)*dot(vv,vv)/mu);
    
    % if IOD orbit is hyperbolic, redo the above with a larger delta-t
    
    if ah_est <= 0
        delta_t = delta_t + 1;
    else
        break
    end
    
end


% iterate for observation #4 and above

for i = 3:obsv_cap
    
    C = combnk(1:i,2);
    comb_n = size(C,1);
    % using sparse matrix to reduce memory consumption
    A = spalloc(comb_n*4,i*2,7*i*(i-1));
    B = zeros(comb_n*4,1);
    
    
    % perform IOD using velocity measurements to get IOD positions
    
    RR = IOD3V_V3(V(1:i,:),A,B,C,mu,'ordered',true);
    R_est(i,:) = RR(i,:);
    
    
    % get IOD orbit parameters
    
    rr = R_est(i,:);
    vv = V(i,:);
    
    
    % propagate to find proper dt
    
    step = 0;
    step_size = 1;
    
    diff_m1 = 0;
    
    while true
        
        % simulate propagating forward using current orbit estimate
        step = step + step_size;
        [~,vt] = TimeProp_V3(rr,vv,mu,step,x);
        
        diff = acos(dot(vv,vt)/norm(vv)/norm(vt));
        
        if diff > da
            delta_t = step;
            break
        end
        
        if abs(diff-da) > abs(diff-diff_m1)
            step_size = step_size * 1.5;
        else
            step_size = step_size * 0.8;
        end
        
        diff_m1 = diff;
        
    end
    
    % propagate time, collect next measurement
    
    T(i+1) = T(i) + delta_t;
    [R(i+1,:),V(i+1,:)] = TimeProp_V3(r,v,mu,T(i+1),x);
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    V(i+1,:) = V(i+1,:) + n_vec;
    
end

% clean up the (i+1)th term

R = R(1:obsv_cap,:);
V = V(1:obsv_cap,:);

end

