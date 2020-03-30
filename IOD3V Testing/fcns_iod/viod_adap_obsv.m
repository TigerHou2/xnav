function [R_est,R,V,T] = viod_adap_obsv(r,v,mu,start_day,obsv_cap,da,noise,x)
%VIOD_ADAP_OBSV Summary of this function goes here
%   Detailed explanation goes here

R = zeros(obsv_cap,3); % position (ground truth)
V = zeros(obsv_cap,3); % velocity (measured, w/ error)
T = zeros(obsv_cap,1); % time since periapsis

% iterate

next_t = start_day;

% symbolic math setup

syms k ff real
assume(k,'positive')

for i = 1:obsv_cap
    
    % propagate time, collect next measurement
    
    T(i) = next_t;
    [R(i,:),V(i,:)] = TimeProp_V3(r,v,mu,T(i),x);
    
    % store ground truth data to use for propagating
    
    rt = R(i,:);
    vt = V(i,:);
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    V(i,:) = V(i,:) + n_vec;
%     vv = V(i,:);
    
    % propagate to find proper dt
    
    dt = get_dt(rt,vt,mu,da,k,ff);
    next_t = T(i) + dt;
    
%     step = 0;
%     step_size = 0.5;
%     
%     angle_prev = 0;
%     
%     while true
%         
%         step = step + step_size;
%         [~,v_next] = TimeProp_Universal_V2(rt,vt,mu,step);
%         
%         angle = acos(dot(vv,v_next)/norm(vv)/norm(v_next));
%         
%         if angle > da
%             next_t = T(i) + step;
%             break
%         end
%         
%         if abs(angle-da) > abs(angle-angle_prev)
%             step_size = step_size * 1.2;
%         else
%             step_size = step_size * 0.8;
%         end
%         
%         angle_prev = angle;
%         
%     end
    
end

% perform IOD

C = combnk(1:obsv_cap,2);
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
A = spalloc(comb_n*4,obsv_cap*2,7*obsv_cap*(obsv_cap-1));
B = zeros(comb_n*4,1);

R_est = IOD3V_V3(V,A,B,C,mu,'ordered',true);

end


%% Auxiliary Functions

function dt = get_dt(r,v,mu,da,k,ff)

[a,e,~,~,~,f] = Get_Orb_Params(r,v,mu);
e = norm(e);

curr_dir = [-sin(f); e+cos(f)];
rotn_mat = [cos(da) -sin(da); sin(da) cos(da)];
next_dir = rotn_mat * curr_dir;

eqn = k*[-sin(ff);e+cos(ff)] == next_dir;
soln = solve(eqn,[k,ff]);

f_next = double(soln.ff);

E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2));
E_next = 2*atan(sqrt((1-e)/(1+e))*tan(f_next/2));

if E < 0
    E = E + 2*pi;
end

if E_next < 0
    E_next = E_next + 2*pi;
end

M = E - e*sin(E);
M_next = E_next - e*sin(E_next);

dt = (M_next-M) / sqrt(mu/a^3) / 3600 / 24;

end