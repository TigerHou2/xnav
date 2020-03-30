function [R_est,R,V,T] = viod_eqsp(r,v,mu,start_day,end_day,obsv_count,noise,x)
%VIOD_EQSP Summary of this function goes here
%   Detailed explanation goes here

R = zeros(obsv_count,3);
V = zeros(obsv_count,3);
T = linspace(start_day,end_day,obsv_count);

for i = 1:obsv_count
    
    % propagate
    
    [R(i,:),V(i,:)] = TimeProp_V3(r,v,mu,T(i),x);
    
    % add noise
    
    n_vec = randn(1,3);
    n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
    V(i,:) = V(i,:) + n_vec;
    
end

C = combnk(1:i,2);           % permutations of two velocity vectors i,j
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
A = spalloc(comb_n*4,i*2,7*i*(i-1)); % LHS matrix of eq. 36, see line 6
B = zeros(comb_n*4,1);               % RHS matrix of eq. 36, see line 6

R_est = IOD3V_V3(V,A,B,C,mu,'ordered',true);

end

