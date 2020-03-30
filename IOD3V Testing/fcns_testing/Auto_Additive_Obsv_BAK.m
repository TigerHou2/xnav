function [err,df,R,V] = Auto_Additive_Obsv_BAK(r,v,mu,dt,noise,start_day,reps)
%AUTO_TEST_EXT Summary of this function goes here
%   dt: interval between measurements (days)
%   noise: velocity measurement accuracy (m/s)
%   start_day: starts all orbits x days after state given by r,v

% change the data collection point and separation
% will affect testing results
t = zeros(1,reps);
for i = 1:reps
    t(i) = start_day + (i-1)*dt;
end

% use Chebyshev points to define data collection intervals
% t = zeros(1,reps);
% dt = DT / (reps-1);
% end_day = start_day + DT;
% for i = 1:reps
%     t(i) = (start_day+end_day)/2 + (start_day-end_day)/2*cos(pi*(i-1)/(reps-1));
% end

R = zeros(reps,3);
V = zeros(reps,3);

for i = 1:reps
    [R(i,:),V(i,:)] = TimeProp_Universal_V2(r,v,mu,t(i));
end

ri = R(1,:);
vi = V(1,:);
rf = R(end,:);
vf = V(end,:);

% add noise to velocity 'measurements'
% assuming variance to be proportional to speed
%   don't know what the error order of magnitude is
%   so we will use this as a test to determine sensor requirements
for i = 1:reps
    n_vec = randn(1,3);
    n_vec = n_vec / norm(n_vec) * noise / 1000;
    V(i,:) = V(i,:) + n_vec;
end

% solve the 3V IOD problem
R = IOD3V_Ext(V,mu,'ordered',true);
r1 = R(1,:);
v1 = V(1,:);

[at,et,it,omgt,wt,ft] = Get_Orb_Params(ri,vi,mu);
[ a, e, i, omg, w, f] = Get_Orb_Params(r1,v1,mu);
[ ~, ~, ~,   ~, ~,ff] = Get_Orb_Params(rf,vf,mu);

df = ff - ft;

err = [ (a-at), ...
        norm(e)-norm(et), ...
        (i-it), ...
        (omg-omgt), ...
        (w-wt), ...
        (f-ft), ...
        (a>0)];

end

