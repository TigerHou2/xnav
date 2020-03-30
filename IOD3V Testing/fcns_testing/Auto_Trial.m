function [err,df,R,V] = Auto_Trial(R,V,A,B,C,mu,noise,reps)
%AUTO_TRIAL Performs VIOD under noise and outputs error of six orbital
% elements of IOD orbit compared to true orbital elements
%
% Created:  01/30/2020
% Modified: 02/03/2020
% Author:   Tiger Hou
%
% Changes: 
%   Previous Version: Auto_Test_Ext.m
%       - Now calls IOD3V_V3.m instead of IOD3V_Ext.m, which has improved
%           by preallocating sparse matrix memory for calculations
%           to reduce computation time.
%   
% Dependencies (in call order):
%   IOD3V_V3.m
%   Get_Orb_Params.m
%
% Arguments:
%   R:      [km]        N-by-3 matrix containing N position measurements
%   V:      [km/s]      N-by-3 matrix containing N velocity measurements
%   A:      [unitless]  [2*N*(N-1)]-by-[2*N] sparse matrix of linear system
%   B:      [unitless]  [2*N*(N-1)]-by-1 matrix of linear system
%   C:      [unitless]  [N*(N-1)/2]-by-2 matrix of combinations of
%                           two numbers selected from 1:N, non-repeating
%   mu:     [km^3/s^2]  gravitational parameter of orbiting body
%   noise:  [m/s]       standard deviation of noise (m/s)
%   reps:   [unitless]  number of points sampled in orbit

ri = R(1,:);
vi = V(1,:);
rf = R(end,:);
vf = V(end,:);

% add noise to velocity 'measurements'
% assuming variance to be proportional to speed
%   don't know what the error order of magnitude is
%   so we will use this as a test to determine sensor requirements
n_vec = randn(reps,3);
n_vec = n_vec ./ vecnorm(n_vec,2,2) * noise / 1000;
V = V + n_vec;

% solve the 3V IOD problem
R = IOD3V_V3(V,A,B,C,mu,'ordered',true);
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

