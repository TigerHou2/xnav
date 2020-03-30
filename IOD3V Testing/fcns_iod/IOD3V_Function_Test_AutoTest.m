function [err,df,r1,v1,v2,v3,r2,r3] = IOD3V_Function_Test_AutoTest(r,v,mu,dt,noise,start_day)
%AUTO_TEST Summary of this function goes here
%   dt: interval between measurements (days)
%   noise: velocity measurement accuracy (m/s)
%   start_day: starts all orbits x days after state given by r,v

% change the data collection point and separation
% will affect testing results
t1 = start_day;
t2 = t1 + dt;
t3 = t2 + dt;

[ri, vi] = TimeProp_Universal_V2(r, v, mu, t1);
[ ~, v2] = TimeProp_Universal_V2(r, v, mu, t2);
[re, ve] = TimeProp_Universal_V2(r, v, mu, t3);

v1 = vi;
v3 = ve;

% add noise to velocity 'measurements'
% assuming variance to be proportional to speed
%   don't know what the error order of magnitude is
%   so we will use this as a test to determine sensor requirements
v1 = v1 + randn(1,3) * noise / 1000;
v2 = v2 + randn(1,3) * noise / 1000;
v3 = v3 + randn(1,3) * noise / 1000;

% solve the 3V IOD problem
% [r1,r2,r3] = IOD3V_V2(v1,v2,v3,mu,'omega',[0,0,1],'prograde',true);
% [r1,r2,r3] = IOD3V_V2(v1,v2,v3,mu,'ordered',true);

r = IOD3V_Ext([v1;v2;v3],mu,'ordered',true);

r1 = r(1,:);
r2 = r(2,:);
r3 = r(3,:);

%================================

% we have been able to show r1,r2,r3 give the same orbits
% so it is no longer necessary to plot all three

[at,et,it,omgt,wt,ft] = Get_Orb_Params(ri,vi,mu);
[ a, e, i, omg, w, f] = Get_Orb_Params(r1,v1,mu);
[ ~, ~, ~,   ~, ~,fe] = Get_Orb_Params(re,ve,mu);

df = fe - ft;

err = [ (a-at), ...
        norm(e)-norm(et), ...
        (i-it), ...
        (omg-omgt), ...
        (w-wt), ...
        (f-ft), ...
        (a>0)];

end