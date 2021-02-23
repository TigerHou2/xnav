%% initialization
close all
clear;clc

% profile on

addpath('../fcns_circ')
addpath('../fcns_orb')

% define orbit
mu = 1;
a = 1e5;
i = deg2rad(0);
omg = deg2rad(0);
w = deg2rad(0);

% define true solution
e = 0.7;
f0 = deg2rad(45);
dM = deg2rad(90);

OPT = [f0,e,dM];

params = [a e i omg w 0];

% define pulsars
P = [0 -1 1;... % pulsar 1
     1 0 1;... % pulsar 2
     -1 0 3]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% number of measurements
numObsv = 7;

% clumped observations?
% ** clumped = 1 means the first 'numObsv' observations are all performed
%    on the first pulsar before moving on to the next. 
%    Conversely, clumped = 0 means for each observation we switch to the
%    next pulsar.
clumped = 1;

% calculate mean and true anomalies of measurements
E0 = 2*atan(sqrt((1-e)/(1+e))*tan(f0/2));
M0 = E0 - e*sin(E0);
M = M0 + linspace(0,dM,size(P,2)*numObsv);
if clumped == 1
    M = reshape(M,numObsv,size(P,2))';
else
    M = reshape(M,size(P,2),numObsv);
end
E = kepler(M,e);
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

% calculate time of measurements
t_true = (M-M0)*sqrt(a^3/mu);

% set random epoch
t_offset = 3.5; % just has to be a number unrelated to anything else
period = 2*pi * sqrt(a^3/mu);
t_offset = mod(t_offset,period);
t_meas = t_true + t_offset;

% perform observations
r = nan(size(f,1),size(f,2),3);
v = nan(size(f,1),size(f,2),3);
obsv = nan(size(f));
pulsar = P;
for i = 1:size(f,1)
    for j = 1:size(f,2)
        params(6) = f(i,j);
        [r(i,j,:),v(i,j,:)] = Get_Orb_Vects(params,mu);
        noise = 0;
%         noise = v(i,j,:) / norm(reshape(v(i,j,:),3,1)) ...
%                          * normrnd(0,2e-6);
        vtemp = v(i,j,:) + noise;
        obsv(i,j) = P(:,i)'*vtemp(:);
    end
end

%% check true solution
% we need to guess six variables - 
%   R = radius of hodograph
%   e = eccentricity
%   f0 = initial true anomaly
%   eul = 1x3 ZYX euler angles to rotate from the inertial frame to LVLH
%       where LVLH is the local-vertical local-horizontal frame with the
%       hodograph in the x-y plane, with its center offset toward the +y
%       axis.

R_true = sqrt(mu/a/(1-e^2));
e_true = e;
f0_true = f0;
eul_true = [0,0,0];

pulsars = P;
time = t_meas;

true_soln = [R_true,e_true,f0_true,eul_true];
obj = obj_6var(true_soln,mu,obsv,pulsars,time);

%% plot objective function value w.r.t. R+e perturbations

len = 40;
TEST = nan(len,len);
R_arr = linspace(R_true*0.6,R_true*1.4,len);
e_arr = linspace(e_true*0.6,e_true*1.4,len);
[RR,EE] = meshgrid(R_arr,e_arr);
for j = 1:len
    for k = 1:len
        input = [R_arr(j),e_arr(k),f0_true,eul_true];
        temp = obj_6var(input,mu,obsv,pulsars,time);
        TEST(j,k) = norm(temp(:));
    end
end
surf(RR,EE,TEST)

%% plot objective function value w.r.t. R+e perturbations

len = 40;
TEST = nan(len,len);
eul_1_arr = linspace(-pi/2,pi/2,len);
eul_2_arr = linspace(-pi/2,pi/2,len);
[ZZ,YY] = meshgrid(R_arr,e_arr);
for j = 1:len
    for k = 1:len
        input = [R_true,e_true,f0_true,eul_1_arr(j),eul_2_arr(k),0];
        temp = obj_6var(input,mu,obsv,pulsars,time);
        TEST(j,k) = norm(temp(:));
    end
end
surf(ZZ,YY,TEST)

%% guess solution
% Rmin = abs(min(obsv(:)))/2;
% R_guess = Rmin*2;
% e_guess = 0.1;
% f0_guess = 0;
% eul_guess = [pi/4,pi/4,pi/4];
% 
% guess = [R_guess,e_guess,f0_guess,eul_guess];
% 
% % fun = @(g) [ vecnorm(obj_6var(g,mu,obsv,pulsars,time),2,2) ; ...
% %              vecnorm(obj_6var(g,mu,obsv,pulsars,time),Inf,2) ];
% fun = @(g) obj_6var(g,mu,obsv,pulsars,time);
% 
% algo = 'trust-region-reflective';
% % algo = 'levenberg-marquardt';
% 
% options = optimoptions('lsqnonlin','Display','iter' ...
%                                   ,'MaxFunctionEvaluations',5000 ...
%                                   ,'StepTolerance', 1e-16 ...
%                                   ,'FunctionTolerance', 1e-16 ...
%                                   ,'OptimalityTolerance', 1e-16 ...
%                                   ,'MaxIterations',750 ...
%                                   ,'Algorithm',algo);
% g_opt = lsqnonlin(fun,guess,[Rmin,0,0,0,0,0],...
%                             [Inf,1,2*pi,2*pi,2*pi,2*pi],options);
% % g_opt = lsqnonlin(fun,guess,[],[],options);