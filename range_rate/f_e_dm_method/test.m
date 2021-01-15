close all
clear;clc

% profile on

addpath('../fcns_circ')
addpath('../fcns_orb')

% define orbit
mu = 1;
a = 1e5;
i = deg2rad(13);
omg = deg2rad(25);
w = deg2rad(90);

% define true solution
e = 0.5;
f0 = deg2rad(90);
dM = deg2rad(30);

OPT = [f0,e,dM];

params = [a e i omg w 0];

% define pulsars
P = [0 -1 1;... % pulsar 1
     1 0 1;... % pulsar 2
     -1 0 3]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% number of measurements
numObsv = 10;

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

% visualize range rate
figure(1)
plot(f(1,:),obsv(1,:),'ro')
hold on
plot(f(2,:),obsv(2,:),'go')
plot(f(3,:),obsv(3,:),'bo')
hold off

% visualize pulsars
figure(100)
quiver3(0,0,0,P(1,1),P(2,1),P(3,1),'DisplayName','Pulsar 1')
hold on
quiver3(0,0,0,P(1,2),P(2,2),P(3,2),'DisplayName','Pulsar 2')
quiver3(0,0,0,P(1,3),P(2,3),P(3,3),'DisplayName','Pulsar 3')
hold off
axis equal
pbaspect([1 1 1])
legend('Location','Best')

% check if perfect input yields perfect output
ref_error = guess(OPT,obsv,pulsar,t_meas,2);
disp(['Error of true soln.: ' num2str(max(ref_error))])

%% global search

res = [30,30,30];
range_f0 = linspace(0,2*pi,res(1)+1);
range_f0(end) = [];
range_e = linspace(0,1,res(2)+1);
range_e(end) = [];
range_dM = linspace(0,2*pi,res(3)+1);
range_dM(1) = [];

fun1 = search_adv(range_f0,range_e,range_dM,res,obsv,pulsar,t_meas,OPT,'plot');

%% find the sine wave

disp('Using fsolve...')

fun = @(x) guess(x,obsv,pulsar,t_meas,0);
options = optimoptions('lsqnonlin','Display','iter' ...
                                  ,'MaxFunctionEvaluations',12000 ...
                                  ,'StepTolerance', 1e-16 ...
                                  ,'FunctionTolerance', 1e-16 ...
                                  ,'OptimalityTolerance', 1e-16 ...
                                  ,'MaxIterations',3000 ...
                                  ,'Algorithm','trust-region-reflective');
% g_opt = fsolve(fun,fun1,options);
g_opt = lsqnonlin(fun,fun1,[0,0,0],[2*pi,1,2*pi],options);

disp('Solve complete!')

soln_opt = g_opt;
soln_opt(1) = mod(g_opt(1),2*pi);
soln_opt(3) = mod(g_opt(3),2*pi);

error = guess(soln_opt,obsv,pulsar,t_meas,1);

disp(['Max error term: ' num2str(max(error))])
disp(['Soln. difference: ' num2str(soln_opt-OPT)])

fun1 = [0 0 0];

%% refined search

res = [20,20,20];
res = res + mod(res,2); % preserves the current best point

delta_min = [2*pi,1,2*pi] * 0.01;

fun2 = soln_opt;
delta = 3*abs(fun2-fun1)+delta_min;
delta(1) = min([fun2(1),2*pi-fun2(1),delta(1)]);
delta(2) = min([fun2(2),   1-fun2(2),delta(2)]);
delta(3) = min([fun2(3),2*pi-fun2(3),delta(3)]);
fun1 = fun2;

lb = soln_opt(1)-delta(1);
rb = soln_opt(1)+delta(1);
range_f0 = linspace(lb,rb,res(1)+1);
range_f0(end) = [];

lb = soln_opt(2)-delta(2);
rb = soln_opt(2)+delta(2);
range_e = linspace(lb,rb,res(2)+1);
range_e(end) = [];

lb = soln_opt(3)-delta(3);
rb = soln_opt(3)+delta(3);
range_dM = linspace(lb,rb,res(3)+1);
range_dM(1) = [];

fun2 = search_adv(range_f0,range_e,range_dM,res,obsv,pulsar,t_meas,OPT,'plot');

fun = @(x) guess(x,obsv,pulsar,t_meas,0);
options = optimoptions('lsqnonlin','Display','iter' ...
                                  ,'MaxFunctionEvaluations',12000 ...
                                  ,'StepTolerance', 1e-18 ...
                                  ,'FunctionTolerance', 1e-18 ...
                                  ,'OptimalityTolerance', 1e-18 ...
                                  ,'MaxIterations',3000 ...
                                  ,'Algorithm','trust-region-reflective');
% g_opt = fsolve(fun,fun1,options);
g_opt = lsqnonlin(fun,fun2,[0,0,0],[2*pi,1,2*pi],options);

soln_opt = g_opt;
soln_opt(1) = mod(g_opt(1),2*pi);
soln_opt(3) = mod(g_opt(3),2*pi);

error = guess(soln_opt,obsv,pulsar,t_meas,1);

disp(['Max error term: ' num2str(max(error))])
disp(['Soln. difference: ' num2str(soln_opt-OPT)])

%% output initial orbit

orb = get_io(soln_opt,t_meas,mu);

disp([newline '----------'])
disp(['Guess SMA = ' num2str(orb.a) newline...
      ' True SMA = ' num2str(a) newline])
disp(['Guess ECC = ' num2str(orb.e) newline...
      ' True ECC = ' num2str(e) newline])
disp(['Guess TA = ' num2str(orb.f0) newline...
      ' True TA = ' num2str(f0) newline])
  
% profile viewer

%% function definitions
function Evect = kepler(Mvect,e)

Evect = nan(size(Mvect));

for i = 1:length(Mvect(:))
    
M = Mvect(i);

if M == 0
    Evect(i) = 0;
    continue
end

d = 1;
n = 5;

if e < 1
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(E-e*sin(E)-M) / ...
            (1-e*cos(E) + sign(1-e*cos(E)) ...
                        * sqrt(abs((n-1)^2*(1-e*cos(E))^2 ...
                                   -n*(n-1)*(E-e*sin(E)-M)*(e*sin(E)))));
        d = abs(E-e*sin(E)-M);
    end
else
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(e*sinh(E)-E-M) / ...
            (e*cosh(E)-1 + sign(e*cosh(E)-1) ...
                         * sqrt(abs((n-1)^2*(e*cosh(E)-1)^2 ...
                                    -n*(n-1)*(e*sinh(E)-E-M)*(e*sinh(E)))));
        d = abs(e*sinh(E)-E-M);
    end
end

Evect(i) = E;

end

end