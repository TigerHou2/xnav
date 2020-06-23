close all
clear;clc

% define orbit
mu = 1;
a = 1;
e = 0.4;
i = pi/9;
omg = pi/4;
w = pi/5;
params = [a e i omg w 0];

% define time offset
t_offset = 0.1378; % just has to be a number unrelated to anything else

% define pulsars
P = [0 0 1;... % pulsar 1
     1 1 1;... % pulsar 2
     3 1 2]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% calculate true position and velocity values
% --- at least 3 measurements for each pulsar
% --- if some pulsars have fewer measurements, fill to end with NaN
f = deg2rad([ 0  40  80  120;... % pulsar 1
             10  50  90  130;... % pulsar 2
             20  60  100 140]); % pulsar 3
E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);
t_true = sqrt(a^3/mu)*M;

t = t_true + t_offset;

% perform observations
r = nan(size(f,1),size(f,2),3);
v = nan(size(f,1),size(f,2),3);
obsv = nan(size(f));
pulsar = P;
for i = 1:size(f,1)
    for j = 1:size(f,2)
        params(6) = f(i,j);
        [r(i,j,:),v(i,j,:)] = Get_Orb_Vects(params,mu);
        vtemp = v(i,j,:);
        obsv(i,j) = P(:,i)'*vtemp(:);
    end
end

% calculate perfect inputs
period = 2*pi * sqrt(a^3/mu);
OPT = [e,period,t_offset];

% check if perfect input yields perfect output
[optDiff,V] = rrFun_sine(OPT,obsv,pulsar,mu,t)

%%
disp(['Max range-rate error: ' num2str(max(err))])
disp(['Max velocity error:   ' num2str(max(max(abs(v-vel))))])

%% global search

res = [12,12,12,12,12];
dat = nan(res);
lb = mean(abs(obsv));
ub = 3*mean(abs(obsv));
xx = linspace(-lb,lb,res(1));
yy = linspace(-lb,lb,res(2));
zz = linspace(-lb,lb,res(3));
tt = linspace(-pi,pi,res(4));
rr = linspace( lb,ub,res(5));
idx = 0;

for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
for m = 1:res(4)
for n = 1:res(5)
    idx = idx + 1;
    fin = [xx(i),yy(j),zz(k),tt(m),rr(n)];
    dat(i,j,k,m,n) = norm(rrFun_ta(fin,obsv,pulsar,mu,dt));
end
end
end
end
end

[~,idx] = min(dat(:));
[i,j,k,m,n] = ind2sub(size(dat),idx);
f0 = [xx(i(1)),yy(j(1)),zz(k(1)),tt(m(1)),rr(n(1))];

%% optimization

disp('Searching for solution...')

% f0 = hodo_true .* [1.02 1.02 1.02 1.02 1.02];

f = @(x) rrFun_ta(x,obsv,pulsar,mu,dt);
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',700 ...
                               ,'Algorithm','levenberg-marquardt' ...
                               ,'StepTolerance',1e-8 ...
                               ,'FunctionTolerance',1e-8);
hodo = fsolve(f,f0,options)

[v_err,vel] = rrFun_ta(hodo,obsv,pulsar,mu,dt,1);

v_diff = (v-vel)./vecnorm(v);
if max(abs(v_diff(:))) > 1e-3
    warning('Velocities do not converge!')
end


%% test solution space

f = @(x) rrFun_ta(x,obsv,pulsar,mu,dt);
options = optimoptions('fsolve','Display','none' ...
                               ,'MaxFunctionEvaluations',3000 ...
                               ,'MaxIterations',500 ...
                               ,'Algorithm','levenberg-marquardt' ...
                               ,'StepTolerance',1e-9);

guess_valid = [];
guess_inval = [];
for i = 1:300
    
    pert = (rand(1,5)-0.5);
    % ----- test section -----
% 	pert(end) = abs(pert(end));
    % ----- end test -----
    f0 = hodo_true + hodo_true .* pert;
    hodo = fsolve(f,f0,options);
    if max(abs(hodo-hodo_true)) < 1e-3
        guess_valid(end+1,:) = pert;
    else
        guess_inval(end+1,:) = pert;
    end
    
end

figure;
scatter(guess_inval(:,4),guess_inval(:,5),20);
hold on
scatter(guess_valid(:,4),guess_valid(:,5),20);
hold off
axis equal
xlabel('x')
ylabel('y')

disp(' ')
disp('Mean Disturbance')
disp(['Inval Case:' mat2str(mean(guess_inval),3)])
disp(['Valid Case:' mat2str(mean(guess_valid),3)])
disp(' ')
disp('Standard Deviation')
disp(['Inval Case:' mat2str(std(guess_inval),3)])
disp(['Valid Case:' mat2str(std(guess_valid),3)])
disp(' ')
disp('Correlation Coefficients:')
disp(corrcoef(guess_valid))