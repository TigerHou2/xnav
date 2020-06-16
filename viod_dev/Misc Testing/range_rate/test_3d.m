close all
clear;clc

mu = 1;
a = 1;
e = 0.3;
i = pi/9;
omg = pi/4;
w = 0;

params = [a e i omg w 0];

% define pulsars
P = [0 0 1;... % pulsar 1
     1 0 1;... % pulsar 2
     1 1 0]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% calculate true position and velocity values
f = deg2rad([330 10 20 40 60 70 80]);
E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);
t = sqrt(a^3/mu)*M;
dt = t(2:end) - t(1:end-1);
r = nan(3,length(f));
v = nan(3,length(f));
for j = 1:length(f)
    params(6) = f(j);
    [r(:,j),v(:,j)] = Get_Orb_Vects(params,mu);
end

v_obsv = nan(3,length(f));
% alternate between three pulsars
idx = 0;
for j = 1:length(f)
    idx = mod(idx,3) + 1;
    v_obsv(:,j) = v(:,j)' * P(:,idx) * P(:,idx);
end

[~,~,V] = svd(v');
k = V(:,end);
% if dot(k,[0;0;1]) < 0
%     k = -k;
% end
ux = cross(v(:,1),k);
ux = ux / norm(ux);
uy = cross(k,ux);
T = [ux';uy';k'];

v_obsv
v

% ----- true hodograph ------
vt = T * v;
[C,R] = fitcircle(vt(1:2,:));
Z = T' * [C;0];
% find theta
uc = Z / norm(Z);
ue = cross([0;0;1],uc);
ue = ue / norm(ue);
k0 = cross(uc,ue);
theta = atan2(norm(cross(k0,k)),dot(k0,k));
if dot(ue,k)>0
    theta = -theta;
end

    % visualize for testing
    figure
    % --- calculate points on hodograph
    angles = linspace(0,2*pi,100);
    cx = R * cos(angles) + C(1);
    cy = R * sin(angles) + C(2);
    circ = T' * [cx;cy;zeros(size(cx))];
    % --- plot hodograph
    plot3(circ(1,:),circ(2,:),circ(3,:),'DisplayName','Hodograph');
    hold on
    % --- plot inertial axes
    quiver3([0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1],'Color','Black')
    % --- plot circle center vector uc
    quiver3(0,0,0,uc(1),uc(2),uc(3),'DisplayName','uc');
    % --- plot eccentricity vector ue
    quiver3(0,0,0,ue(1),ue(2),ue(3),'DisplayName','ue')
    % --- plot hodograph normal at theta = 0
    quiver3(0,0,0,k0(1),k0(2),k0(3),'DisplayName','k0')
    % --- plot true hodograph normal
    quiver3(0,0,0,k(1),k(2),k(3),'DisplayName','k')
    % --- plot pulsar vector
    quiver3(0,0,0,P(1,1),P(2,1),P(3,1),'-.','LineWidth',1,'DisplayName','P1')
    quiver3(0,0,0,P(1,2),P(2,2),P(3,2),'--','LineWidth',1,'DisplayName','P2')
    quiver3(0,0,0,P(1,3),P(2,3),P(3,3),':' ,'LineWidth',1,'DisplayName','P3')
    % --- plot measurement points
    scatter3(v(1,:),v(2,:),v(3,:))
    hold off
    xlabel('x');ylabel('y');zlabel('z');legend
    axis equal
    view([1 1 1])

% test if objective function returns 0 given perfect input
hodo_true = [Z',theta,R];
% [err,val] = rrFun([Z',theta,R],v_obsv,mu,dt)
[err,val] = rrFun_ta([Z',theta,R],v_obsv,mu,dt)

%% global search

res = 15;
rGuess = 5*mean(vecnorm(v_obsv));
dat = nan(res,res,res);
xx = linspace(-0.5,0.5,res);
yy = linspace(-0.5,0.5,res);
zz = linspace(-0.5,0.5,res);
tt = linspace(-pi,pi,res);
idx = 0;

for i = 1:res
    for j = 1:res
        for k = 1:res
            for m = 1:res
                idx = idx + 1;
                dat(i,j,k) = norm(rrFun_ta([xx(i),yy(j),zz(k),tt(m),rGuess],v_obsv,mu,dt));
            end
        end
    end
end

[~,idx] = min(dat(:));
[i,j,k,m] = ind2sub(size(dat),idx);
f0 = [xx(i(1)),yy(j(1)),zz(k(1)),tt(m(1)),rGuess];

%% optimization

% f0 = hodo_true .* [1.02 1.02 1.02 1.02 1.02];

f = @(x) rrFun_ta(x,v_obsv,mu,dt);
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',500 ...
                               ,'Algorithm','levenberg-marquardt');
hodo = fsolve(f,f0,options)

[t_err,vel] = rrFun_ta(hodo,v_obsv,mu,dt,1)


%% test solution space

f = @(x) rrFun_ta(x,v_obsv,mu,dt);
options = optimoptions('fsolve','Display','none' ...
                               ,'MaxFunctionEvaluations',3000 ...
                               ,'MaxIterations',500 ...
                               ,'Algorithm','levenberg-marquardt' ...
                               ,'StepTolerance',1e-9);

guess_valid = [];
guess_inval = [];
for i = 1:300
    
    pert = (rand(1,5)-0.5)/10;
    % ----- test section -----
% 	pert(end) = abs(pert(end));
    % ----- end test -----
    f0 = hodo_true + hodo_true .* pert;
    [~,~,exitFlag] = fsolve(f,f0,options);
    if exitFlag > 0
        guess_valid(end+1,:) = pert;
    else
        guess_inval(end+1,:) = pert;
    end
    
end

figure;
scatter3(guess_inval(:,1),guess_inval(:,4),guess_inval(:,5),20);
hold on
scatter3(guess_valid(:,1),guess_valid(:,4),guess_valid(:,5),20);
hold off
axis equal

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