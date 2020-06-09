close all
clear;clc

mu = 1;
a = 1;
e = 0.3;
i = 0;
omg = 0;
w = 0;

params = [a e i omg w 0];

up = [1;0;0]; % pulsar direction
up = up / norm(up);
k = [0;0;1]; % orbit normal

% calculate true values
f = deg2rad([30 70 110 150]);
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

v_obsv = v'*up;


% get the actual hodograph
[Z,R] = fitcircle(v(1:2,:));


f0 = [Z' R/100];
fun(f0,up,v_obsv,mu,dt,k)

f0 = [0.1 0.1 0.02];

f = @(x) fun(x,up,v_obsv,mu,dt,k);
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',300);
hodo = fsolve(f,f0,options)

fun(hodo,up,v_obsv,mu,dt,k)

%% function definitions
function t_diff = fun(C,up,v,mu,dt,k)

    % we will assume 2D for now

    p = C(1); % x-offset of hodograph
    q = C(2); % y-offset of hodograph
    R = C(3)*100; % radius of hodograph
    
    e = sqrt(p^2 + q^2) / R; % eccentricity
    uu = cross(k,up); % vector orthogonal to pulsar direction
    
    % solve intersection on the circle
    
    % --- calculate point on each solution set
    P1 = up * v';
    P2 = P1 + uu;
    
    I1 = nan(2,length(v));
    I2 = nan(2,length(v));
    for j = 1:length(v)
        dx = P2(1,j) - P1(1,j);
        dy = P2(2,j) - P1(2,j);
        dr = sqrt(dx^2+dy^2);
        D = (P1(1,j)-p)*(P2(2,j)-q) - (P2(1,j)-p)*(P1(2,j)-q);
        sgn = sign(dy);
        if sgn == 0
            sgn = 1;
        end
        delta = R^2*dr^2-D^2;
        if delta < 0
            t_diff = ones(size(dt))*10000;
            return
        end
        I1(1,j) = ( D*dy + sgn*dx*sqrt(delta) ) / dr^2 + p;
        I1(2,j) = (-D*dx + sgn*dy*sqrt(delta) ) / dr^2 + q;
        I2(1,j) = ( D*dy - sgn*dx*sqrt(delta) ) / dr^2 + p;
        I2(2,j) = (-D*dx - sgn*dy*sqrt(delta) ) / dr^2 + q;
    end
    
    % append zeros for finding true anomaly
    I1(end+1,:) = 0;
    I2(end+1,:) = 0;
    
    I = I1;
    
    for j = 2:size(I,2)
        if all(ismembertol(I(:,j),I(:,1:j-1),1e-3))
            I(:,j) = I2(:,j);
        end
    end
    
    % solve for true anomaly at each point
    c_vect = [p;q;0]; % vector for center of hodograph
    r_vect = I - c_vect;
    
    f = nan(1,size(r_vect,2));
    for j = 1:size(r_vect,2)
        f(j) = atan2(norm(cross(c_vect,r_vect(:,j))),dot(c_vect,r_vect(:,j)));
    end
    
    a = mu/(I(:,1)'*I(:,1)) * (1 + e^2 + 2*e*cos(f(1))) / (1 - e^2);
    
    E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M = E - e*sin(E);
    t = sqrt(a^3/mu)*M;
    dt_guess = t(2:end) - t(1:end-1);
    
    t_diff = dt_guess - dt;
    
end
