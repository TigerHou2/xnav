close all
clear;clc

mu = 1;
a = 1;
e = 0.3;
i = pi/5;
omg = pi/4;
w = 0;

params = [a e i omg w 0];

up = [1;0;1]; % pulsar direction
up = up / norm(up);

% calculate true position and velocity values
f = deg2rad([330 30 90 150 210 270]);
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

[~,~,V] = svd(v');
k = V(:,end);
ux = cross(v(:,1),k);
ux = ux / norm(ux);
uy = cross(k,ux);
T = [ux';uy';k'];

% check that pulsar isn't normal to orbit plane
if dot(k,up) == 1
    error('Degenerate case: the pulsar is normal to the orbital plane.')
end

v

% ----- true hodograph ------
vt = T * v;
[Z,R] = fitcircle(vt(1:2,:));
Z = T' * [Z;0];
% find theta
uc = Z / norm(Z);
ue = cross([0;0;1],uc);
ue = ue / norm(ue);
k0 = cross(uc,ue);
theta = atan2(norm(cross(k0,k)),dot(k0,k));
if dot(cross(uc,k0),k)<0
    theta = -theta;
end
test = [uc,ue,k0];
%     % visualize for testing
%     figure
%     quiver3(0,0,0,uc(1),uc(2),uc(3));hold on
%     quiver3(0,0,0,ue(1),ue(2),ue(3))
%     quiver3(0,0,0,k0(1),k0(2),k0(3))
%     quiver3(0,0,0,k(1),k(2),k(3))
%     quiver3(0,0,0,up(1),up(2),up(3))
%     quiver3([0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1],'Color','Black')
%     scatter3(v(1,:),v(2,:),v(3,:))
%     hold off
%     xlabel('x');ylabel('y');zlabel('z');
%     legend('uc','ue','k0','k','pulsar')
%     axis equal
%     view([1 1 1])

% test if objective function returns 0
fun([Z',theta,R],up,v_obsv,mu,dt)

res = 20;
rGuess = R;
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
                dat(i,j,k) = norm(fun([xx(i),yy(j),zz(k),tt(m),rGuess],up,v_obsv,mu,dt))^0.25;
            end
        end
    end
end

[~,idx] = min(dat(:));
[i,j,k,m] = ind2sub(size(dat),idx);
f0 = [xx(i(1)),yy(j(1)),zz(k(1)),tt(m(1)),rGuess];

%%

f = @(x) fun(x,up,v_obsv,mu,dt);
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',420);
hodo = fsolve(f,f0,options)

fun(hodo,up,v_obsv,mu,dt)

%% function definitions
function t_diff = fun(C,up,v,mu,dt)

    x = C(1); % x-offset of hodograph
    y = C(2); % y-offset of hodograph
    z = C(3); % z-offset of hodograph
    theta = C(4); % hodograph rotation w.r.t line thru origin and center
                  % where theta = 0 is defined by the point where the
                  % orbit normal crosses the z-axis in the inertial frame
    R = C(5); % radius of hodograph
    
    e = sqrt(x^2 + y^2 + z^2) / R; % eccentricity
    
    % find k, the orbit normal vector
    uc = [x;y;z] / norm([x;y;z]);
    ue = cross([0;0;1],uc);
    ue = ue / norm(ue);
    k0 = cross(uc,ue);
    k  = rotVec(k0,uc,-theta);
    
    % calculate the rotation matrix T that moves the hodograph to 2D
    T = [uc';ue';k';];
    
    % calculate the pulsar direction in new hodograph frame
    up = T * up;
    
    % calculate the new hodograph center
    c_vect = T * [x;y;z];
    cx = c_vect(1);
    cy = c_vect(2);
    
    % calculate the angle between the pulsar vector and the orbital plane
    gamma = atan2(norm(cross(up,[0;0;1])),dot(up,[0;0;1]));
    gamma = pi/2 - gamma;
    gamma = mod(gamma,pi/2);
    
    uu = cross([0;0;1],up); % vector orthogonal to pulsar direction
    
    % solve intersection on the circle
    
    % --- convert the out of plane pulsar observation to an in-plane value
    v = v / cos(gamma);
    
    % --- calculate intercepts on each solution set
    
        % we only care about the in-plane components of up now
        % so make a new variable upp with only the x and y components
        % and make it a unit vector.
        % alternatively, simply divide up by cos^2 in the above step
    upp = [up(1:2);0] / norm([up(1:2);0]);
    P1 = upp * v';
    P2 = P1 + uu;
    
    I1 = nan(2,length(v));
    I2 = nan(2,length(v));
    for j = 1:length(v)
        dx = P2(1,j) - P1(1,j);
        dy = P2(2,j) - P1(2,j);
        dr = sqrt(dx^2+dy^2);
        D = (P1(1,j)-cx)*(P2(2,j)-cy) - (P2(1,j)-cx)*(P1(2,j)-cy);
        sgn = sign(dy); if sgn == 0, sgn = 1; end
        delta = R^2*dr^2-D^2;
        if delta < 0
            t_diff = dt;
            return
        end
        I1(1,j) = ( D*dy + sgn*dx*sqrt(delta) ) / dr^2 + cx;
        I1(2,j) = (-D*dx + sgn*dy*sqrt(delta) ) / dr^2 + cy;
        I2(1,j) = ( D*dy - sgn*dx*sqrt(delta) ) / dr^2 + cx;
        I2(2,j) = (-D*dx - sgn*dy*sqrt(delta) ) / dr^2 + cy;
    end
    
    % append zeros for finding true anomaly
    I1(end+1,:) = 0;
    I2(end+1,:) = 0;
    II = cat(3,I1,I2);
    
    I = I1;
    
    % iterate through solution pairs to find ordered true anomalies
    % for now we assume measurements do not cross periapsis (f=0)

    % --- solve for true anomaly at each point
    r1_vect = I1 - c_vect;
    r2_vect = I2 - c_vect;
    
    f1 = nan(1,size(r1_vect,2));
    for j = 1:size(r1_vect,2)
        f1(j) = atan2(norm(cross(c_vect,r1_vect(:,j))),dot(c_vect,r1_vect(:,j)));
        if dot(k,cross(c_vect,r1_vect(:,j)))<0
            f1(j) = 2*pi-f1(j);
        end
    end
    
    f2 = nan(1,size(r2_vect,2));
    for j = 1:size(r2_vect,2)
        f2(j) = atan2(norm(cross(c_vect,r2_vect(:,j))),dot(c_vect,r2_vect(:,j)));
        if dot(k,cross(c_vect,r2_vect(:,j)))<0
            f2(j) = 2*pi-f2(j);
        end
    end
    
    ff = [f1;f2];
    
    % --- find ordered true anomaly
    indices = dec2bin(0:1:2^length(v)-1)-'0' + 1;
    f_sorted = [];
    idx_sorted = [];
    for i = 1:size(indices,1)
        idx = sub2ind(size(ff),indices(i,:),1:length(v));
        if isMono(ff(idx))
            f_sorted(end+1,:) = ff(idx);
            idx_sorted(end+1,:) = indices(i,:);
        end
    end
    
    % --- get rid of rows with duplicate values
    f_tol = 1e-3;
    idx_sorted(sum(abs(diff(sort(f_sorted,2),1,2))>f_tol,2)+1~=length(v),:) = [];
    f_sorted(sum(abs(diff(sort(f_sorted,2),1,2))>f_tol,2)+1~=length(v),:) = [];
    
    % --- if no sorted true anomalies exist, pick a random one
    if isempty(f_sorted)
        f = f1;
        I = I1;
        
        a = mu/(I(:,1)'*I(:,1)) * (1 + e^2 + 2*e*cos(f(1))) / (1 - e^2);
    
        E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
        M = E - e*sin(E);
        t = sqrt(a^3/mu)*M;
        dt_guess = t(2:end) - t(1:end-1);
        dt_guess = mod(dt_guess,sqrt(a^3/mu));

        t_diff = dt_guess - dt;
        
        return
    end
    
    % --- otherwise, iterate through all sorted sets and find minimum error
    t_diff_vect = [];
    for j = 1:size(f_sorted,1)
        f = f_sorted(j,:);
        index = idx_sorted(j,:);
        for i = 1:length(v)
            I(:,i) = II(:,i,index(i));
        end
        
        a = mu/(I(:,1)'*I(:,1)) * (1 + e^2 + 2*e*cos(f(1))) / (1 - e^2);
    
        E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
        M = E - e*sin(E);
        t = sqrt(a^3/mu)*M;
        dt_guess = t(2:end) - t(1:end-1);

        t_diff_vect(end+1,:) = dt_guess - dt;
    end

    t_diff = t_diff_vect(vecnorm(t_diff_vect,2,2)==min(vecnorm(t_diff_vect,2,2)),:);
    
end

function [tf,idx] = isMono(S)
% ISMONO returns true if array S is loop-monotonic, and returns the 
%   index of the starting element in S going to the right
%
% Definition:
%   An array is loop-monotonic if it can be split into two segments that
%       can be concatenated to form a monotonic sequence.
%   Example:
%       [4,5,6,1,2,3] is loop-monotonic, because it can be split into
%       [4,5,6] and [1,2,3], and then concatenated to form [1,2,3,4,5,6].
%   The starting element is the first element in the concatenated array.
%       The index returned is the position of the starting element in the
%       original array S.

    S1 = [S(:);S(1)];
    S2 = [S(end);S(:)];
    
    diff = S1 - S2;
    
    % monotonically increasing
    if sum(double(diff<0)) == 1
        tf = true;
        idx = find(S==S(diff(1:end-1)<0));
        return
    elseif sum(double(diff<0)) == 0
        tf = true;
        idx = 1;
        return
    end
    
    % monotonically decreasing
    if sum(double(diff>0)) == 1
        tf = true;
        idx = find(S==S(diff(1:end-1)>0));
        return
    elseif sum(double(diff>0)) == 0
        tf = true;
        idx = 1;
        return
    end
    
    tf = false;

end
