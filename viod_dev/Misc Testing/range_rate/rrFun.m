function [t_diff,V] = rrFun(C,obsv,mu,dt,testing)
%RRFUN Is the objective function for range-rate hodograph fitting.
%
% Author:
%   Tiger Hou
%
% Variables:
%   C     - 1x5 array, optimization variable
%           |- C(1) - x-offset of hodograph
%           |- C(2) - y-offset of hodograph
%           |- C(3) - z-offset of hodograph
%           |- C(4) - hodograph rotation along axis uc, where
%                     uc is defined by the origin and hodograph center
%           |- C(5) - radius of hodograph
%   obsv  - 3xN array of N range-rate measurements to arbitrary pulsars
%   mu    - gravitational parameter of central body
%   dt    - 1x(N-1) array of true time intervals between measurements
%
% Outputs:
%   t_diff - time different between guess intervals and true intervals
%   V      - velocities corresponding to hodograph defined by C

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
k  = rotVec(k0,uc,theta);
ue = cross(k,uc);

% calculate the rotation matrix T that moves the hodograph to 2D
% with uc along the x-axis
T = [uc';ue';k';];

% calculate the pulsar direction in new hodograph frame
obsv = T * obsv;

% calculate the new hodograph center
c_vect = T * [x;y;z];
cx = c_vect(1);
cy = c_vect(2);

V1 = nan(2,size(obsv,2));
V2 = nan(2,size(obsv,2));
for j = 1:size(obsv,2)
    
    % calculate perpendicular line for each pulsar
    P1 = obsv(:,j) * norm(obsv(:,j))^2 / norm(obsv(1:2,j))^2;
    P2 = P1 + cross([0;0;1],P1);
    
    % intercept line and circle
    dx = P2(1) - P1(1);
    dy = P2(2) - P1(2);
    dr = sqrt(dx^2+dy^2);
    D = (P1(1)-cx)*(P2(2)-cy) - (P2(1)-cx)*(P1(2)-cy);
    sgn = sign(dy); if sgn == 0, sgn = 1; end
    delta = R^2*dr^2-D^2;
    % no intercept, return objective function as a metric of how close the
    % line is to the circle
    if delta < 0
        t_diff = ones(size(dt)) * sqrt(-delta) * 20 + dt * 2;
        V = nan;
        return
    end
    V1(1,j) = ( D*dy + sgn*dx*sqrt(delta) ) / dr^2 + cx;
    V1(2,j) = (-D*dx + sgn*dy*sqrt(delta) ) / dr^2 + cy;
    V2(1,j) = ( D*dy - sgn*dx*sqrt(delta) ) / dr^2 + cx;
    V2(2,j) = (-D*dx - sgn*dy*sqrt(delta) ) / dr^2 + cy;
    
end

% append zeros for finding true anomaly
V1(end+1,:) = 0;
V2(end+1,:) = 0;
VV = cat(3,V1,V2);

V = V1;

% iterate through solution pairs to find ordered true anomalies
% for now we assume measurements do not cross periapsis (f=0)

% --- solve for true anomaly at each point
r1_vect = V1 - c_vect;
r2_vect = V2 - c_vect;

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
indices = dec2bin(0:1:2^size(obsv,2)-1)-'0' + 1;
f_sorted = [];
idx_sorted = [];
for i = 1:size(indices,1)
    idx = sub2ind(size(ff),indices(i,:),1:size(obsv,2));
    if isMono(ff(idx))
        f_sorted(end+1,:) = ff(idx);
        idx_sorted(end+1,:) = indices(i,:);
    end
end

% --- get rid of rows with duplicate values
f_tol = 1e-3;
idx_sorted(sum(abs(diff(sort(f_sorted,2),1,2))>f_tol,2)+1~=size(obsv,2),:) = [];
f_sorted(sum(abs(diff(sort(f_sorted,2),1,2))>f_tol,2)+1~=size(obsv,2),:) = [];

% --- if no sorted true anomalies exist, pick a random one
if isempty(f_sorted)
    f = f1;
    V = V1;

    a = mu/(V(:,1)'*V(:,1)) * (1 + e^2 + 2*e*cos(f(1))) / (1 - e^2);

    E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M = E - e*sin(E);
    t = sqrt(a^3/mu)*M;
    dt_guess = t(2:end) - t(1:end-1);
    dt_guess = mod(dt_guess,sqrt(a^3/mu));

    t_diff = dt_guess - dt;

% --- otherwise, iterate through all sorted sets and find minimum error
else
    t_diff_mat = [];
    Is = nan(size(V,1),size(V,2),size(f_sorted,1));
    for j = 1:size(f_sorted,1)
        f = f_sorted(j,:);
        index = idx_sorted(j,:);
        for i = 1:size(obsv,2)
            V(:,i) = VV(:,i,index(i));
        end

        Is(:,:,j) = V;

        a = mu/(V(:,1)'*V(:,1)) * (1 + e^2 + 2*e*cos(f(1))) / (1 - e^2);

        E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
        M = E - e*sin(E);
        t = sqrt(a^3/mu)*M;
        dt_guess = t(2:end) - t(1:end-1);

        t_diff_mat(end+1,:) = dt_guess - dt;
    end

    t_diff_vect = vecnorm(t_diff_mat,2,2);
    idx = find(t_diff_vect==min(t_diff_vect));
    idx = idx(1);
    t_diff = t_diff_mat(idx,:);
    V = Is(:,:,idx);

end

%------------ TESTING --------------
if exist('testing','var')
    % convert everything to 2D
    uc = T * uc;
    ue = T * ue;
    k = T * k;
    k0 = T * k0;
    obsv = obsv ./ vecnorm(obsv,2,1);

    figure
    % --- calculate points on hodograph
    angles = linspace(0,2*pi,100);
    cX = R * cos(angles) + c_vect(1);
    cY = R * sin(angles) + c_vect(2);
    circ = [cX;cY;zeros(size(cX))];
    % --- plot hodograph
    plot3(circ(1,:),circ(2,:),circ(3,:),'DisplayName','Hodograph');
    hold on
    % --- plot inertial axes
    quiver3([0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1],'Color','Black','DisplayName','Axes')
    % --- plot circle center vector uc
    quiver3(0,0,0,uc(1),uc(2),uc(3),'DisplayName','uc');
    % --- plot eccentricity vector ue
    quiver3(0,0,0,ue(1),ue(2),ue(3),'DisplayName','ue')
    % --- plot hodograph normal at theta = 0
    quiver3(0,0,0,k0(1),k0(2),k0(3),'DisplayName','k0')
    % --- plot true hodograph normal
    quiver3(0,0,0,k(1),k(2),k(3),'DisplayName','k')
    % --- plot pulsar vector
    quiver3(0,0,0,obsv(1,1),obsv(2,1),obsv(3,1),'-.','LineWidth',1,'DisplayName','P1')
    quiver3(0,0,0,obsv(1,2),obsv(2,2),obsv(3,2),'--','LineWidth',1,'DisplayName','P2')
    quiver3(0,0,0,obsv(1,3),obsv(2,3),obsv(3,3),':' ,'LineWidth',1,'DisplayName','P3')
    % --- plot guess velocities
    scatter3(V(1,:),V(2,:),V(3,:),'DisplayName','Guess Obsv Positions')
    hold off
    xlabel('x');ylabel('y');zlabel('z');legend
    axis equal
    view([1 1 1])
end
%------------ END TEST -------------

% final processing of outputs
V = T' * V; % conversion back to 3D

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