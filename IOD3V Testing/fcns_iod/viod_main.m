%% viod_main.m
%
% Created:  03/29/2020
% Modified: 03/29/2020
% Author:   Tiger Hou
%
% 

%% function definition

function [r,k] = viod_main(N,A,B,C,mu)
%VIOD_MAIN Solves the three-velocity initial orbit determination problem.
%
%
% The spacecraft's orbit is calculated from velocity vectors
%   collected at different points in the orbit. The velocity measurements
%   must be provided in order to properly determine orbit orientation.
%
%
% Arguments:
%   N:      [km/s]      N-by-3 matrix containing N velocity vectors
%   A:      [unitless]  [2*N*(N-1)]-by-[2*N] sparse matrix of linear system
%   B:      [unitless]  [2*N*(N-1)]-by-1 matrix of linear system
%   C:      [unitless]  [N*(N-1)/2]-by-2 matrix of combinations of
%                           two numbers selected from 1:N, non-repeating
%   mu:     [km^3/s^2]  gravitational parameter of central body
%
% Sources:
%   - [1] Initial Orbit Determination from Three Velocity Vectors

[~,~,V] = svd(N);
k = V(:,end);

if not(all(p.Results.omega == [0 0 0]))
    % orbit is known to be prograde or retrograde
    % because omega (orbiting body angular velocity) is provided
    if p.Results.prograde
        k = sign(p.Results.omega*k)*k;
    else
        k = -sign(p.Results.omega*k)*k;
    end
    r = decomp(N,A,B,C,k,mu,false);
elseif p.Results.ordered == true
    % the velocty vectors are ordered
    [r,k] = decomp(N,A,B,C,k,mu,true);
else
    r = decomp(N,A,B,C,k,mu,false);
end

end


function [r,k] = decomp(N,A,B,C,k,mu,order)
    reps = size(N,1);
    u = zeros(size(N')); % [3xN] matrix of u values defined in [1] eq.11
    w = zeros(size(N')); % [3xN] matrix of w values defined in [1] eq.12
    z = zeros(size(N')); % [3xN] matrix of z values defined in [1] eq.21
    for i = 1:reps
        u(:,i) = N(i,:)' / norm(N(i,:));
        w(:,i) = cross(u(:,i),k) / norm(cross(u(:,i),k));
        z(:,i) = w(:,i) / norm(N(i,:));
    end
    
    nperms = size(C,1);
    for i = 1:nperms
        % insert z and u values in top 3/4 of LHS matrix
        A((i-1)*3+1:i*3,2*C(i,1)-1) =  z(:,C(i,1));
        A((i-1)*3+1:i*3,2*C(i,1))   =  u(:,C(i,1));
        A((i-1)*3+1:i*3,2*C(i,2)-1) = -z(:,C(i,2));
        A((i-1)*3+1:i*3,2*C(i,2))   = -u(:,C(i,2));
        % insert 1 and -1 in bottom 1/4 of LHS matrix
        A((i-1)*3+nperms*3+1,2*C(i,1)-1) =  1;
        A((i-1)*3+nperms*3+1,2*C(i,2)-1) = -1;
        % insert vectors in top 3/4 of RHS matrix
        B((i-1)*3+1:i*3) = norm(N(C(i,1),:))*w(:,C(i,1)) - ...
                           norm(N(C(i,2),:))*w(:,C(i,2));
        % insert vectors in bottom 1/4 of RHS matrix
        B((i-1)*3+nperms*3+1) =  norm(N(C(i,1),:))^2/2 - ...
                                    norm(N(C(i,2),:))^2/2;
    end
    % solve the linear system
    g = A\B;
    
    a = zeros(1,reps);
    b = zeros(1,reps);
    
    a(:) = g(1:2:end);
    b(:) = g(2:2:end);

    r_norm = mu ./ a;
    
    h = zeros(1,reps);
    for i = 1:reps
        h(i) = sqrt( r_norm(i)^2*dot(N(i,:),N(i,:)) / ...
                (1 + dot(N(i,:),N(i,:))*(b(i)/a(i))^2) );
    end
    h = mean(h);

    rw = zeros(size(N'));
    ru = zeros(size(N'));
    r  = zeros(size(N'));
    
    for i = 1:reps
        rw(:,i) = h*z(:,i);
        ru(:,i) = b(i)*h*r_norm(i)/mu*u(:,i);
        r(:,i)  = rw(:,i) + ru(:,i);
    end
    r = r';
    
    if order
        % single recursion
        delta_E = zeros(1,reps-1);
        for i = 1:reps-1
            delta_E(i) = get_ecc(r(i+1,:),N(i+1,:),mu) - ...
                         get_ecc(r(i,:),N(i,:),mu);
        end
        if any(delta_E<0)
            k = -k;
            [r,k] = decomp(N,A,B,C,k,mu,false);
        end
    end
        
end


function E = get_ecc(r,v,mu)
    a = norm(r)/(2 - norm(r)*dot(v,v)/mu);
    e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu;
    e = norm(e_vec);
    f = acos((a*(1-e^2)/norm(r)-1)/e);
    if dot(r,v) < 0
        f = 2*pi - f;
    end
    if e < 1
        E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2));
    else
        E = 2*atanh(sqrt((e-1)/(e+1))*tan(f/2));
    end
end