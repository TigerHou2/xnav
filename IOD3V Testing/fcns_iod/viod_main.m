%% viod_main.m
% Author: Tiger Hou
%
%% function definition

function [r,k] = viod_main(V,A,B,C,mu)
%VIOD_MAIN Solves the three-velocity initial orbit determination problem.
%
% The spacecraft's orbit is calculated from velocity vectors
%   collected at different points in the orbit. The velocity measurements
%   must be provided in order to properly determine orbit orientation.
%
% Arguments:
%   N:      [km/s]      N-by-3 matrix containing N velocity vectors
%   A:      [nd]        [2*N*(N-1)]-by-[2*N] sparse matrix of linear system
%   B:      [nd]        [2*N*(N-1)]-by-1 matrix of linear system
%   C:      [nd]        [N*(N-1)/2]-by-2 matrix of combinations of
%                           two numbers selected from 1:N, non-repeating
%   mu:     [km^3/s^2]  gravitational parameter of central body
%
% Outputs:
%   r:      [km]        N-by-3 matrix containing N position vectors
%                           corresponding to the N velocity vectors in N
%   k:      [nd]        3-by-1 orbit plane normal vector
%
% Sources:
%   - [1] Initial Orbit Determination from Three Velocity Vectors

[~,~,VV] = svd(V);
k = VV(:,end);

[r,k] = decomp(V,A,B,C,k,mu,true);

end


function [r,k] = decomp(V,A,B,C,k,mu,order)
    reps = size(V,1);
    u = zeros(size(V')); % [3xN] matrix of u values defined in [1] eq.11
    w = zeros(size(V')); % [3xN] matrix of w values defined in [1] eq.12
    z = zeros(size(V')); % [3xN] matrix of z values defined in [1] eq.21
    for i = 1:reps
        u(:,i) = V(i,:)' / norm(V(i,:));
        w(:,i) = cross(u(:,i),k) / norm(cross(u(:,i),k));
        z(:,i) = w(:,i) / norm(V(i,:));
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
        B((i-1)*3+1:i*3) = norm(V(C(i,1),:))*w(:,C(i,1)) - ...
                           norm(V(C(i,2),:))*w(:,C(i,2));
        % insert vectors in bottom 1/4 of RHS matrix
        B((i-1)*3+nperms*3+1) =  norm(V(C(i,1),:))^2/2 - ...
                                    norm(V(C(i,2),:))^2/2;
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
        h(i) = sqrt( r_norm(i)^2*dot(V(i,:),V(i,:)) / ...
                (1 + dot(V(i,:),V(i,:))*(b(i)/a(i))^2) );
    end
    h = mean(h);

    rw = zeros(size(V'));
    ru = zeros(size(V'));
    r  = zeros(size(V'));
    
    for i = 1:reps
        rw(:,i) = h*z(:,i);
        ru(:,i) = b(i)*h*r_norm(i)/mu*u(:,i);
        r(:,i)  = rw(:,i) + ru(:,i);
    end
    r = r';
    
    if order
        % single recursion
        % check if calculated position corresponds to velocity order
        % this works because we know that the velocity vector rotates
        %   by 180 deg exactly when the position vector rotates by 180 deg
        % so by comparing the sign of the cross product, we can determine
        %   whether the orientation of the orbit is correct.
        % furthermore, we only need to check two vectors instead of all
        
        norm_v = cross(V(1,:),V(2,:));
        norm_r = cross(r(1,:),r(2,:));
        if dot(norm_v,norm_r) < 0
            k = -k;
            [r,k] = decomp(V,A,B,C,k,mu,false);
        end
    end
        
end