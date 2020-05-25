function [r,k] = IOD3V_Ext(N,mu,varargin)
%IOD3V_EXT Solves the three-velocity initial orbit determination problem.
%
%   The spacecraft's orbit is calculated from three velocity vectors
%   collected at different points in the orbit. The orbit is specified by
%   either providing retrograde/prograde information, or by providing the
%   order of measurements.
%
%   LIMITATIONS: assumes measurements are taken within one revolution.
%
%   Changes since IOD_V2:
%       - fixed hyperbolic orbit case eccentric anomaly equation
%       - added ability to process any number of velocity inputs
%
%   N:  velocity vectors [km/s]
%   mu: gravitational parameter of central body
%   Name-Value Pairs:
%       'ordered',  [boolean]: if true, assume v1,v2,v3 measured in order
%       'omega'  ,  [1x3]:     angular velocity vector of central body
%       'prograde', [boolean]: true -> prograde, false -> retrograde
%       'display',  [boolean]: if true, prints assumptions on command line
%
%   PAPER: Initial Orbit Determination from Three Velocity Vectors

p = inputParser;
addParameter(p,'ordered' ,false  ,@islogical);
addParameter(p,'omega'   ,[0 0 0],@(x) isnumeric(x)&&all(size(x)==[1,3]));
addParameter(p,'prograde',false  ,@islogical);
addParameter(p,'display' ,false  ,@islogical);
parse(p,varargin{:});

[~,~,V] = svd(N);
k = V(:,end);

if not(all(p.Results.omega == [0 0 0]))
    if p.Results.display
        disp('Assuming prograde / retrograde orbit')
        disp(newline)
    end
    % orbit is known to be prograde or retrograde
    if p.Results.prograde
        k = sign(p.Results.omega*k)*k;
    else
        k = -sign(p.Results.omega*k)*k;
    end
    r = decomp(N,k,mu,false);
elseif p.Results.ordered == true
    if p.Results.display
        disp('Assuming velocity vectors are ordered')
        disp(newline)
    end
    % the velocty vectors are ordered
    [r,k] = decomp(N,k,mu,true);
else
    r = decomp(N,k,mu,false);
end

end


function [r,k] = decomp(N,k,mu,order)
    reps = size(N,1);
    u = zeros(size(N'));
    w = zeros(size(N'));
    z = zeros(size(N'));
    for i = 1:reps
        u(:,i) = N(i,:)' / norm(N(i,:));
        w(:,i) = cross(u(:,i),k) / norm(cross(u(:,i),k));
        z(:,i) = w(:,i) / norm(N(i,:));
    end
    
    C = combnk(1:reps,2);           % permutations of two velocity vectors
    A = zeros(size(C,1)*4,reps*2);  % LHS matrix of eq. 36, see line 23
    B = zeros(size(C,1)*4,1);       % RHS matrix of eq. 36, see line 23
    for i = 1:size(C,1)
        % insert z and u values in top 3/4 of LHS matrix
        A((i-1)*3+1:i*3,2*C(i,1)-1) =  z(:,C(i,1));
        A((i-1)*3+1:i*3,2*C(i,1))   =  u(:,C(i,1));
        A((i-1)*3+1:i*3,2*C(i,2)-1) = -z(:,C(i,2));
        A((i-1)*3+1:i*3,2*C(i,2))   = -u(:,C(i,2));
        % insert 1 and -1 in bottom 1/4 of LHS matrix
        A((i-1)*3+size(C,1)*3+1,2*C(i,1)-1) =  1;
        A((i-1)*3+size(C,1)*3+1,2*C(i,2)-1) = -1;
        % insert vectors in top 3/4 of RHS matrix
        B((i-1)*3+1:i*3) = norm(N(C(i,1),:))*w(:,C(i,1)) - ...
                           norm(N(C(i,2),:))*w(:,C(i,2));
        % insert vectors in bottom 1/4 of RHS matrix
        B((i-1)*3+size(C,1)*3+1) =  norm(N(C(i,1),:))^2/2 - ...
                                    norm(N(C(i,2),:))^2/2;
    end
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
            [r,k] = decomp(N,k,mu,false);
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