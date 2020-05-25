function [r1,r2,r3,k] = IOD3V_V2(v1,v2,v3,mu,varargin)
%IOD3V_V2 Solves the three-velocity initial orbit determination problem.
%
%   The spacecraft's orbit is calculated from three velocity vectors
%   collected at different points in the orbit. The orbit is specified by
%   either providing retrograde/prograde information, or by providing the
%   order of measurements.
%
%   LIMITATIONS: assumes measurements are taken within one revolution.
%
%   v:  velocity vectors [km/s]
%   mu: gravitational parameter of central body
%   Name-Value Pairs:
%       'ordered',  [boolean]: if true, assume v1,v2,v3 measured in order
%       'omega'  ,  [1x3]:     angular velocity vector of central body
%       'prograde', [boolean]: true -> prograde, false -> retrograde
%       'display',  [boolean]: if true, prints assumptions on command line

p = inputParser;
addParameter(p, 'ordered' , false  , @islogical);
addParameter(p, 'omega'   , [0 0 0], @(x) isnumeric(x) && all(size(x)==[1,3]));
addParameter(p, 'prograde', false  , @islogical);
addParameter(p, 'display' , false  , @islogical);
parse(p,varargin{:});

N = [v1; v2; v3];
[~,~,V] = svd(N);
k = V(:,3);

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
    [r1,r2,r3] = decomp(v1,v2,v3,k,mu);
elseif p.Results.ordered == true
    if p.Results.display
        disp('Assuming velocity vectors are ordered')
        disp(newline)
    end
    % the velocty vectors are ordered
    [r1,r2,r3] = decomp(v1,v2,v3,k,mu);
    
    delta_E1 = get_ecc(r2,v2,mu) - get_ecc(r1,v1,mu);
    if delta_E1 < 0
        delta_E1 = delta_E1 + 2*pi;
    end
    
    delta_E2 = get_ecc(r3,v3,mu) - get_ecc(r2,v2,mu);
    if delta_E2 < 0
        delta_E2 = delta_E2 + 2*pi;
    end
    
    if (delta_E1 + delta_E2) > 2*pi
        k = -k;
    end
    [r1,r2,r3] = decomp(v1,v2,v3,k,mu);
else
    [r1,r2,r3] = decomp(v1,v2,v3,k,mu);
end

end


function [r1,r2,r3] = decomp(v1,v2,v3,k,mu)
    u1 = v1' / norm(v1);
    w1 = cross(u1,k) / norm(cross(u1,k));
    u2 = v2' / norm(v2);
    w2 = cross(u2,k) / norm(cross(u2,k));
    u3 = v3' / norm(v3);
    w3 = cross(u3,k) / norm(cross(u3,k));

    z1 = w1 / norm(v1);
    z2 = w2 / norm(v2);
    z3 = w3 / norm(v3);

    A = [z1         u1         -z2         -u2         zeros(3,1)  zeros(3,1); ...
         z1         u1         zeros(3,1)  zeros(3,1)  -z3         -u3; ...
         zeros(3,1) zeros(3,1) z2          u2          -z3         -u3; ...
         1          0          -1          0           0           0; ...
         1          0          0           0           -1          0; ...
         0          0          1           0           -1          0];
    B = [norm(v1)*w1 - norm(v2)*w2; ...
         norm(v1)*w1 - norm(v3)*w3; ...
         norm(v2)*w2 - norm(v3)*w3; ...
         norm(v1)^2/2 - norm(v2)^2/2; ...
         norm(v1)^2/2 - norm(v3)^2/2; ...
         norm(v2)^2/2 - norm(v3)^2/2];
    g = A\B;

    a1 = g(1);
    b1 = g(2);
    a2 = g(3);
    b2 = g(4);
    a3 = g(5);
    b3 = g(6);

    r1_norm = mu / a1;
    r2_norm = mu / a2;
    r3_norm = mu / a3;

    h1 = sqrt( r1_norm^2*dot(v1,v1) / (1 + dot(v1,v1)*(b1/a1)^2) );
    h2 = sqrt( r2_norm^2*dot(v2,v2) / (1 + dot(v2,v2)*(b2/a2)^2) );
    h3 = sqrt( r3_norm^2*dot(v3,v3) / (1 + dot(v3,v3)*(b3/a3)^2) );

    h = mean([h1,h2,h3]);

    rw1 = h*z1;
    ru1 = b1*h*r1_norm/mu*u1;

    rw2 = h*z2;
    ru2 = b2*h*r2_norm/mu*u2;

    rw3 = h*z3;
    ru3 = b3*h*r3_norm/mu*u3;

    r1 = rw1 + ru1;
    r1 = r1';

    r2 = rw2 + ru2;
    r2 = r2';

    r3 = rw3 + ru3;
    r3 = r3';
end


function E = get_ecc(r,v,mu)
    e_vec = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu;
    e = norm(e_vec);
    f = acos(dot(e_vec,r)/e/norm(r));
    if dot(r,v) < 0
        f = 2*pi - f;
    end
    if e < 1
        E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2));
    else
        E = 2*atanh(sqrt((1-e)/(1+e))*tan(f/2));
    end
end