function [r1] = IOD2V(v1,v2,u1,u2,mu)
%IOD2V Summary of this function goes here
%   Detailed explanation goes here

if abs(norm(v1)/norm(v2)-1) > 0.02
    r1 = -2 * mu * ...
        (norm(cross(v1,u1)) - norm(cross(v2,u2))) / ...
        ((v1*v1' - v2*v2') * norm(cross(v1,u1))) * ...
        u1;
else
    r1 = -mu / dot(v1,v1) * u1;
end

end

