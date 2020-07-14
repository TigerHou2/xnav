function rv_dot = gravity(rv,varargin)
%GRAVITY Summary of this function goes here
%   Detailed explanation goes here

m = varargin{1};
G = varargin{2};

N = size(rv,2);
rv_dot = zeros(size(rv));

for i = 1:N
    % velocity
    rv_dot(1:3,i) = rv(4:6,i);
    % acceleration
    for j = 1:N
        if i ~= j
        rv_dot(4:6,i) = rv_dot(4:6,i) + ...
                        G*m(j)*(rv(1:3,j)-rv(1:3,i)) / norm(rv(1:3,j)-rv(1:3,i))^3;
        end
    end
end

end

