function [u0] = rk4(fcn,u0,T,dt,varargin)
%RK4 Summary of this function goes here
%   Detailed explanation goes here

for jj = 1 : T/dt
    
    % ============= RK4 =============
    k1 = dt * fcn( u0,varargin{:});
    k2 = dt * fcn((u0 + k1/2),varargin{:});
    k3 = dt * fcn((u0 + k2/2),varargin{:});
    k4 = dt * fcn((u0 + k3)  ,varargin{:});
    u0 = u0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
end

end

