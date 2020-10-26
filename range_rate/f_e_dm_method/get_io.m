function orb = get_io(soln_opt,time,mu)
%GET_IO Summary of this function goes here
%   Detailed explanation goes here

orb.f0 = soln_opt(1);
orb.e  = soln_opt(2);
orb.period = 2*pi/soln_opt(3)*( max(time(:))-min(time(:)) );
orb.a = ( (orb.period/2/pi)^2*mu )^(1/3);

end

