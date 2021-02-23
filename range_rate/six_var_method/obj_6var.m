function obj = obj_6var(guess,mu,obsv,pulsars,time)
%OBJ_6VAR objective function for the 6-variable RROD technique

R = guess(1);
e = guess(2);
f0 = guess(3);

c = e*R;
n = 1/mu*(R^2-c^2)^(3/2);

E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
M0 = E0 - e*sin(E0);
M = n * (time - time(1)) + M0;
E = kepler(M,e);
F = 2 * atan( tan(E/2) * sqrt((1+e)/(1-e)) );

rotm = eul2rotm(guess(4:6));

P = rotm * pulsars;
A = atan2(P(2,:),P(1,:))';
B = atan2(P(3,:),vecnorm(P(1:2,:),2,1))';

obj = R*cos(B).*(sin(A-F)+e*sin(A)) - obsv;

obj = reshape(obj,size(obsv));

end

