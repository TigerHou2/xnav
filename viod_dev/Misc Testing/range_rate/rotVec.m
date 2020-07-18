function v_rot = rotVec(v,k,theta)
    k = k/sqrt(k'*k); % normalize rotation axis
    crosskv = v(:,1); % initialize cross product k and v with right dim.
        crosskv(1) = k(2)*v(3) - k(3)*v(2);
        crosskv(2) = k(3)*v(1) - k(1)*v(3); 
        crosskv(3) = k(1)*v(2) - k(2)*v(1);
        v_rot = cos(theta)*v + (crosskv)*sin(theta)...
                        + k*(k'*v)*(1 - cos(theta));
end
