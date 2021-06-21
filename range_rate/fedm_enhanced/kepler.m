function E = kepler(M,e)
%KEPLER solves Kepler's equation using the Conway-Laguerre algorithm.

d = 1;
n = 5;

if e < 1
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(E-e*sin(E)-M) ./ ...
            (1-e*cos(E) + sign(1-e*cos(E)) ...
                       .* sqrt(abs((n-1)^2*(1-e*cos(E)).^2 ...
                                   -n*(n-1)*(E-e*sin(E)-M).*(e*sin(E)))));
        d = max(abs(E-e*sin(E)-M),[],'all');
    end
else
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(e*sinh(E)-E-M) ./ ...
            (e*cosh(E)-1 + sign(e*cosh(E)-1) ...
                        .* sqrt(abs((n-1)^2*(e*cosh(E)-1).^2 ...
                                    -n*(n-1)*(e*sinh(E)-E-M).*(e*sinh(E)))));
        d = max(abs(e*sinh(E)-E-M),[],'all');
    end
end

end