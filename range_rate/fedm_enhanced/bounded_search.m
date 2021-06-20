function initGuess = bounded_search(fun, r_f0, r_e, r_dM, res)
%BOUNDED_SEARCH searches a given space for fmin at a fixed resolution.

dat = nan(res);

for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
    in = [r_f0(i),r_e(j),r_dM(k)];
    out = fun(in);
    dat(i,j,k) = norm(out(:));
end
end
end

[~,idx] = min(dat(:));
[i,j,k] = ind2sub(size(dat),idx);
initGuess = [r_f0(i(1)),r_e(j(1)),r_dM(k(1))];

end

