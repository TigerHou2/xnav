% Destiny Fawley
% 11/6/2019

function V = quatVectorRotation(q, v)

[m, n] = size(q);
if m ~= 4 && n ~= 4
    error('Quaternion must be given as [4 X N] or [N X 4].');
end

[p, r] = size(v);
if p ~= 3 && r ~= 3
    error('Vector must be given as [3 X N] or [N X 3].');
end

if m == 4
    q = q';
end
if p == 3
    v = v.';
end

V = quatRotate(q,v);

if p == 3
    V = V.';
end


end

