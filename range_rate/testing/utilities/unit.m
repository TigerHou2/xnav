% Destiny Fawley
% 11/6/2019

function unitVec = unit(v)

[m, n] = size(v);
if m == 1
    v = v';
end

unitVec = zeros(size(v));

I = find(normVectorColumn(v) ~= 0);
unitVec(:,I) = v(:,I)./normVectorColumn(v(:,I));

end