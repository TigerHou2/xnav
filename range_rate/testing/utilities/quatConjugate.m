% Destiny Fawley
% 11/6/2019

function quatConj = quatConjugate(q)
% Takes conjugate of scalar-first quaternion q. Can be given as [4 X N] or
% [N X 4]. A [4 X 4] is assumed to have quaternion columns.

[m, n] = size(q);

% chekc that one dimension is 4

if m ~= 4 && n ~= 4
    error('Quaternion must be given as [4 X N] or [N X 4].');
end

if m == 4
    quatConj = q;
    quatConj(2:end,:) = -quatConj(2:end,:);
else
    quatConj = q;
    quatConj(:,2:end) = -quatConj(:,2:end);
    
end

