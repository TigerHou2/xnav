function [initGuess,fVal] = bounded_search(fun, r_f0, r_e, r_dM, res)
%BOUNDED_SEARCH searches a given space for fmin at a fixed resolution.

% dat = nan(res);
% 
% for i = 1:res(1)
% for j = 1:res(2)
% for k = 1:res(3)
%     in = [r_f0(i),r_e(j),r_dM(k)];
%     out = fun(in);
%     dat(i,j,k) = out;
% end
% end
% end
% 
% [fVal,idx] = min(dat(:));
% [i,j,k] = ind2sub(size(dat),idx);
% initGuess = [r_f0(i(1)),r_e(j(1)),r_dM(k(1))];



options_fminsearch = optimset(  'Display','none', ...
                                'MaxFunEvals', 50, ...
                                'MaxIter', 50, ...
                                'TolFun',1e-32, ...
                                'TolX', 1e-16);
fVal = Inf;
initGuess = [pi,0.5,pi];
for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
    in = [r_f0(i),r_e(j),r_dM(k)];
    [out,fValThis] = fminsearch(fun,in,options_fminsearch);
    if fValThis < fVal
        fVal = fValThis;
        initGuess = out;
    end
end
end
end

% [Y,X] = meshgrid(r_dM,r_f0);
% Z = reshape(dat(:,j(1),:),res(1),res(3));
% figure;
% scatter3(X(:),Y(:),Z(:),30,Z(:),'filled');
% xlabel('X');
% ylabel('Y');

end

