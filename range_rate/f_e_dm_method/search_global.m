function init_guess = search_global(r_f0,r_e,r_dM,res,obsv,pulsar,mu,t_meas,OPT,plot)
%SEARCH_GLOBAL searches a given space for fmin at a fixed resolution.

disp('Searching for initial guess...')

dat = nan(res);

F = nan(res(1),res(2),res(3));
E = nan(res(1),res(2),res(3));
M = nan(res(1),res(2),res(3));

for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
    fin = [r_f0(i),r_e(j),r_dM(k)];
    out = guess(fin,obsv,pulsar,mu,t_meas);
    dat(i,j,k) = norm(out(:));
    F(i,j,k) = r_f0(i);
    E(i,j,k) = r_e(j);
    M(i,j,k) = r_dM(k);
end
end
end

% mins = islocalmin(dat);
% [~,idx] = min(dat(:)./(mins(:)+0.01));
[~,idx] = min(dat(:));
[i,j,k] = ind2sub(size(dat),idx);
init_guess = [r_f0(i(1)),r_e(j(1)),r_dM(k(1))];

disp('Search complete!')

if ~exist('plot','var')
    return
end

% scatter to see solution space
clf(figure(79))
Fmesh = F(:);
Emesh = E(:);
Mmesh = M(:);
D = dat(:);
lb = quantile(D,0);
ub = quantile(D,0.002);
Fmesh = Fmesh(D>lb&D<ub);
Emesh = Emesh(D>lb&D<ub);
Mmesh = Mmesh(D>lb&D<ub);
D = D(D>lb&D<ub);
mins = islocalmin(dat);
Floc = F(mins);
Eloc = E(mins);
Mloc = M(mins);
hold on
scatter3(OPT(1),OPT(2),OPT(3),24,'cyan',...
        'DisplayName','True Soln',...
        'LineWidth',1.5)
scatter3( init_guess(1), init_guess(2), init_guess(3),24,'green',...
        'DisplayName','Init Guess',...
        'LineWidth',1.5)
fig = scatter3(Fmesh,Emesh,Mmesh,18,D,'filled',...
        'DisplayName','Objective Function');
fig.MarkerFaceAlpha = 0.6;
% scatter3(Floc(:),Eloc(:),Mloc(:),10,'red','filled',...
%         'DisplayName','Local Minima')
hold off

colormap(bone(200))
xlabel('f0')
ylabel('e')
zlabel('dM')
pbaspect([1 1 1])
view([1 1 1])
colorbar
legend('Location','Best')

end

