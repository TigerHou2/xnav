function init_guess = search_adv(r_f0,r_e,r_dM,res,obsv,pulsar,mu,t_meas,OPT,plot)
%SEARCH_ADV searches a given space for fmin at a fixed resolution.

disp('Searching for initial guess...')

%% coarse mesh to find convex hull of low objective function values

disp('---- Calculating coarse mesh...')

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

disp('---- Done!')

%% find convex hull

disp('---- Finding convex hull...')

% find region of low function values
Fmesh = F(:);
Emesh = E(:);
Mmesh = M(:);
D = dat(:);
lb = quantile(D,0);
ub = quantile(D,0.02); % find the lowest 2% of values
Fmesh = Fmesh(D>lb&D<ub);
Emesh = Emesh(D>lb&D<ub);
Mmesh = Mmesh(D>lb&D<ub);

% create convex hull of the region
ndPoints = [Fmesh,Emesh,Mmesh];
k = convhulln(ndPoints);

% create refined mesh
adv_res = res * 5; % refine the mesh by 5 times in each dimension
adv_r_f0 = linspace(r_f0(1),r_f0(end),adv_res(1));
adv_r_e  = linspace( r_e(1), r_e(end),adv_res(2));
adv_r_dM = linspace(r_dM(1),r_dM(end),adv_res(3));
[aEE,aMM,aFF] = meshgrid(adv_r_e,adv_r_dM,adv_r_f0);
adv_grid = [aFF(:),aEE(:),aMM(:)];

% check for points in convex hull
tol = 1.e-6*mean(abs(ndPoints(:)));
in = inhull(adv_grid,ndPoints,k,tol);
in = reshape(in,adv_res(1),adv_res(2),adv_res(3));

disp('---- Done!')

%% refined search

disp('---- Performing refined search...')

adv_dat = Inf(adv_res);
aF = nan(adv_res(1),adv_res(2),adv_res(3));
aE = nan(adv_res(1),adv_res(2),adv_res(3));
aM = nan(adv_res(1),adv_res(2),adv_res(3));

% iterate over refined mesh and only calculate points within hull
for i = 1:adv_res(1)
for j = 1:adv_res(2)
for k = 1:adv_res(3)
if in(i,j,k)
    fin = [adv_r_f0(i),adv_r_e(j),adv_r_dM(k)];
    out = guess(fin,obsv,pulsar,mu,t_meas);
    adv_dat(i,j,k) = norm(out(:));
    aF(i,j,k) = adv_r_f0(i);
    aE(i,j,k) = adv_r_e(j);
    aM(i,j,k) = adv_r_dM(k);
end
end
end
end

[~,idx] = min(adv_dat(:));
[i,j,k] = ind2sub(size(adv_dat),idx);
init_guess = [adv_r_f0(i(1)),adv_r_e(j(1)),adv_r_dM(k(1))];

disp('---- Done!')
disp('Search complete!')

%% plot results

if ~exist('plot','var')
    return
end

% scatter to see solution space
clf(figure(79))
Fmesh = aF(:);
Emesh = aE(:);
Mmesh = aM(:);
D = adv_dat(:);
lb = quantile(D,0);
ub = quantile(D,0.00005);
Fmesh = Fmesh(D>lb&D<ub);
Emesh = Emesh(D>lb&D<ub);
Mmesh = Mmesh(D>lb&D<ub);
D = D(D>lb&D<ub);
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

