function lims = LimQuery
%LIMQUERY Gets current figure axis limits
%   Stores the subplot configuration and axis limits in a matrix

h = gcf;
N = numel(h.Children);
for n = 1:N
    pos1(n) = h.Children(n).Position(1);
    pos2(n) = h.Children(n).Position(2);
end
Ncols = numel(unique(pos1));
Nrows = numel(unique(pos2));
lims = [Ncols, Nrows, 0, 0];

for i = 1:N
    h = subplot(Nrows,Ncols,i);
    lims = [lims; [xlim ylim]];
end

end