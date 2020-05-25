function LimSet(lims)
%LIMSET Applies axis limits to current figure
%   Takes axis limit matrix from the LimQuery function

h = gcf;
N = numel(h.Children);
for n = 1:N
    pos1(n) = h.Children(n).Position(1);
    pos2(n) = h.Children(n).Position(2);
end
Ncols = numel(unique(pos1));
Nrows = numel(unique(pos2));

if not( N == size(lims,1)-1 && Ncols == lims(1,1) && Nrows == lims(1,2) )
    error(['Figure layout inconsistent with axes limits format;' newline...
           'Check if current figure has the proper number of subplots.']);
end

for i = 1:N
    h = subplot(Nrows,Ncols,i);
    xlim(h,lims(i+1,1:2));
    ylim(h,lims(i+1,3:4));
end

end

