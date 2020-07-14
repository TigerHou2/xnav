function Sort_Orbit(ell,hyp)
%SORT_ORBIT Summary of this function goes here
%   Detailed explanation goes here

% count the number of elliptic and hyperbolic orbits for each interval
if isempty(ell)
    elliptic = [];
else
    elliptic = ell(:,1)';
end
if isempty(hyp)
    hyperbolic = [];
else
    hyperbolic = hyp(:,1)';
end
% sort the intervals and create bins for the histogram
bins = unique([elliptic hyperbolic]);
bins_hist = [bins bins(end)];
% count the occurrences of each bin for elliptic and hyperbolic
E = histcounts(elliptic,bins_hist);
H = histcounts(hyperbolic,bins_hist);
% format for the bar() function
bars = [E; H]';
% create figure
figure('Name','Number of Elliptic vs. Hyperbolic Orbits')
bar(bins,bars,'stacked')
xlabel('Measurement Interval (days)')
ylabel('Number of Cases')
legend('Elliptic','Hyperbolic')

end

