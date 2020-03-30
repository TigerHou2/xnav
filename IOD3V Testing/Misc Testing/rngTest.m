clear;clc
close all hidden
diff = 0;
sum1 = 0;
sum2 = 0;
for i = 1:1000
    n1 = norm((2*rand(1,3)-1) * 6 / 1000);
    n2   = rand * norm([6,6,6]) / 1000;
    diff = diff + n1 - n2;
    sum1 = sum1 + n1;
    sum2 = sum2 + n2;
end
diff

N1 = [];
N2 = [];
for i = 1:10000
    n1 = (2*rand(1,3)-1) * 6 / 1000;
    N1 = [N1; n1];
    n_r   = rand * norm([6,6,6]) / 1000;
    theta = rand * 2*pi;
    phi   = rand * 2*pi;
    n2 = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)] * n_r;
    N2 = [N2; n2];
end

org = zeros(size(N2));

hold on
quiver3(org(:,1),org(:,2),org(:,3),N1(:,1),N1(:,2),N1(:,3),'LineWidth',0.5)
quiver3(org(:,1),org(:,2),org(:,3),N2(:,1),N2(:,2),N2(:,3),'LineWidth',0.5)
hold off
pbaspect([1,1,1])
axis equal
grid on