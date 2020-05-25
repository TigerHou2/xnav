function disp_orbit(R,V,mu,res,dur,start_day,ref)
%DISP_ORBIT Summary of this function goes here
%   Detailed explanation goes here

if nargin == 7 && ref == 1
    color = 'Red';
else
    color = 'Blue';
end

pos = Get_Orb_Points(R(end,:),V(end,:),mu,res,dur,start_day);
plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',1,'Color',color);

view(cross(R(end,:),V(end,:)))
pbaspect([1 1 1])
axis equal
grid on

end

