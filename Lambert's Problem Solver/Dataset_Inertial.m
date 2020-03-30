close all
clear;clc
Parameters

% tolerance for position error
% this tolerance is proportional to the distance between points
%   to ensure closer points require higher resolution
tol = 0.01;

r = [1 0.2 0.05] * AU;
v = [8, 30, 7];

t1 = 10;
t2 = 50;
t3 = 150;

[r1, v1] = TimeProp(r, v, sun.mu, t1);
[r2, v2] = TimeProp(r, v, sun.mu, t2);
[r3, v3] = TimeProp(r, v, sun.mu, t3);

% add noise to velocity 'measurements'
% assuming variance to be proportional to speed
%   don't know what the error order of magnitude is
%   so we will use this as a test to determine sensor requirements
noise_factor = 0.001;
v1 = v1 + randn(1,3) .* (noise_factor * v1);
v2 = v2 + randn(1,3) .* (noise_factor * v2);
v3 = v3 + randn(1,3) .* (noise_factor * v3);

% randomly generate fake position vectors for the three velocities
num_fakes = 500;
r1_set = (rand(num_fakes,3)-0.5)*5*AU;
r1_set = [r1; r1_set];
r2_set = (rand(num_fakes,3)-0.5)*5*AU;
r2_set = [r2; r2_set];
r3_set = (rand(num_fakes,3)-0.5)*5*AU;
r3_set = [r3; r3_set];

% grab a copy of the sets in the absolute frame
r1_set_abs = r1_set;
r2_set_abs = r2_set;
r3_set_abs = r3_set;

% generate rotation
% rx = rand*90;
% ry = rand*90;
% rz = rand*90;
% R = rotz(rz)*roty(ry)*rotx(rx);

% no rotation: we can use the stars for attitude determination!
R = eye(3);

% generate translation
L = (rand(1,3)-0.5)*0.5*AU;

r1_set = (R*r1_set')' + L;
r2_set = (R*r2_set')' + L;
r3_set = (R*r3_set')' + L;

% display the true solution
r_true = [r1;r2;r3];
disp(r_true)

% solve the 3V IOD problem
[r1,r2,r3,K] = IOD3V(v1,v2,v3,sun.mu,'omega',[0,0,1],'prograde',true);

% use results of IOD3V to parse candidate points
r12 = r2 - r1;
r23 = r3 - r2;

counter(1).r = [];
counter(1).L = [];
sol_counter = 1;

% permutate the three sets and compare orbit to ground truth
for i = 1:num_fakes+1
% iterate over every point in first set

    for j = 1:num_fakes+1
    % look for points in a spherical shell centered on the first point
    
        r12_cand = r2_set(j,:) - r1_set(i,:);
        if abs(norm(r12_cand)-norm(r12))<tol*norm(r12)
            for k = 1:num_fakes+1
            % look for point in a sphere in a fixed vector
            % from the first two points
            
                r23_cand = r3_set(k,:) - r2_set(j,:);
                if abs(norm(r23_cand)-norm(r23)) < ...
                        tol * norm(r23) && ...
                   abs((dot(r23_cand,r12_cand)/norm(r23_cand)/norm(r12_cand) - ...
                        dot(r23,r12)/norm(r23)/norm(r12))) < ...
                        tol * abs(dot(r23,r12)/norm(r23)/norm(r12))
                    
                    counter(sol_counter).r(1,:) = r1_set(i,:);
                    counter(sol_counter).r(2,:) = r2_set(j,:);
                    counter(sol_counter).r(3,:) = r3_set(k,:);
                    
                    % find the translation between two frames
                    counter(sol_counter).L = mean([...
                        counter(sol_counter).r(1,:) - r1;
                        counter(sol_counter).r(2,:) - r2;
                        counter(sol_counter).r(3,:) - r3]);
                    
                    counter(sol_counter).r = ...
                        counter(sol_counter).r - counter(sol_counter).L;
                    sol_counter = sol_counter + 1;
                    
                    disp('Finished with valid solution!')
                end
            end
        end
    end
end

for i = 1:sol_counter-1
    disp(counter(i).r)
end

% plot in absolute frame
figure
hold on
% show all fake points
plot3(r1_set_abs(:,1), r1_set_abs(:,2), r1_set_abs(:,3), 'ro')
plot3(r2_set_abs(:,1), r2_set_abs(:,2), r2_set_abs(:,3), 'ko')
plot3(r3_set_abs(:,1), r3_set_abs(:,2), r3_set_abs(:,3), 'bo')
% show discovered solutions
for i = 1:sol_counter-1
    r = counter(i).r;
    plot3(r(1,1), r(1,2), r(1,3), 'ro','MarkerFaceColor','r')
    plot3(r(2,1), r(2,2), r(2,3), 'ko','MarkerFaceColor','k')
    plot3(r(3,1), r(3,2), r(3,3), 'bo','MarkerFaceColor','b')
end
% overlay true solution
plot3(r_true(:,1), r_true(:,2), r_true(:,3), 'g+')
plot3(r_true(:,1), r_true(:,2), r_true(:,3), 'g+')
plot3(r_true(:,1), r_true(:,2), r_true(:,3), 'g+')
% fix axes
pbaspect([1 1 1])
axis equal
grid on
hold off

% % plot in inertial frame
% figure
% hold on
% plot3(r1_set(:,1), r1_set(:,2), r1_set(:,3), 'ro')
% plot3(r2_set(:,1), r2_set(:,2), r2_set(:,3), 'ko')
% plot3(r3_set(:,1), r3_set(:,2), r3_set(:,3), 'bo')
% pbaspect([1 1 1])
% axis equal
% grid on
% hold off