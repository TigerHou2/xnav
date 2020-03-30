close all
clear;clc
Parameters

r = [1 0.2 0.05] * AU;
v = [8, 30, 7];

t1 = 10;
t2 = 50;
t3 = 150;

[r1, v1] = TimeProp(r, v, sun.mu, t1);
[r2, v2] = TimeProp(r, v, sun.mu, t2);
[r3, v3] = TimeProp(r, v, sun.mu, t3);

num_fakes = 30;

% randomly generate fake position vectors for the three velocities
r1_set = (rand(num_fakes,3)-0.5)*5*AU;
r1_set = [r1; r1_set];
r2_set = (rand(num_fakes,3)-0.5)*5*AU;
r2_set = [r2; r2_set];
r3_set = (rand(num_fakes,3)-0.5)*5*AU;
r3_set = [r3; r3_set];

r_sol = zeros(3,3);

% permutate the three sets and compare orbit to ground truth
for i = 1:num_fakes+1
    if norm(r1_set(i,:)-r1) < 0.05*norm(r1)
        r_sol(1,:) = r1_set(i,:);
        for j = 1:num_fakes+1
            if norm(r2_set(j,:)-r2) < 0.05*norm(r2)
                r_sol(2,:) = r2_set(j,:);
                for k = 1:num_fakes+1
                    if norm(r3_set(k,:)-r3) < 0.05*norm(r3)
                        r_sol(3,:) = r3_set(k,:);
                        disp('Finished with valid solution!')
                    end
                end
            end
        end
    end
end

r_true = [r1;r2;r3];
disp(r_true)
disp(r_sol)

hold on
plot3(r1_set(:,1), r1_set(:,2), r1_set(:,3), 'ro')
plot3(r2_set(:,1), r2_set(:,2), r2_set(:,3), 'ko')
plot3(r3_set(:,1), r3_set(:,2), r3_set(:,3), 'bo')
plot3(r_true(1,1), r_true(1,2), r_true(1,3), 'ro','MarkerFaceColor','r')
plot3(r_true(2,1), r_true(2,2), r_true(2,3), 'ko','MarkerFaceColor','k')
plot3(r_true(3,1), r_true(3,2), r_true(3,3), 'bo','MarkerFaceColor','b')
pbaspect([1 1 1])
axis equal
grid on
hold off