function points = Cut_Noodle(r,v,mu,DT,noise,start_day,reps,samples)
%CUT_NOODLE Summary of this function goes here
%   Detailed explanation goes here

points = [];
prog_bar = waitbar(0,'Noodle Progress: 0%');

for k = 1:samples
    waitbar(k/samples,prog_bar,...
            ['Noodle Progress: ' num2str(k/samples*100,'%5.2f') '%']);
    
    t = zeros(1,reps);
    dt = DT / (reps-1);
    for i = 1:reps
        t(i) = start_day + (i-1)*dt;
    end

    R = zeros(reps,3);
    V = zeros(reps,3);

    for i = 1:reps
        [R(i,:),V(i,:)] = TimeProp_Universal_V2(r,v,mu,t(i));
    end

    ri = R(1,:);
    vi = V(1,:);
    rf = R(end,:);
    vf = V(end,:);

    % add noise to velocity 'measurements'
    % assuming variance to be proportional to speed
    %   don't know what the error order of magnitude is
    %   so we will use this as a test to determine sensor requirements
    for i = 1:reps
        n_r   = rand * noise / 1000;
        theta = rand * 2*pi;
        phi   = rand * 2*pi;
        n_vec = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)] * n_r;
        V(i,:) = V(i,:) + n_vec;
    end

    % solve the 3V IOD problem
    R = IOD3V_Ext(V,mu,'ordered',true);
    rf = R(end,:);
    points = [points; rf];
end

scatter3(points(:,1),points(:,2),points(:,3),28,'b');
axis equal
grid on
pbaspect([1,1,1])

end

