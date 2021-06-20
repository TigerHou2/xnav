% Destiny Fawley
% 11/6/2019

function fig = drawRocket(result, rocket, models, varargin)

p = inputParser;
addOptional(p,'VideoFile','NewMovie',@ischar);
addOptional(p,'makeMovie',false,@islogical);

parse(p,varargin{:});

args = p.Results;
if args.makeMovie
    v = VideoWriter(args.VideoFile);
    v.FrameRate = models.dataRate; % real time
    open(v);    
end

% Clear the current figure.
figure(1);
clf;
fig.text.axis = axes('position',[0 0 1 1]);
axis([0 1 0 1]);
hold on
% axis off
grid on
fs = 14;

velocity = norm(result.traj.velI(:,1));
alt = result.traj.posI(3,1);

fig.text.vel = text(0.05,0.975,...
    sprintf('Velocity: %0.2f ',velocity),...
    'fontweight','bold','fontsize',fs,...
    'color','k','verticalalignment','top');
fig.text.alt = text(0.05,0.9,...
    sprintf('Altitude: %0.2f',alt),...
    'fontweight','bold','fontsize',fs,...
    'color','k','verticalalignment','top');

set(gcf,'renderer','opengl');
set(gcf,'color','w');

fig.view0.axis = axes('position',[-0.75 -0.2 2.5 1.5]);
grid on
axis equal;
fig.view0.dx = .5;
fig.view0.dy = .5;
fig.view0.dz1 = .8;
fig.view0.dz2 = .4;
set(fig.view0.axis,'xlim',[-10 10]);%[result.traj.posI(1,1)-fig.view0.dx, result.traj.posI(1,1)+fig.view0.dx]);
set(fig.view0.axis,'ylim',[-10 10]);%[result.traj.posI(2,1)-fig.view0.dy, result.traj.posI(2,1)+fig.view0.dy]);
set(fig.view0.axis,'zlim',[0, 15]);
% set(fig.view0.axis,'zlim',[0-fig.view0.dz, result.traj.posI(3,1)+fig.view0.dz]);
axis manual;
hold on;
axis off;
box on;

% view([180-37.5,20]);
view([1 1 1]);
set(gca,'projection','perspective');
set(gca,'clipping','on','clippingstyle','3dbox');
lighting gouraud
fig.view0.light = light('position',[-1;1;1],'style','local');

[pointsRocket, facesRocket, cDataRocket, pointsThruster, facesThruster, cDataThruster] = getRocketModel(rocket);

%%
for i = 1:length(result.traj.time)-1
    
    xvec = quatVectorRotation(quatConjugate(result.traj.qi2b(:,i)), [1 0 0]);
    yvec = quatVectorRotation(quatConjugate(result.traj.qi2b(:,i)), [0 1 0]);
    zvec = quatVectorRotation(quatConjugate(result.traj.qi2b(:,i)), [0 0 2]);
    
    X = plot3([result.traj.posI(1,i), xvec(1)+result.traj.posI(1,i)], ...
        [result.traj.posI(2,i), xvec(2)+result.traj.posI(2,i)], ...
        [result.traj.posI(3,i), xvec(3)+result.traj.posI(3,i)], 'r->');
    Y = plot3([result.traj.posI(1,i), yvec(1)+result.traj.posI(1,i)], ...
        [result.traj.posI(2,i), yvec(2)+result.traj.posI(2,i)], ...
        [result.traj.posI(3,i), yvec(3)+result.traj.posI(3,i)], 'g->');
    Z = plot3([result.traj.posI(1,i), zvec(1)+result.traj.posI(1,i)], ...
        [result.traj.posI(2,i), zvec(2)+result.traj.posI(2,i)], ...
        [result.traj.posI(3,i), zvec(3)+result.traj.posI(3,i)], 'b->');
    
    velocity = norm(result.traj.velI(:,i));
    if result.traj.velI(3,i) > 1
        velSign = 1;
    else
        velSign = -1;
    end
    alt = result.traj.posI(3,i);
    
    
%     set(fig.view0.axis,'xlim',[-1 1]);%[result.traj.posI(1,1)-fig.view0.dx, result.traj.posI(1,1)+fig.view0.dx]);
%     set(fig.view0.axis,'ylim',[-1 1]);%[result.traj.posI(2,1)-fig.view0.dy, result.traj.posI(2,1)+fig.view0.dy]);
%     set(fig.view0.axis,'zlim',[0, 15]);
    
    set(fig.view0.axis,'xlim',[result.traj.posI(1,i)-fig.view0.dx, result.traj.posI(1,i)+fig.view0.dx]);
    set(fig.view0.axis,'ylim',[result.traj.posI(2,i)-fig.view0.dy, result.traj.posI(2,i)+fig.view0.dy]);
    set(fig.view0.axis,'zlim',[result.traj.posI(3,i)-fig.view0.dz2, result.traj.posI(3,i)+fig.view0.dz1]);
    
    % - transformations
    pointsRocketI = quatVectorRotation(quatConjugate(result.traj.qi2b(:,i)), pointsRocket);
    pointsRocketI(1,:) = pointsRocketI(1,:) + result.traj.posI(1,i);
    pointsRocketI(2,:) = pointsRocketI(2,:) + result.traj.posI(2,i);
    pointsRocketI(3,:) = pointsRocketI(3,:) + result.traj.posI(3,i);
    
    pointsThruster(3,1) = -norm(result.traj.thrustI(:,i)/100);
    % Thrust Vector
    R1 = [1 0 0; 
          0 cos(result.ctrl.beta(i)) sin(result.ctrl.beta(i));
          0 -sin(result.ctrl.beta(i)) cos(result.ctrl.beta(i))];
    R2 = [cos(result.ctrl.alpha(i)) 0 sin(result.ctrl.alpha(i));
          0 1 0;
          -sin(result.ctrl.alpha(i)) 0 cos(result.ctrl.alpha(i))];
    R = R1*R2;
%     if i > 115
%         i
%     end
    pointsThrusterB = pointsThruster;
    pointsThrusterB(:,1) = R*pointsThrusterB(:,1);
    pointsThrusterI = quatVectorRotation(quatConjugate(result.traj.qi2b(:,i)), pointsThrusterB);
    pointsThrusterI(1,:) = pointsThrusterI(1,:) + result.traj.posI(1,i);
    pointsThrusterI(2,:) = pointsThrusterI(2,:) + result.traj.posI(2,i);
    pointsThrusterI(3,:) = pointsThrusterI(3,:) + result.traj.posI(3,i);
    
    fig.rocket = patch('Vertices',pointsRocketI','Faces',facesRocket,'FaceColor','flat',...
        'FaceVertexCData',cDataRocket,'FaceAlpha',1,'EdgeAlpha',0,...
        'backfacelighting','reverselit','AmbientStrength',0.6);
    
    if norm(result.traj.thrustI(:,i)) ~= 0
        fig.thruster =  patch('Vertices',pointsThrusterI','Faces',facesThruster,'FaceColor','flat',...
        'FaceVertexCData',cDataThruster,'FaceAlpha',1,'EdgeAlpha',0,...
        'backfacelighting','reverselit','AmbientStrength',0.6);
    end
    
    set(fig.text.vel,'string',sprintf('Velocity: %0.2f m/s',velSign*velocity));
    set(fig.text.alt,'string',sprintf('Altitude: %0.2f m',alt));
    drawnow;
    pause(0.01)
    
    if args.makeMovie
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
    delete(fig.rocket);
    delete(X);
    delete(Y);
    delete(Z);
    
    if isfield(fig,'thruster')
        delete(fig.thruster);
    end
end

% final iteration
velocity = norm(result.traj.velI(:,end));
    if result.traj.velI(3,end) > 1
        velSign = 1;
    else
        velSign = -1;
    end
    alt = result.traj.posI(3,end);
    
    
%     set(fig.view0.axis,'xlim',[-1 1]);%[result.traj.posI(1,1)-fig.view0.dx, result.traj.posI(1,1)+fig.view0.dx]);
%     set(fig.view0.axis,'ylim',[-1 1]);%[result.traj.posI(2,1)-fig.view0.dy, result.traj.posI(2,1)+fig.view0.dy]);
%     set(fig.view0.axis,'zlim',[0, 15]);
    
    set(fig.view0.axis,'xlim',[result.traj.posI(1,end)-fig.view0.dx, result.traj.posI(1,end)+fig.view0.dx]);
    set(fig.view0.axis,'ylim',[result.traj.posI(2,end)-fig.view0.dy, result.traj.posI(2,end)+fig.view0.dy]);
    set(fig.view0.axis,'zlim',[result.traj.posI(3,end)-fig.view0.dz2, result.traj.posI(3,end)+fig.view0.dz1]);
%     
    % - transformations
    pointsRocketI = quatVectorRotation(quatConjugate(result.traj.qi2b(:,end)), pointsRocket);
    pointsRocketI(1,:) = pointsRocketI(1,:) + result.traj.posI(1,end);
    pointsRocketI(2,:) = pointsRocketI(2,:) + result.traj.posI(2,end);
    pointsRocketI(3,:) = pointsRocketI(3,:) + result.traj.posI(3,end);
    
    pointsThruster(3,1) = -norm(result.traj.thrustI(:,end)/100);
    pointsThrusterI = quatVectorRotation(quatConjugate(result.traj.qi2b(:,end)), pointsThruster);
    pointsThrusterI(1,:) = pointsThrusterI(1,:) + result.traj.posI(1,end);
    pointsThrusterI(2,:) = pointsThrusterI(2,:) + result.traj.posI(2,end);
    pointsThrusterI(3,:) = pointsThrusterI(3,:) + result.traj.posI(3,end);
    
    fig.rocket = patch('Vertices',pointsRocketI','Faces',facesRocket,'FaceColor','flat',...
        'FaceVertexCData',cDataRocket,'FaceAlpha',1,'EdgeAlpha',0,...
        'backfacelighting','reverselit','AmbientStrength',0.6);
    
    if norm(result.traj.thrustI(:,end)) ~= 0
        fig.thruster =  patch('Vertices',pointsThrusterI','Faces',facesThruster,'FaceColor','flat',...
        'FaceVertexCData',cDataThruster,'FaceAlpha',1,'EdgeAlpha',0,...
        'backfacelighting','reverselit','AmbientStrength',0.6);
    end
    
    set(fig.text.vel,'string',sprintf('Velocity: %0.2f m/s',velSign*velocity));
    set(fig.text.alt,'string',sprintf('Altitude: %0.2f m',alt));
    drawnow;
    
    if args.makeMovie
        frame = getframe(gcf);
        writeVideo(v,frame);
        close(v);
    end
end