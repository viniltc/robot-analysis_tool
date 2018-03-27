function [traj]=onetraj(traj,col,angles)
if traj.angle>4
    notarget=find(abs(angles-(traj.angle-2*pi)) < 0.1);
else
    notarget=find(abs(angles-traj.angle) < 0.1);
end
the_title=['Target ' num2str(notarget) ];

traj{notarget}=traj.pos;


% figure
% line(traj.pos(:,1),traj.pos(:,2))
% set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')
% title(the_title)