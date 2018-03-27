function [trajtgt,notarget,t]=get_onetraj(traj,tgt,angles)
if traj.angle>4
    notarget=find(abs(angles-(traj.angle-2*pi)) < 0.1);
else
    notarget=find(abs(angles-traj.angle) < 0.1);
end
 
trajtgt=[cos(2*pi-traj.angle+traj.rotation(1)) -sin(2*pi-traj.angle+traj.rotation(1));sin(2*pi-traj.angle+traj.rotation(1)) cos(2*pi-traj.angle+traj.rotation(1)) ]*traj.pos';
% traj.rotation*pi/6
tgtr=[cos(2*pi-traj.angle+traj.rotation(1)) -sin(2*pi-traj.angle+traj.rotation(1));sin(2*pi-traj.angle+traj.rotation(1)) cos(2*pi-traj.angle+traj.rotation(1)) ]*tgt';
t=tgtr(:,notarget)'

% figure
% line(traj.pos(:,1),traj.pos(:,2))
% set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')
% title(the_title)