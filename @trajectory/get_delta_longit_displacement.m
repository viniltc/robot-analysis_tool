function [mean_long_dev] = get_delta_longit_displacement(traj)
% Calculates mean longitudinal deviation for a given bimanual trajectory
% (C) V. Sanguineti, 2009

maxrind = max(traj.rind);
minlind = min(traj.lind);

maxrind = max(traj.t1);
minlind = min(traj.t2);
%size(traj.pos)
%traj.rind
%traj.lind
%figure(1)
%clf
%plot(traj.pos(:,1),traj.pos(:,2),'r',traj.pos(:,3),traj.pos(:,4),'b')
%hold on
%plot(traj.pos(traj.rind(1),1),traj.pos(traj.rind(1),2),'r*',...
%     traj.pos(traj.rind(2),3),traj.pos(traj.rind(2),4),'b*')
%plot(traj.pos(traj.lind(1),1),traj.pos(traj.lind(1),2),'ro',...
%     traj.pos(traj.lind(2),3),traj.pos(traj.lind(2),4),'bo')
% drawnow
if maxrind>minlind
    minlind=maxrind;
end
thetraj = traj.pos(maxrind:minlind,:);

d_thetraj = thetraj(:,3:4)-thetraj(:,1:2);

targetdir = [cos(traj.angle); sin(traj.angle)];
mean_long_dev = mean(d_thetraj*targetdir); 
 