function [max_horiz_dev,mean_horiz_dev] = get_delta_horiz_deviation(traj)
% Calculates maximum and mean horizontal deviation for a given bimanual trajectory
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
[ma, mi]=max(abs(thetraj(:,3)-thetraj(:,1)));

max_horiz_dev = thetraj(mi,3)-thetraj(mi,1);
max(thetraj(:,3)-thetraj(:,1));
mean_horiz_dev = mean(thetraj(:,3)-thetraj(:,1)); 
 