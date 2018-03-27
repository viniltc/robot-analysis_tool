function [deltat] = get_final_deltat(traj)
% Calculates the difference between times at which trajectories first enter
% target area...
%
% (C) V. Sanguineti, 2009

% for hand = 1:traj.nhands
%  error(:,1) = traj.xyT(traj.lind(hand),1) - traj.pos(:,2*(hand-1)+1);
%  error(:,2) = traj.xyT(traj.lind(hand),2) - traj.pos(:,2*(hand-1)+2);
%  
%  find(sqrt(error(:,1).^2 + error(:,2).^2) <1e-2)
%  
%  mi(hand) =min( find(sqrt(error(:,1).^2 + error(:,2).^2) <1e-2));
%     
% end
% deltat = (mi(2)-mi(1))/traj.fc;
deltat = traj.rtime(2)+traj.tdur(2) -traj.rtime(1)-traj.tdur(1) 
