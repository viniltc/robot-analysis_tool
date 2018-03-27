function [psp, asp] = get_speed(traj)
% Computes peak and average speed
% (C) V. Sanguineti, 2008
 speed = sqrt(traj.vel(:,1).^2+traj.vel(:,2).^2);
 psp = traj.pkspeed;
 asp = mean(speed(traj.rind:traj.lind));

%  speed1 = mean(sqrt(sum(traj.vel(traj.interval{1},2:3).^2,2)));
%  speed2 = mean(sqrt(sum(traj.vel(traj.interval{2},5:6).^2,2)));
