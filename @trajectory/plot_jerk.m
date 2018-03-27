function plot_jerk(traj,col)
% plots trajectory jerk (norm) in specified colour
% 
% (C) V. Sanguineti, 2008


jerk = sqrt(traj.jerk(:,1).^2+traj.jerk(:,2).^2);
line(traj.time(traj.interval,1)-traj.time(traj.interval(1),1),jerk(traj.interval),'col',col,'lines','-')

