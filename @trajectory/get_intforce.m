function [forcenorm]=get_intforce(traj,par)


forcenorm = mean(sqrt(traj.intforce(:,1).^2+traj.intforce(:,2).^2)); 

