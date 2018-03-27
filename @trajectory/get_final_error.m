function [error_norm, err] = get_final_error(traj)
% Calculates the difference between target and final position
%
% (C) V. Sanguineti, 2008

error = traj.xyT(traj.lind,:) - traj.pos(traj.lind,:);
error_norm = sqrt(err(1).^2+err(2).^2);
