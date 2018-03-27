function [lat_dev,n_lat_dev] = get_lateral_deviation(traj)
% Calculates maximum lateral deviation (absolute, normalized) for a given trajectory
% (C) V. Sanguineti, 2008

%rangle = traj.rotation;

%R = [cos(rangle) -sin(rangle);...
%     sin(rangle) cos(rangle)];

direz = [cos(traj.angle) sin(traj.angle)];

n_direz = [-direz(2) direz(1)];



thetraj = traj.pos(traj.rind:traj.lind,:);

% correct for visual rotation: want to calculate latdev of cursor trajectory not
% hand trajectory
%thetraj = traj.pos(traj.rind:traj.lind,:)*R';
 
% project traj over direz
t_norm = (sum((n_direz(1:2)'*ones(1,length(thetraj))).*thetraj'))';

[m,im]=max(abs(t_norm));
lat_dev = m*sign(t_norm(im));


ampl = sqrt((thetraj(end,1)-thetraj(1,1)).^2+ (thetraj(end,2)-thetraj(1,2)).^2);
n_lat_dev = lat_dev/ampl;
 