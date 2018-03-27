function [early_cacc,late_cacc]=get_centripetal_acceleration(traj)
% Calculates average acceleration in orthogonal direction(a measure of tremor)
% in acceleration and deceleration parts of trajectory
% (C) V. Sanguineti, 2008

dir_vers = traj.vel./(eps+sqrt(traj.vel(:,1).^2+traj.vel(:,2).^2)*ones(1,2));
tacc = traj.acc(:,1).*dir_vers(:,1)+traj.acc(:,2).*dir_vers(:,2); % vector quantity of the tangential accleration
pacc = traj.acc - (tacc*ones(1,2)).*dir_vers; % Residual acceleration
cacc = -pacc(:,1).*dir_vers(:,2)+pacc(:,2).*dir_vers(:,1); %Centripetal Acceleration

early_cacc = mean(cacc(traj.rind:traj.pkind));
late_cacc = mean(cacc(traj.pkind:traj.lind));