function [S, L] = get_length(traj)
% Calculates length of trajectory (S) and movement amplitude (L)
% (C) V. Sanguineti, 2008

nfvel = traj.vel(traj.rind:traj.lind,:);
nfspeed = sqrt(nfvel(:,1).^2+nfvel(:,2).^2);

dt = mean(diff(traj.time));
S = sum(nfspeed)*dt;

%l'indice di linearità è dato dal rapporto tra lunghezza della traiettoria (s) e la
%lunghezza della congiungente tra punto iniziale e finale (L) meno 1
% indice = s/L - 1
start_traj = traj.pos(traj.rind,:);
end_traj = traj.pos(traj.lind,:);
L = sqrt((end_traj(1)-start_traj(1)).^2+ (end_traj(2)-start_traj(2)).^2);
