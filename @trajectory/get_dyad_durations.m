function [ts1,ts2,te1,te2]=get_duration(traj)
% Calculates significant movement-related durations for a dyad
% (C) V. Sanguineti, TC Vinil 2015
%keyboard
ts1 = traj.time(traj.rind(1));
ts2 = traj.time(traj.rind(2));
te1 = traj.time(traj.lind(1));
te2 = traj.time(traj.lind(2));
