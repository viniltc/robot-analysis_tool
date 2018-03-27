function [tdur,adur,ddur,rtime,rind,lind]=get_duration(traj)
% Calculates significant movement-related durations 
% Return values are total duration, acceleration duration, 
% deceleration duration, reaction time, and 
% index of movement onset and movement end
% (C) V. Sanguineti, 2008

tdur = traj.tdur;
adur = traj.adur;
ddur = traj.ddur;

rtime = traj.rtime;
rind = traj.rind;
lind = traj.lind;