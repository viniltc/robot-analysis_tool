function [teutotjerk,teuprejerk,teupostjerk]=get_jerk(traj)
% Calculates Teulings' index (normalized jerk) for a given trajectory
% Return values are total, acceleration part, deceleration part
% (C) V. Sanguineti, 2008

[S,L] = get_length(traj);
nfjerk = traj.jerk(traj.rind:traj.lind,:);
mjerk =((nfjerk(:,1).^2+nfjerk(:,2).^2))*(traj.tdur.^6)/S.^2;
totjerk = mean(mjerk);

rmsjerk =(nfjerk(:,1).^2+nfjerk(:,2).^2);

teutotjerk = sqrt(totjerk)./2;

firsthalf = 1:(traj.pkind-traj.rind+1);
secondhalf = (traj.pkind-traj.rind+1):(traj.lind-traj.rind);
prejerk =  mean((mjerk(firsthalf,:))).*(traj.adur/traj.tdur);
postjerk = mean((mjerk(secondhalf,:))).*(traj.ddur/traj.tdur);
teuprejerk = sqrt(prejerk+eps)./2;
teupostjerk = sqrt(postjerk+eps)./2;

