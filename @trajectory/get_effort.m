function [eff1,eff2]=get_effort(traj,par)
% Calculates effort (mean squared force) for partner
% 'par' of a dyad
% (C) TC Vinil 2015
interval_common = max(traj.interval{1}(1),traj.interval{2}(1)):min(traj.interval{1}(end),traj.interval{2}(end));

% eff1 = sqrt(mean(sum(traj.intforce(traj.interval{1},2:3).^2,2)));
% eff2 = sqrt(mean(sum(traj.intforce(traj.interval{2},5:6).^2,2)));
    
eff1 = sqrt(mean(sum(traj.intforce(interval_common,2:3).^2,2)));
eff2 = sqrt(mean(sum(traj.intforce(interval_common,5:6).^2,2)));

%eff1= mean(traj.intforce(traj.interval{1}).^2);
%eff2= mean(traj.intforce(traj.interval{2}).^2);