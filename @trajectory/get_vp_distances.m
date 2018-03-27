function [tc1,tc2,md1,md2]=get_vp_distances(traj,par)
% Calculates time of crossing and minimum distance from VP for partner
% 'par' of a dyad
% (C) TC Vinil 2015
% keyboard

% interval_common = max(traj.interval{1}(1),traj.interval{2}(1)):min(traj.interval{1}(end),traj.interval{2}(end));
[md1,i1]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,2)).^2+...
                   (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(1,3)).^2));
[md2,i2]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(2,2)).^2+...
                   (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,3)).^2));
               

%   keyboard
               
tc1 = traj.time(traj.interval{par}(i1));
tc2 = traj.time(traj.interval{par}(i2));