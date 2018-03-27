%function [pow1,pow2,tp1,tp2] = get_power(traj,par)
function [pow1,pow2] = get_power(traj,par)

% interval_common = max(traj.interval{1}(1),traj.interval{2}(1)):min(traj.interval{1}(end),traj.interval{2}(end));
% pow1= mean(sum(traj.vel(interval_common,2:3).*traj.intforce(interval_common,2:3)));
% pow2= mean(sum(traj.vel(interval_common,5:6).*traj.intforce(interval_common,5:6)));


% keyboard
interval = max(traj.interval{1}(1),traj.interval{2}(1)):min(traj.interval{1}(end),traj.interval{2}(end));

pow1 = mean(sumabs(traj.vel(interval,2:3).*traj.intforce(interval,2:3))); %mean value on trial
pow2 = mean(sumabs(traj.vel(interval,5:6).*traj.intforce(interval,5:6)));

% pow1 = mean(sum(traj.vel(traj.interval{1},2:3)).*sum(traj.intforce(traj.interval{1},2:3)));
% pow2 = mean(sum(traj.vel(traj.interval{2},5:6)).*sum(traj.intforce(traj.interval{2},5:6)));
% 
% pow1 = mean(sumabs(traj.vel(traj.interval{1},2:3).*traj.intforce(traj.interval{1},2:3)));
% pow2 = mean(sumabs(traj.vel(traj.interval{2},5:6).*traj.intforce(traj.interval{2},5:6)));
% keyboard
%  tp1 = traj.time(interval_common);
%  tp2 = traj.time(interval_common);
%  tp1 = find(traj.time(traj.vel(interval_common,2:3).*traj.intforce(interval_common,2:3))>0);
%  tp2 = find(traj.time(traj.vel(interval_common,5:6).*traj.intforce(interval_common,5:6))>0);
%  p1 = traj.vel(interval_common,2:3)*traj.intforce(interval_common,2:3)
%  p2 = traj.vel(interval_common,5:6)*traj.intforce(interval_common,5:6)
%  tp1 = mean((traj.vel(interval_common,2:3))>0);
%  tp2 = mean((traj.vel(interval_common,5:6))>0);