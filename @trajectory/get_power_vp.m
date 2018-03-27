%function [pow1,pow2,tp1,tp2] = get_power(traj,par)
function [pow1, pow2, tc1, tc2] = get_power_vp(traj,par)

% keyboard
% 
[md1,i1]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,2)).^2+...
                   (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(1,3)).^2));
[md2,i2]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(2,2)).^2+...
                   (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,3)).^2));
               

%   keyboard
               
tc1 = traj.time(traj.interval{par}(i1));
tc2 = traj.time(traj.interval{par}(i2));


% tc11 = tc1(1,1);
% tc22 = tc2(2,1);
% tc12 = tc1(2,1);
% tc21 = tc2(1,1);

%  keyboard
interval = max(traj.interval{1}(1),traj.interval{2}(1)):min(traj.interval{1}(end),traj.interval{2}(end));

[pow1, l1] = mean(sumabs(traj.vel(interval,2:3).*traj.intforce(interval,2:3))); %mean value on trial
[pow2, l2] = mean(sumabs(traj.vel(interval,5:6).*traj.intforce(interval,5:6)));

% tc1 = traj.time(traj.interval{par}(l1));
% tc2 = traj.time(traj.interval{par}(l2));

% pow_vp1 = mean(sumabs(traj.vel(tc11, 2:3).*traj.inforce(tc11, 2:3)))-mean(sumabs(traj.vel(tc21, 5:6).*traj.inforce(tc21, 5:6)));
% pow_vp2 = mean(sumabs(traj.vel(tc22, 5:6).*traj.inforce(tc22, 5:6)))-mean(sumabs(traj.vel(tc12, 2:3).*traj.inforce(tc12, 2:3)));
