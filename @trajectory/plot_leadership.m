function [pow1, pow2, tc11, tc22, tc12, tc21] = plot_leadership(traj,col)
% plots trajectory power in specified colour
% 
% (C) V. Sanguineti, 2008

col = [1 0 0; 0 0 1]; % P1 is red, P2 is blue
% find minimum reaction time and maximum termination time
if traj.nhands==2
    first = [traj.interval{1}(1)  traj.interval{2}(1)];  % these are the reaction times 
    last = [traj.interval{1}(end)  traj.interval{2}(end)];  % these are the termination times
    [fi,hf]=min(first);
    [la,hl]=max(last);
    interval = fi:la;
end
pp = [];
for ha=1:traj.nhands
    switch traj.trsize
        case 3, % in 3D trajectories, display [YZ] projection
            speed = sqrt(traj.vel(:,3*(ha-1)+2).^2+traj.vel(:,3*(ha-1)+3).^2);
            forcenorm = sqrt(traj.intforce(:,3*(ha-1)+2).^2+traj.intforce(:,3*(ha-1)+3).^2);
%           power = (traj.vel(:,3*(ha-1)+2).*traj.intforce(:,3*(ha-1)+2)+...
%                      traj.vel(:,3*(ha-1)+3).*traj.intforce(:,3*(ha-1)+3))./speed./forcenorm;
            power = traj.vel(:,3*(ha-1)+2).*traj.intforce(:,3*(ha-1)+2)+...
                     traj.vel(:,3*(ha-1)+3).*traj.intforce(:,3*(ha-1)+3);
%            power = dot(traj.vel(:,3*(ha-1)+2),traj.intforce(:,3*(ha-1)+2))+... 
%                      dot(traj.vel(:,3*(ha-1)+3),traj.intforce(:,3*(ha-1)+3));
     
                 
        otherwise, 
            speed = sqrt(traj.vel(:,2*(ha-1)+1).^2+traj.vel(:,2*(ha-1)+2).^2);
            forcenorm = sqrt(traj.intforce(:,2*(ha-1)+1).^2+traj.intforce(:,2*(ha-1)+2).^2);
%           power = (traj.vel(:,2*(ha-1)+1).*traj.intforce(:,2*(ha-1)+1)+...
%                      traj.vel(:,2*(ha-1)+2).*traj.intforce(:,2*(ha-1)+2))./speed./forcenorm;
            power = traj.vel(:,2*(ha-1)+1).*traj.intforce(:,2*(ha-1)+1)+...
                      traj.vel(:,2*(ha-1)+2).*traj.intforce(:,2*(ha-1)+2);
                  
%            power = dot(traj.vel(:,2*(ha-1)+1),traj.intforce(:,2*(ha-1)+1))+... 
%                      dot(traj.vel(:,2*(ha-1)+2),traj.intforce(:,2*(ha-1)+2))
                  
    end
    
    pow(:,ha)=power;
    
    [md1,i1]= min(sqrt((traj.pos(traj.interval{ha},3*(ha-1)+2)-traj.viapoints(1,2)).^2+...
                   (traj.pos(traj.interval{ha},3*(ha-1)+3)-traj.viapoints(1,3)).^2));
    [md2,i2]= min(sqrt((traj.pos(traj.interval{ha},3*(ha-1)+2)-traj.viapoints(2,2)).^2+...
                       (traj.pos(traj.interval{ha},3*(ha-1)+3)-traj.viapoints(2,3)).^2));

    tc1_1(ha,1) = traj.time(traj.interval{ha}(i1));  % TC wrt target onset
    tc1_2(ha,1) = traj.time(traj.interval{ha}(i2));

end

% get intersection of two intervals
interval = min(traj.interval{1}(1),traj.interval{2}(1)):max(traj.interval{1}(end),traj.interval{2}(end));    

% keyboard

% This is TC wrt movement start
for ha=1:traj.nhands
    tc1_1(ha,1) = tc1_1(ha,1)-traj.time(interval(1));
    tc1_2(ha,1) = tc1_2(ha,1)-traj.time(interval(1));
end
% 
%   keyboard

pow1=pow(interval,1);
pow2=pow(interval,2);
% pow1=pow(traj.interval{1},1);
% pow2=pow(traj.interval{2},2);
tc11 = tc1_1(1,1);
tc22 = tc1_2(2,1);
tc12 = tc1_1(2,1);
tc21 = tc1_2(1,1);

% pow11 = pow(tc11,1);
% pow12 = pow(tc12,1);
% pow21 = pow(tc21,2);
% pow22 = pow(tc22,2);
% pow_vp1 = pow(tc11,1)-pow(tc21,2);
% pow_vp2 = pow(tc22,2)-pow(tc12,1);
end



