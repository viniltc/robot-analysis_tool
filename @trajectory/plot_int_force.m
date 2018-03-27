function plot_int_force(traj,col)
% plots trajectory interaction force in specified colour
% 
% (C) V. Sanguineti, 2008,2016

cols = [1 0 0; 0 0 1];

% find minimum reaction time and maximum termination time
if traj.nhands==2
    first = [traj.interval{1}(1)  traj.interval{2}(1)];  % these are the reaction times 
    last = [traj.interval{1}(end)  traj.interval{2}(end)];  % these are the termination times
    [fi,hf]=min(first);
    [la,hl]=max(last);
    interval = fi:la;
end


% keyboard
for ha=1:1
    switch traj.trsize
        case 3, % in 3D trajectories, display [YZ] projection
           
           forcenorm = sqrt(traj.intforce(:,3*(ha-1)+2).^2+traj.intforce(:,3*(ha-1)+3).^2);          
        otherwise, 
           forcenorm = sqrt(traj.intforce(:,2*(ha-1)+1).^2+traj.intforce(:,2*(ha-1)+2).^2);
    end
    if traj.nhands ==1
%       line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),forcenorm(traj.interval{ha}),'col',col,'lines','-')
      line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),forcenorm(traj.interval{1}),'col','k','lines','-')
    else
%        line(traj.time(interval)-traj.time(interval(1)),forcenorm(interval),'col',cols(ha,:),'lines','-')
       line(traj.time(interval)-traj.time(interval(1)),forcenorm(interval),'col','k','lines','-')
    end
end