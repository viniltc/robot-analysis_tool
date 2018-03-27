function plot_speed_correlation(traj,col)
% plots trajectory speed in specified colour
% 
% (C) V. Sanguineti, 2008

cols = [1 0 0; 0 0 1];

% find minimum reaction time and maximum termination time
if traj.nhands==2
    first = [traj.interval{1}(1)  traj.interval{2}(1)];  % these are the reaction times 
    last = [traj.interval{1}(end)  traj.interval{2}(end)];  % these are the termination times
    [fi,hf]=min(first);
    [la,hl]=max(last);
    interval = fi:la;
end

correlation = xcorr(traj.vel(:,3*.....

for ha=1:traj.nhands
    switch traj.trsize
        case 3, % in 3D trajectories, display [YZ] projection
            speed = sqrt(traj.vel(:,3*(ha-1)+2).^2+traj.vel(:,3*(ha-1)+3).^2);  
            %[a,b]= speed;
            correlation = xcorr((traj.vel(:,3*(ha-1)+2),(traj.vel(:,3*(ha-1)+3)))
        otherwise, 
            speed = sqrt(traj.vel(:,2*(ha-1)+1).^2+traj.vel(:,2*(ha-1)+2).^2);
    end
    if traj.nhands ==1
      line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),speed(traj.interval{ha}),'col',col,'lines','-')
    else
       line(traj.time(interval)-traj.time(interval(1)),speed(interval),'col',cols(ha,:),'lines','-')
    end
end

