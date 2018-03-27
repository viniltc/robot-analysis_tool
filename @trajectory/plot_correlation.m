function plot_correlation(traj,col)
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

for ha=1:traj.nhands
    switch traj.trsize
        case 3, % in 3D trajectories, display [YZ] projection
      speed1= sqrt(traj.vel(:,1).^2 + traj.vel(:,2).^2+ traj.vel(:,3).^2);
      speed2= sqrt(traj.vel(:,4).^2 + traj.vel(:,5).^2+ traj.vel(:,6).^2);
      [R,lags] = xcov(speed1(interval),speed2(interval));
     %R = R./(std(speed1(interval))*std(speed2(interval)));
       R = (speed1(interval)-speed2(interval)).^2;
       
       
    %[R , L] = xcorr(traj.vel(:,1),traj.vel(:,2))./(std(traj.vel(:,1))+std(traj.vel(:,2)));
   % R = xcorr(traj.vel(:,1),traj.vel(:,2))./(std(traj.vel(:,1))+std(traj.vel(:,2)));
%     
%     sample_frequency=1e3; %sample frequency in Hz
%     time_window=2e-3; %time window length in seconds
%     time_step=1e-3; %time step in seconds
%     max_lag_time=1e-3; %maximum lag time in seconds initially 1e-3
%     [lag_time,twin,xcl]=...
%     timewindow_xcorr(traj.vel(:,1),traj.vel(:,2),sample_frequency,...
%     time_window,time_step,max_lag_time);
% 
%     Rv = sum(traj.vel(:,1).*traj.vel(:,2))./sqrt(sum(traj.vel(:,1).^2).*traj.vel(:,2).^2)
%   
   
      otherwise, 
         %  R = xcorr(traj.vel(:,2*(ha-1)+2),traj.vel(:,2*(ha-1)+3));
         
        
    end
    if traj.nhands ==1
      % line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),speed(traj.interval{ha}),'col',col,'lines','-') 
      line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),R(traj.interval{ha}),'col',col,'lines','-')
      line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),R(traj.interval{ha}))
  
    else
 
    %line(traj.time(interval),corr(interval))
    %line(traj.time(interval)-traj.time(interval(1)),R(interval),'col',cols(ha,:),'lines','-')
    %line(traj.time(interval)-traj.time(interval(1)),Rv(interval))
    %plot(L/1000,R,'k')
    line(traj.time(interval)-traj.time(interval(1)),R)
    %line(lags/traj.fc,R)
    %line(L,R(interval))
    %line(traj.time(interval)-traj.time(interval(1)),twin(interval),'col',cols(ha,:),'lines','-'
    %line(twin(interval));
    %pcolor(twin,lag_time*1e3,xcl')
    %size(lag_time)
    %size(xcl')

    %plot(lag_time,mean(xcl'))

    end
end