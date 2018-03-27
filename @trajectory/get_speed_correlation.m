%function [speed1,speed2,lag_time] = get_speed_correlation(traj,par)
function [speed1,speed2,R] = get_speed_correlation(traj,par)
  

    first=[traj.interval{1}(1)  traj.interval{2}(1)];  % these are the reaction times 
    last =[traj.interval{1}(end)  traj.interval{2}(end)];  % these are the termination times
    [fi,hf]=min(first);
    [la,hl]=max(last);
    interval = fi:la;
   
    speed1 = sqrt(sum(traj.vel(interval,2:3).^2,2));
    speed2 = sqrt(sum(traj.vel(interval,5:6).^2,2));



% speed1 = sqrt(sum(traj.vel(traj.interval{1},2:3).^2));
% speed2 = sqrt(sum(traj.vel(traj.interval{2},5:6).^2));
%R= mean(xcorr2(speed1,speed2)./(std(speed1)+std(speed2)));
%R= mean(xcorr2(speed1,speed2));
rr = corrcoef(speed1,speed2);
R= rr(1,2).^2;



