function plot_power(traj,col)
% plots trajectory power in specified colour
% 
% (C) V. Sanguineti, 2008

cols = [0 0 1; 1 0 0];
cols1 = [.5 .5 0; 0 .5 .5];
cols2 = [0 .5 .5; .5 .5 0];

% find minimum reaction time and maximum termination time
if traj.nhands==2
    first = [traj.interval{1}(1)  traj.interval{2}(1)];  % these are the reaction times 
    last = [traj.interval{1}(end)  traj.interval{2}(end)];  % these are the termination times
    [fi,hf]=min(first);
    [la,hl]=max(last);
    interval = fi:la;
end
 %keyboard
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
          ind1 = (power(interval)>0);
          %ind2 = power(power(interval)<0);
          
            
        otherwise, 
            speed = sqrt(traj.vel(:,2*(ha-1)+1).^2+traj.vel(:,2*(ha-1)+2).^2);
            forcenorm = sqrt(traj.intforce(:,2*(ha-1)+1).^2+traj.intforce(:,2*(ha-1)+2).^2);
%           power = (traj.vel(:,2*(ha-1)+1).*traj.intforce(:,2*(ha-1)+1)+...
%                      traj.vel(:,2*(ha-1)+2).*traj.intforce(:,2*(ha-1)+2))./speed./forcenorm;
             power = traj.vel(:,2*(ha-1)+1).*traj.intforce(:,2*(ha-1)+1)+...
                      traj.vel(:,2*(ha-1)+2).*traj.intforce(:,2*(ha-1)+2);
                  
 
    end
    
    
%     keyboard
    if traj.nhands ==1
      line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),power(traj.interval{ha}),'col',col)
    else
      line(traj.time(interval)-traj.time(interval(1)),power(interval),'col',cols(ha,:))
%         line(traj.time(interval)-traj.time(interval(1)),sign(power(interval)),'col',cols(ha,:),'lines','-') %% sign power
%       line(traj.time(interval)-traj.time(interval(1)),ind1,'col',cols1(ha,:),'lines','-')
       
%     plot(traj.time(interval)-traj.time(interval(1)),ind1,'*')
       %plot(traj.time(interval(find(power<0))-traj.time(interval(1)), sign(power)));
      %line(traj.time(interval)-traj.time(interval(1)),ind2,'col',cols2(ha,:),'lines','-')
      
       
    end
end