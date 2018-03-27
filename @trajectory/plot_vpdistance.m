function plot_vpdistance(traj,col,par)
% plots trajectory speed in specified colour
% (C) V. Sanguineti, 2008

cols = [0 0 1; 1 0 0];

% [md1,i1]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,2)).^2+...
%                    (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(1,3)).^2));
% [md2,i2]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(2,2)).^2+...
%                    (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,3)).^2));


% find minimum reaction time and maximum termination time

if traj.nhands==2
    first = [traj.interval{1}(1)  traj.interval{2}(1)];  % these are the reaction times 
    last  = [traj.interval{1}(end)  traj.interval{2}(end)];  % these are the termination times
    [fi,hf]=min(first);
    [la,hl]=max(last);
    interval = fi:la;
    
    %%%%%%%%%%%%%%%%%%%%%%%
%     tc1 = traj.time(traj.interval{par}(i1));
%     tc2 = traj.time(traj.interval{par}(i2));
end

for vp=1:size(traj.viapoints,1)
    switch traj.trsize
        case 3, % in 3D trajectories, display [YZ] projection
            ha = 1;
            dist1_vp{vp} = sqrt((traj.pos(traj.interval{ha},3*(ha-1)+2)-traj.viapoints(vp,2)).^2+...
                                (traj.pos(traj.interval{ha},3*(ha-1)+3)-traj.viapoints(vp,3)).^2);
            ha=2;
            dist2_vp{vp} = sqrt((traj.pos(traj.interval{ha},3*(ha-1)+2)-traj.viapoints(vp,2)).^2+...
                                (traj.pos(traj.interval{ha},3*(ha-1)+3)-traj.viapoints(vp,3)).^2);
        otherwise, 
            ha=1;
            dist1_vp{vp}= sqrt((traj.pos(traj.interval{ha},2*(ha-1)+1)-traj.viapoints(vp,1)).^2+...
                       (traj.pos(traj.interval{ha},2*(ha-1)+2)-traj.viapoints(vp,2)).^2);
            ha=2;
            dist2_vp{vp} = sqrt((traj.pos(traj.interval{ha},2*(ha-1)+1)-traj.viapoints(vp,1)).^2+...
                       (traj.pos(traj.interval{ha},2*(ha-1)+2)-traj.viapoints(vp,2)).^2);
       
    end
end
    if traj.nhands ==1
           % line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),speed(traj.interval{ha}),'col',col,'lines','-')
    else
            ha=1;
            
            line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),dist1_vp{1},'col',cols(1,:),'lines','-')   
            line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),dist1_vp{2},'col',cols(1,:),'lines','--')  

            
%            [md1,i1]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,2)).^2+...
%                      (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(1,3)).^2));
%            [md2,i2]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(2,2)).^2+...
%                    (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,3)).^2));
            
            ha=2;
            line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),dist2_vp{2},'col',cols(2,:),'lines','-')  
            line(traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1)),dist2_vp{1},'col',cols(2,:),'lines','--') 
         %   line(traj.time(traj.interval{par}(1))*[1 1],[0 0.1],'col',[.8 .8 .8],'lines','--')
         %   line(traj.time(traj.interval{par}(2))*[1 1],[0 0.1],'col',[.8 .8 .8],'lines','--')
         
   
    end
    



