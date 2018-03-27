function plot_traj(traj,col)
% plots trajectory in specified colour
% 
% (C) V. Sanguineti, 2008,2015


for ha=1:traj.nhands
    tmp = rgb2hsv(col);
    tmp(2) = tmp(2)-(ha-1)*0.5; % reduce saturation 
    %col = hsv2rgb(tmp);
    col = [0 0 0.7; 0.7 0 0]; % P1 is blue, P2 is red
    col = [0 0.4 0.9; 0.9 0 0.4]; % P1 is blue, P2 is red
    
            
    switch traj.trsize
        case 2, % planar robot
            if traj.xyT(1,1)==0 && traj.xyT(1,2)==0
                line(traj.pos(traj.interval{ha},2*(ha-1)+1)-traj.pos(traj.interval{ha},2*(ha-1)+1),...
                     traj.pos(traj.interval{ha},2*(ha-1)+2)-traj.pos(traj.interval{ha},2*(ha-1)+2),'col',col,'lines','none','marker','.','markersize',3)
            else
                line(traj.pos(traj.interval{ha},2*(ha-1)+1),traj.pos(traj.interval{ha},2*(ha-1)+2),'col',col,'lines','none','marker','.','markersize',3)
            end
            
        case 3, % 3D robot... in this case (dyad), show in [YZ] space
            if traj.xyT(1,1)==0 && traj.xyT(1,2)==0 && traj.xyT(1,3)==0
                line(traj.pos(traj.interval{ha},3*(ha-1)+2),...
                     traj.pos(traj.interval{ha},3*(ha-1)+3),...
                    'col',col(ha,:),'lines','none','marker','.','markersize',3)
%                 line(traj.pos(:,3*(ha-1)+2)-traj.pos(:,3*(ha-1)+2),...
%                      traj.pos(:,3*(ha-1)+3)-traj.pos(:,3*(ha-1)+3),...
%                      'col',col,'lines','none','marker','.','markersize',3)
            else
                line(traj.pos(traj.interval{ha},3*(ha-1)+2),...
                     traj.pos(traj.interval{ha},3*(ha-1)+3),...
                     'col',col(ha,:),'lines','none','marker','.','markersize',3)
                %line(traj.pos(:,3*(ha-1)+2),traj.pos(:,3*(ha-1)+3),...
                %     'col',col(ha,:),'lines','none','marker','.','markersize',3)
            end
    end
end