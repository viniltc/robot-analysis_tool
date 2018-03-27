function plot_trajs(tset,phaseno,varargin)
bwfigs = 0;
if bwfigs
    col_tgt = [0 0 0];
    col_for = [0 0 0];
    col_nul = [0 0 0];
    col_rot = [0 0 0];
    col_frot = [0 0 0]; 
    col_sta = [0 0 0];
else
    col_tgt=  [0 1 0];
    col_tgt=  [1 1 0]; % yellow
    col_frot = [0 1 1]; %ciano
    col_nul = [0 0 0]; %rosso
    col_rot = [1 0 0]; %rosso
    col_for = [0 0 1]; %magenta
    col_sta = [1 1 1];
end

%null=nero=[0 0 0]
%rotazione=rosso=[1 0 0]
%forza=blu=[0 0 1]
%=[1 0 1]=magenta forza + rotazione 

switch tset.type
    case 'dyad',
        lx=[-0.02 0.12];
%         lx=[-0.03 0.25];
%         ly= [-0.01 0.01]; % zoom image
        ly= [-0.065 0.065];
        munit = 'm';
        scalebar.center=[0.10 -0.06];
        scalebar.xdim=0.01;
        scalebar.ydim=0.01;
        trajstrval=['YZ Trajectory - ', num2str(phaseno)];
    otherwise,
        lx = [-0.160 0.160];
        ly = lx;
        munit = 'm';
        scalebar.center=[0.120 -0.120];
        scalebar.xdim=0.020;
        scalebar.ydim=0.020;
        trajstrval=['XY Trajectory - ', num2str(phaseno)];
end

figure
set(gcf,'Name',trajstrval)
% set(gcf,'pos',[0 250 300 300]);
set(gcf,'pos',[0+(phaseno-1)*100 250 300 300]);
axis('equal');
set(gca,'xlim',lx,'ylim',ly);
    
xlabel(['X [', munit,']'], 'FontSize',15);
ylabel(['Y [', munit,']'], 'FontSize',15);
    
% line(scalebar.center(1)+[0 scalebar.xdim],scalebar.center(2)+[0 0],'col','k');  % xy scale line in bottom
% line(scalebar.center(1)+[0 0],scalebar.center(2)+[0 scalebar.ydim],'col','k');
set(gca,'visible','off');

thr = 0:pi/20:2*pi;
for i = 1:tset.Ntargets
       %line(tset.targets(i,1)+tset.tgsize.*cos(thr),...
       %     tset.targets(i,2)+tset.tgsize.*sin(thr),'col',col_tgt);
       patch(tset.targets(i,1)+tset.tgsize.*cos(thr),...
             tset.targets(i,2)+tset.tgsize.*sin(thr),col_tgt,'edgecol','k');
end
patch(tset.tgsize.*cos(thr),tset.tgsize.*sin(thr),col_sta,'edgecol','k');




 rotations = get_rotations(tset);
 forces = get_forces(tset);
     rotation_trial = not(sum(rotations+forces)==36);%any(rotations)%
     force_trial=any(forces);
 %    keyboard;
for trial = 1:tset.Ntrials
     isf = mean(isforce(tset.traj{trial}))>0.5;
     isr = mean(isrotation(tset.traj{trial}))>0.5;

    trial
     isfor = mean(isforce(tset.traj{trial}))
     isrot = mean(isrotation(tset.traj{trial}))

   if force_trial & rotation_trial
        if (forces(trial) & rotations(trial))
             cols = col_frot;
        elseif forces(trial) & ~rotations(trial)  %%ct sulla rot!
             cols = col_for;
        elseif ~forces(trial) & rotations(trial) %%ct sulla for!
             cols = col_rot; 
        else
             cols = col_nul;
        end
    end
    if ~force_trial & rotation_trial %TUTTE TRANNE BLOCCO 3 DALLA 6 ALLA 13
        if rotations(trial)
             cols = col_rot; 
        else
             cols = col_nul; 
        end  
    end
   
    if ~force_trial & ~rotation_trial  %
             cols = col_nul; 
    end
% % %     if isf && ~isr
% % %         cols = col_for;    
% % %     elseif isr && ~isf 
% % %         cols = col_rot; %null
% % %     elseif isf && isr
% % %         cols = col_frot;
% % %     else
% % %         cols = col_nul;
% % %     end
%     target_angles = atan2(tset.targets(:,2), tset.targets(:,1))';
%     plot_onetraj(tset.traj{trial},cols,target_angles)
    %if tset.targetno(trial)==2
    
    
       plot_traj(tset.traj{trial},cols);
    %end
    
    
    
 if nargin==3
     feval(varargin{1})
 end
    
    
    
   
end