function plot_1trajs(tset,phaseno)
bwfigs = 0;
if bwfigs
    col_tgt = [0 0 0];
    col_for = [0 0 0];
    col_nul = [0 0 0];
    col_rot = [0 0 0];
    col_frot = [0 0 0]; 
else
    col_tgt=  [0 1 0];
    col_frot = [0 1 1]; %ciano
    col_nul = [0 0 0]; %rosso
    col_rot = [1 0 0]; %rosso
    col_for = [0 0 1]; %magenta
end


 rotations = get_rotations(tset);
 forces = get_forces(tset);
     rotation_trial = not(sum(rotations+forces)==36);%any(rotations)%
     force_trial=any(forces);
     
     %calcolo learning index (smith and shadmehr, 2005)
     
 %    keyboard;
for trial = 1:tset.Ntrials
     isf = mean(isforce(tset.traj{trial}))>0.5;
     isr = mean(isrotation(tset.traj{trial}))>0.5;

    trial;
     isfor = mean(isforce(tset.traj{trial}));
     isrot = mean(isrotation(tset.traj{trial}));

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
    
    target_angles = atan2(tset.targets(:,2), tset.targets(:,1))';
    tgtraj{trial}=[];
    
     [pos,notarget,tgt]=get_onetraj(tset.traj{trial},tset.targets,target_angles);   
    switch notarget
        case 1
            tgtraj{trial}=pos;
        case 2
            tgtraj{trial}=pos;
        case 3
            tgtraj{trial}=pos; 
        case 4
            tgtraj{trial}=pos; 
        case 5
            tgtraj{trial}=pos; 
        case 6
            tgtraj{trial}=pos;  
    end

    numtgt(trial)=notarget;
    
   
end

figure
 set(gcf,'pos',[10+10*phaseno 250 200 200]);
for i=0:length(target_angles)-1
    ind=find(numtgt==i+1)
%     figure(length(target_angles)*phaseno+i)
%     keyboard
    subplot(2,3,i+1)
   
    axis('equal');
    the_title=['Target ' num2str(i+1) ' phase ' num2str(phaseno) ];
    line(tgt(1),tgt(2),'col',col_tgt,'Marker','.','MarkerSize',20,'MarkerFaceColor','auto');
    title(the_title)
    for j=1:length(ind)
        xmax(j)=max(tgtraj{ind(j)}(1,:));
%         ymax(j)=max(tgtraj{ind(j)}(:,2))
%         xmin(j)=min(tgtraj{ind(j)}(:,1))
%         ymin(j)=min(tgtraj{ind(j)}(:,2))
        
        line(tgtraj{ind(j)}(1,:),tgtraj{ind(j)}(2,:),'col',cols);
%         xlim([0 xmax(j)]);
    end
    
end