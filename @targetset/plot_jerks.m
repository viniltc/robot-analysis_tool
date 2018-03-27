function plot_jerks(tset,phaseno)
bwfigs = 0;
if bwfigs
    col_avg='k';
    col_fig='k';
    col_tgt='k';
    col_for = 'k';
    col_nul = 'k';
else
    col_avg='b';
    col_fig='r';
    col_tgt='g';
    col_for = 'm';
    col_nul = 'k';
end

tlim = [-0.1 1];
slim = [0 60];
munit = 'm/s^2';
scalebar.center=[1 0.05];
scalebar.xdim=0.20;
scalebar.ydim=10;
spdstrval=['XY Jerk - ', num2str(phaseno)];

th = atan2(tset.targets(:,2),tset.targets(:,1));

figure
set(gcf,'Name',spdstrval)
set(gcf,'pos',[350 250 350 350]);

for i = 1:length(th)
            hva(i)=axes('position',[0.4+0.3*cos(th(i)),0.4+0.3*sin(th(i)),0.3,0.3]);
set(gca,'xlim',tlim,'ylim',slim);
xlabel(['TIME [s]']);
ylabel(['JERK [', munit,']']);
            %set(gca,'vis','off');
end


thr = 0:pi/20:2*pi;

cnt_nul = 0;
cnt_for = 0;
for trial = 1:tset.Ntrials
    for i = 1:length(th)
            axes(hva(i));
               

    if mean(isforce(tset.traj{trial}))>0.5 % force is on
     plot_jerk(tset.traj{trial},col_for);
%      cnt_for = cnt_for+1;
%      if(cnt_for>(2*tset.Ntargets)) break
%      end
    else
     plot_jerk(tset.traj{trial},col_nul);
%      cnt_nul = cnt_nul+1;
%      if(cnt_nul>(2*tset.Ntargets)) break
%      end
    end

   end
end     