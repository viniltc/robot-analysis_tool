function plot_intforce(tset,phaseno,mode)
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

tlim = [-0.1 3];
slim = [0 14];
munit = 'N';
scalebar.center=[1 0.05];
scalebar.xdim=0.20;
scalebar.ydim=0.10;
spdstrval=['XY Int_force - ', num2str(phaseno)];

th = atan2(tset.targets(:,2),tset.targets(:,1));

figure
set(gcf,'Name',spdstrval)

set(gcf,'pos',[0+(phaseno-1)*100 250 300 300]);
if length(th)>1
    for i = 1:length(th)
            hva(i)=axes('position',[0.4+0.3*cos(th(i)),0.4+0.3*sin(th(i)),0.3,0.3]);
            %set(gca,'xlim',tlim,'ylim',slim);
            set(gca,'xlim',tlim);
            xlabel(['time [s]']);
            ylabel(['Interaction force [', munit,']'],'FontSize',15);
            set(gca,'vis','off');
    end
else
    hva=gca;
    set(gca,'xlim',tlim,'ylim',slim);
    %set(gca,'xlim',tlim);
    xlabel(['time [s]'], 'FontSize',15);
    ylabel(['Interaction force [', munit,']'], 'FontSize',15);
    %set(gca,'vis','off');
end    

%line(scalebar.center(1)+[0 scalebar.xdim],scalebar.center(2)+[0 0],'col','k');
%line(scalebar.center(1)+[0 0],scalebar.center(2)+[0 scalebar.ydim],'col','k');
%set(gca,'visible','off');

thr = 0:pi/20:2*pi;

cnt_nul = 0;
cnt_for = 0;
for trial = 1:tset.Ntrials
for i = 1:length(th)
            axes(hva(i));
               
           if mean(isforce(tset.traj{trial}))>0.5 % force is on
               plot_int_force(tset.traj{trial},col_for);
           else
               plot_int_force(tset.traj{trial},col_nul);
           end
end   
end