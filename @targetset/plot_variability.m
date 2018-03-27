function plot_variability(tset,phaseno,mode)
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

cols = ['b' 'r'];
%tlim = [-0.1 3];
tlim = [0 1];
slim = [0 0.7];

munit = 'm';
scalebar.center=[1 0.05];
scalebar.xdim=0.20;
scalebar.ydim=0.10;
spdstrval=['Variability - ', num2str(phaseno)];

th = atan2(tset.targets(:,2),tset.targets(:,1));

figure
set(gcf,'Name',spdstrval)
set(gcf,'pos',[0+(phaseno-1)*100 250 300 300]);

if length(th)>1
    for i = 1:length(th)
            hva(i)=axes('position',[0.4+0.3*cos(th(i)),0.4+0.3*sin(th(i)),0.3,0.3]);
            %set(gca,'xlim',tlim,'ylim',slim);
            %set(gca,'xlim',tlim);
            xlabel(['time [s]']);
           % ylabel(['variability [', munit,']']);
            set(gca,'vis','off');
    end
else
    hva=gca;
    %set(gca,'xlim',tlim,'ylim',slim);
    set(gca,'xlim',tlim);
    xlabel(['time [s]']);
    ylabel(['variability [', munit,']']);
    %set(gca,'vis','off');
end

cnt_nul = 0;
cnt_for = 0;
for trial = 1:tset.Ntrials % for each trial...
    for ha=1:2
        for i = 1:length(th) % for each target...
            axes(hva(i));
            
            ns=100;
           [ntime,npos] = get_npos( tset.traj{trial},ns );
                
            if mean(isforce(tset.traj{trial}))>0.5 % force is on
                normpos_for{ha}(trial,:,:)=npos{ha};
            else
                normpos_for{ha}(trial,:,:)=npos{ha};
            end
            
        end
    end
    
 end

for t=1:ns
        [V1,D1]=eig(cov(reshape(normpos_for{1}(:,t,2:3),tset.Ntrials,2)));
        [V2,D2]=eig(cov(reshape(normpos_for{2}(:,t,2:3),tset.Ntrials,2)));
        
       variability{1}(t) = sqrt(sqrt(det(D1)));
       variability{2}(t) = sqrt(sqrt(det(D2)));
    end
    line(ntime,variability{1},'col',cols(1),'marker','*');
    line(ntime,variability{2},'col',cols(2),'marker','*');
    
ylim([0 2e-2])