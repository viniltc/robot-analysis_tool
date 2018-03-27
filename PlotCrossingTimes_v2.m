%% Plot via crossing time, tc11, tc12, tc21, tc22 for the two groups

clear all
close all

 % Set location of data
datadir = '../../Data/';
% Details on where are results
resdir = [datadir, 'results/'];
subH = 1;
subVH = 2;  % 2 removed

load([resdir, 'tc11_statmat.txt']); % S1 crossing VP1
load([resdir, 'tc12_statmat.txt']); % S1 crossing VP2
load([resdir, 'tc21_statmat.txt']); % S2 crossing VP1
load([resdir, 'tc22_statmat.txt']); % S2 crossing VP2

load([resdir, 'md12_statmat.txt']); % S1 crossing VP2
load([resdir, 'md21_statmat.txt']); % S2 crossing VP1

load([resdir, 'tc11_statmat_nofor.txt']);
load([resdir, 'tc12_statmat_nofor.txt']);
load([resdir, 'tc21_statmat_nofor.txt']);
load([resdir, 'tc22_statmat_nofor.txt']);

load([resdir, 'ts1_statmat.txt']);
load([resdir, 'ts2_statmat.txt']);
load([resdir, 'te1_statmat.txt']);
load([resdir, 'te2_statmat.txt']);

group = tc11_statmat(:,1);
tc11 = tc11_statmat(:,2:end);
tc12 = tc12_statmat(:,2:end);
tc21 = tc21_statmat(:,2:end);
tc22 = tc22_statmat(:,2:end);

md12 = md12_statmat(:,2:end);
md21 = md21_statmat(:,2:end);

%  notok21 = find(md21>0.02); % if md is less than 3 cm there is no crossing at all
%  notok12 = find(md12>0.02); % if md is less than 3 cm there is no crossing at all
%  tc21(notok21) = nan;
%  tc12(notok12) = nan;


ts1 = ts1_statmat(:,2:end);
ts2 = ts2_statmat(:,2:end);
te1 = te1_statmat(:,2:end);
te2 = te2_statmat(:,2:end);

% load per-trial data
load ([resdir,'tc11.mat']);
tc11_ts = (ind_ts);
load ([resdir,'tc12.mat']);
tc12_ts = (ind_ts);
load ([resdir,'tc21.mat']);
tc21_ts = (ind_ts);
load ([resdir,'tc22.mat']);
tc22_ts = (ind_ts);

% This is to assess goodness of crossing...
load ([resdir,'md12.mat']);
md12_ts = (ind_ts);
load ([resdir,'md21.mat']);
md21_ts = (ind_ts);

mdavg_ts = 0.5*(md12_ts+md21_ts)


%  notok21 = find(md21_ts>0.02); % if md is less than 3 cm there is no crossing at all
%  notok12 = find(md12_ts>0.02); % if md is less than 3 cm there is no crossing at all
%  tc21_ts(notok21) = nan;
%  tc12_ts(notok12) = nan;
    
% Load per-trial data
load([resdir, 'forces.mat']);

load ([resdir,'ts1.mat']);
ts1_ts = ind_ts;
load ([resdir,'ts2.mat']);
ts2_ts = ind_ts;
load ([resdir,'te1.mat']);
te1_ts = ind_ts;
load ([resdir,'te2.mat']);
te2_ts = ind_ts;



subjs=size(tc11,1);
tsets=size(tc22,2);
plotlabel = 'Via-point crossing time [s]';
yrange = [0 4];

allvals11 = [];
allvals22 = [];
allvals12 = [];
allvals21 = [];
rtime = 0*ts1;
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 250 350])
    patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel(plotlabel)
    title(sprintf('S%d',subj))

    % reaction time and movement time
    rtime(subj,:) = min([ts1(subj,:); ts2(subj,:)]);
    mtime(subj,:) = max([te1(subj,:)-rtime(subj,:); te2(subj,:)-rtime(subj,:)]);
    
    % these denote collaboration level...
    deltatc1(subj,:) = tc21(subj,:)-tc11(subj,:);
    deltatc2(subj,:) = tc12(subj,:)-tc22(subj,:);
    
    %patch([1:tsets tsets:-1:1],[ts2(subj,:) fliplr(te2(subj,:))],'b');
    %patch([1:tsets tsets:-1:1],[ts1(subj,:) fliplr(te1(subj,:))],'r');
   w1= line(1:tsets,(tc11(subj,:)-rtime(subj,:)),'col','b');
   w2= line(1:tsets,(tc21(subj,:)-rtime(subj,:)),'col','r','lines',':');
   w3= line(1:tsets,(tc12(subj,:)-rtime(subj,:)),'col','b','lines',':');
   w4= line(1:tsets,(tc22(subj,:)-rtime(subj,:)),'col','r');
%       line(1:tsets,tc12(subj,:),'col','k','lines','-.');
%     line(1:tsets,(tc11(subj,:)-rtime(subj,:))./mtime(subj,:),'col','b','lines',':');
%     line(1:tsets,(tc21(subj,:)-rtime(subj,:))./mtime(subj,:),'col','r');
%     line(1:tsets,(tc12(subj,:)-rtime(subj,:))./mtime(subj,:),'col','b');
%     line(1:tsets,(tc22(subj,:)-rtime(subj,:))./mtime(subj,:),'col','r','lines',':');
%     

%line(1:tsets,mtime(subj,:),'col','r','lines',':');
    
    legend([w1,w2,w3,w4], {'S1-VP1','S2-VP1','S1-VP2','S2-VP2'})
    
    allvals11 = [allvals11; (tc11(subj,:)-rtime(subj,:))];
    allvals22 = [allvals22; (tc22(subj,:)-rtime(subj,:))];
    allvals12 = [allvals12; (tc12(subj,:)-rtime(subj,:))];
    allvals21 = [allvals21; (tc21(subj,:)-rtime(subj,:))];
    
    
    xlim([1 tsets])
    ylim(yrange)
    
    %patch([1:tsets tsets:-1:1],[ts2(subj,:) fliplr(te2(subj,:))],'b');
end
    
%% Plot average, per-epoch 

% group VH
figure
set(gcf,'pos',[100 100 250 350])
patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group VH')
hold on

e1=errorbar(1:tsets, mean(allvals11(subVH,:)), std(allvals11(subVH,:))./sqrt(length(subVH)),'bs:')
e2=errorbar(1:tsets, mean(allvals21(subVH,:)), std(allvals21(subVH,:))./sqrt(length(subVH)),'rs-')
hold on
e5=errorbar(1:tsets, mean(allvals12(subVH,:)), std(allvals12(subVH,:))./sqrt(length(subVH)),'bs-')
e6=errorbar(1:tsets, mean(allvals22(subVH,:)),std(allvals22(subVH,:))./sqrt(length(subVH)),'rs:')
e7=errorbar(1:tsets,mean(mtime(subVH,:)),...
                    std(mtime(subVH,:))./sqrt(length(subVH)),'col','k','lines',':');
l1=legend([e1,e2,e5, e6],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')
set(l1,'Location','north')
legend boxoff
xlim([1 tsets])
ylim(yrange)
% catch trials:
%errorbar(1:tsets,nanmean(tc21_statmat_nofor(subVH,2:end)),...
%    nanstd(tc21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'r:')
%errorbar(1:tsets,nanmean(tc12_statmat_nofor(subVH,2:end)),...
%    nanstd(tc12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'b:')

% group H 
figure
set(gcf,'pos',[100 100 250 350])
patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group H')
hold on

e1=errorbar(1:tsets, mean(allvals11(subH,:)), std(allvals11(subH,:))./sqrt(length(subH)),'bs:')
e2=errorbar(1:tsets, mean(allvals21(subH,:)), std(allvals21(subH,:))./sqrt(length(subH)),'rs-')
hold on
e5=errorbar(1:tsets, mean(allvals12(subH,:)), std(allvals12(subH,:))./sqrt(length(subH)),'bs-')
e6=errorbar(1:tsets, mean(allvals22(subH,:)), std(allvals22(subH,:))./sqrt(length(subH)),'rs:')
e7=errorbar(1:tsets,mean(mtime(subH,:)),...
    std(mtime(subH,:))./sqrt(length(subH)),'col','k','lines',':');
l1=legend([e1,e2,e5, e6],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')

%set(l1,'Location','north')
legend boxoff
xlim([1 tsets])
ylim(yrange)
% catch trials:
%errorbar(1:tsets,nanmean(tc21_statmat_nofor(subH,2:end)),...
%    nanstd(tc21_statmat_nofor(subH,2:end))./sqrt(length(subH)),'r:')
%errorbar(1:tsets,nanmean(tc12_statmat_nofor(subH,2:end)),...
%    nanstd(tc12_statmat_nofor(subH,2:end))./sqrt(length(subH)),'b:')

%% Plot average, per-epoch , relative to MT

% group VH
figure
set(gcf,'pos',[100 100 250 350])
yrrange = [0 3]
patch([1.5 11.5 11.5 1.5],yrrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group VH')
hold on

e1=errorbar(1:tsets, mean(allvals11(subVH,:)./mtime(subVH,:)), std(allvals11(subVH,:)./mtime(subVH,:))./sqrt(length(subVH)),'bs:')
e2=errorbar(1:tsets, mean(allvals21(subVH,:)./mtime(subVH,:)), std(allvals21(subVH,:)./mtime(subVH,:))./sqrt(length(subVH)),'rs-')
hold on
e5=errorbar(1:tsets, mean(allvals12(subVH,:)./mtime(subVH,:)), std(allvals12(subVH,:)./mtime(subVH,:))./sqrt(length(subVH)),'bs-')
e6=errorbar(1:tsets, mean(allvals22(subVH,:)./mtime(subVH,:)),std(allvals22(subVH,:)./mtime(subVH,:))./sqrt(length(subVH)),'rs:')
l1=legend([e1,e2,e5, e6],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')
set(l1,'Location','north')
legend boxoff
xlim([1 tsets])
ylim(yrrange)
% catch trials:
%errorbar(1:tsets,nanmean(tc21_statmat_nofor(subVH,2:end)),...
%    nanstd(tc21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'r:')
%errorbar(1:tsets,nanmean(tc12_statmat_nofor(subVH,2:end)),...
%    nanstd(tc12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'b:')

% group H 
figure
set(gcf,'pos',[100 100 250 350])
patch([1.5 11.5 11.5 1.5],yrrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group H')
hold on

e1=errorbar(1:tsets, mean(allvals11(subH,:)./mtime(subH,:)), std(allvals11(subH,:)./mtime(subH,:))./sqrt(length(subH)),'bs:')
e2=errorbar(1:tsets, mean(allvals21(subH,:)./mtime(subH,:)), std(allvals21(subH,:)./mtime(subH,:))./sqrt(length(subH)),'rs-')
hold on
e5=errorbar(1:tsets, mean(allvals12(subH,:)./mtime(subH,:)), std(allvals12(subH,:)./mtime(subH,:))./sqrt(length(subH)),'bs-')
e6=errorbar(1:tsets, mean(allvals22(subH,:)./mtime(subH,:)), std(allvals22(subH,:)./mtime(subH,:))./sqrt(length(subH)),'rs:')
%e7=errorbar(1:tsets,mean(mtime(subH,:)),...
%    std(mtime(subH,:))./sqrt(length(subH)),'col','k','lines',':');
l1=legend([e1,e2,e5, e6],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')

%set(l1,'Location','north')
legend boxoff
xlim([1 tsets])
ylim(yrrange)
% catch trials:
%errorbar(1:tsets,nanmean(tc21_statmat_nofor(subH,2:end)),...
%    nanstd(tc21_statmat_nofor(subH,2:end))./sqrt(length(subH)),'r:')
%errorbar(1:tsets,nanmean(tc12_statmat_nofor(subH,2:end)),...
%    nanstd(tc12_statmat_nofor(subH,2:end))./sqrt(length(subH)),'b:')


%% Plot average, per  null, trng, washout , relative to MT
figure
set(gcf,'pos',[100 100 300 400])
title(plotlabel)
errorbar(1+[-0.15 0.15], mean(allvals11(subVH,[1 10 ])./mtime(subVH,[1 10 ])), ...
                         std(allvals11(subVH,[1 10 ])./mtime(subVH,[1 10 ]))./sqrt(length(subVH)),'k.')
hold on
errorbar(2+[-0.15 0.15], mean(allvals22(subVH,[1 10 ])./mtime(subVH,[1 10 ])), ...
                         std(allvals22(subVH,[1 10 ])./mtime(subVH,[1 10 ]))./sqrt(length(subVH)),'k.')

h=bar([mean(allvals11(subVH,1)./mtime(subVH,1)) mean(allvals11(subVH,10)./mtime(subVH,10)); ...
     mean(allvals22(subVH,1)./mtime(subVH,1)) mean(allvals22(subVH,10)./mtime(subVH,10))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'VP_{1}', 'VP_{2}'})
xlabel('group')
ylabel(plotlabel)
ylim([0 1])
legend(h,{'early','late'})
legend boxoff
title('VH')
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])


figure
set(gcf,'pos',[100 100 300 400])
errorbar(1+[-0.15 0.15], mean(allvals11(subH,[1 10])./mtime(subH,[1 10 ])), ...
                         std(allvals11(subH,[1 10])./mtime(subH,[1 10]))./sqrt(length(subH)),'k.')
hold on
errorbar(2+[-0.15 0.15], mean(allvals22(subH,[1 10])./mtime(subH,[1 10 ])), ...
                         std(allvals22(subH,[1 10])./mtime(subH,[1 10 ]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals11(subH,1)./mtime(subH,1)) mean(allvals11(subH,10)./mtime(subH,10)) ; ...
     mean(allvals22(subH,1)./mtime(subH,1)) mean(allvals22(subH,10)./mtime(subH,10))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'VP_{1}', 'VP_{2}'})
xlabel('group')
ylabel(plotlabel)
ylim([0 1])
legend(h,{'early','late'})
legend boxoff
title('H')
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
%%%%%%%%%%%%%%%%%%%%%%%%%%************************************
%% Plot per-subject, per-trial
alltsvals11 = [];
alltsvals22 = [];
alltsvals12 = [];
alltsvals21 = [];
trials = size(tc11_ts,2);
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('trials')
    ylabel(plotlabel)
    title(sprintf('S%d',subj))
    
    rtime_ts(subj,:) = min([ts1_ts(subj,:); ts2_ts(subj,:)]); 
    mtime_ts(subj,:) = max([te1_ts(subj,:)-rtime_ts(subj,:); te2_ts(subj,:)-rtime_ts(subj,:)]);
    
    % these denote collaboration level...
   % deltatc1_ts(subj,:) = abs(tc21_ts(subj,:)-tc11_ts(subj,:));
   % deltatc2_ts(subj,:) = abs(tc12_ts(subj,:)-tc22_ts(subj,:));
    deltatc1_ts(subj,:) = tc21_ts(subj,:)-tc11_ts(subj,:);
    deltatc2_ts(subj,:) = tc12_ts(subj,:)-tc22_ts(subj,:);
    
    
    line(1:trials,tc11_ts(subj,:)-rtime_ts(subj,:),'col','b','lines',':');
    line(1:trials,tc12_ts(subj,:)-rtime_ts(subj,:),'col','b','marker','s');
    line(1:trials,tc21_ts(subj,:)-rtime_ts(subj,:),'col','r','marker','s');
    line(1:trials,tc22_ts(subj,:)-rtime_ts(subj,:),'col','r','lines',':');
    line(1:trials,mtime_ts(subj,:),'col','k','lines',':');
    
    xlim([1 trials])
    ylim(yrange)
    
    alltsvals11 = [alltsvals11; tc11_ts(subj,:)-rtime_ts(subj,:)];
    alltsvals22 = [alltsvals22; tc22_ts(subj,:)-rtime_ts(subj,:)];
    alltsvals12 = [alltsvals12; tc12_ts(subj,:)-rtime_ts(subj,:)];
    alltsvals21 = [alltsvals21; tc21_ts(subj,:)-rtime_ts(subj,:)];
end

%% Plots average per-trial

% VH group
figure
set(gcf,'pos',[100 100 250 350])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel(plotlabel)
title('VH group')
hold on
%s1=errorbar(1:trials, nanmean(alltsvals11(subVH,:)), nanstd(alltsvals11(subVH,:))./sqrt(sum(isfinite(alltsvals11(subVH,:)))),'bs-')
hold on
s2=errorbar(1:trials, nanmean(tc21_ts(subVH,:)), nanstd(tc21_ts(subVH,:))./sqrt(sum(isfinite(alltsvals21(subVH,:)))),'r')
hold on
s5=errorbar(1:trials, nanmean(tc12_ts(subVH,:)), nanstd(tc12_ts(subVH,:))./sqrt(sum(isfinite(alltsvals12(subVH,:)))),'b')
% s6=errorbar(1:trials, nanmean(mtime_ts(subVH,:)), ...
%                       nanstd(mtime_ts(subVH,:))./sqrt(length(subVH)),'ks:')


%s6=errorbar(1:trials, nanmean(alltsvals22(subVH,:)), nanstd(alltsvals22(subVH,:))./sqrt(sum(isfinite(alltsvals22(subVH,:)))),'rs-')
%t1=legend([s1,s2,s5,s6],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')
t1=legend([s2,s5],'S1-VP2','S2-VP1')
set(t1,'Location','north')
legend boxoff
xlim([1 trials])
%ylim(yrange)

%catchtrials:
center_trials = 6+12*((1:tsets)-1);
hold on
errorbar(center_trials,nanmean(tc21_statmat_nofor(subVH,2:end)),...
    nanstd(tc21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'k')
errorbar(center_trials,nanmean(tc12_statmat_nofor(subVH,2:end)),...
    nanstd(tc12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'k')

%%
% H group
figure
set(gcf,'pos',[100 100 250 350])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel(plotlabel)
title('H group')
%s3=errorbar(1:trials, nanmean(alltsvals11(subH,:),1), nanstd(alltsvals11(subH,:),1)./sqrt(sum(isfinite(alltsvals11(subH,:)))),'cs-')
%hold on
s4=errorbar(1:trials, nanmean(tc21_ts(subH,:),1), nanstd(tc21_ts(subH,:),1)./sqrt(sum(isfinite(alltsvals21(subH,:)))),'rs:')
hold on
patch([1:trials trials:-1:1],[nanmean(tc21_ts(subH,:),1)-nanstd(tc21_ts(subH,:),1) ...
                        fliplr(nanmean(tc21_ts(subH,:),1)+nanstd(tc21_ts(subH,:),1))],0.5+[0.5 0 0],'edgecol',0.5+[0.5 0 0]);
s4=line(1:trials, nanmean(tc21_ts(subH,:),1), 'col','m','marker','s','lines','-')


s7=errorbar(1:trials, nanmean(tc12_ts(subH,:),1), nanstd(tc12_ts(subH,:),1)./sqrt(sum(isfinite(alltsvals12(subH,:)))),'bs:')

s8=errorbar(1:trials, nanmean(mtime_ts(subH,:)), nanstd(mtime_ts(subH,:))./sqrt(length(subH)),'ks:')

%s8=errorbar(1:trials, nanmean(alltsvals22(subH,:),1), nanstd(alltsvals22(subH,:),1)./sqrt(sum(isfinite(alltsvals22(subH,:)))),'ms-')
%t2=legend([s3,s4,s7,s8],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')
t2=legend([s4,s7],'S2-VP1','S1-VP2')
set(t2,'Location','north')
legend boxoff
xlim([1 trials])
ylim(yrange)

%%%%catch trial try!
% catch trials:
errorbar(center_trials,nanmean(tc21_statmat_nofor(subH,2:end)),...
    nanstd(tc21_statmat_nofor(subH,2:end))./sqrt(length(subH)),'k')
errorbar(center_trials,nanmean(tc12_statmat_nofor(subH,2:end)),...
    nanstd(tc12_statmat_nofor(subH,2:end))./sqrt(length(subH)),'k')


%% Plot delta crossing time per epoch (makes little sense...)
% Group H: 
figure
set(gcf,'pos',[100 100 250 350])
patch([1.5 11.5 11.5 1.5],[-2 -2 2 2],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group H')
hold on

e1=errorbar(1:tsets, nanmean(deltatc1(subH,:)), nanstd(deltatc1(subH,:))./sqrt(sum(isfinite(deltatc1(subH,:)))),'bs-')
hold on
e2=errorbar(1:tsets, nanmean(deltatc2(subH,:)), nanstd(deltatc2(subH,:))./sqrt(sum(isfinite(deltatc2(subH,:)))),'rs-')
l2=legend([e1,e2],'VP_{1}','VP_{2}')
line([1 tsets],[0 0],'col','k')
set(l1,'Location','north')
legend boxoff
xlim([1 tsets])
ylim([-2 2])

% Group VH: 
figure
set(gcf,'pos',[100 100 250 350])
patch([1.5 11.5 11.5 1.5],[-2 -2 2 2],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group VH')
hold on

e1=errorbar(1:tsets, nanmean(deltatc1(subVH,:)), nanstd(deltatc1(subVH,:))./sqrt(sum(isfinite(deltatc1(subVH,:)))),'bs-')
hold on
e2=errorbar(1:tsets, nanmean(deltatc2(subVH,:)), nanstd(deltatc2(subVH,:))./sqrt(sum(isfinite(deltatc2(subVH,:)))),'rs-')
l2=legend([e1,e2],'VP_{1}','VP_{2}')
line([1 tsets],[0 0],'col','k')
set(l1,'Location','north')
legend boxoff
xlim([1 tsets])

ylim([-2 2])


%% Plot delta crossing time per trial

figure
set(gcf,'pos',[100 100 250 350])
patch([12.5 132.5 132.5 12.5],[-2 -2 2 2],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group H')
hold on

e1=errorbar(1:trials, nanmean(deltatc1_ts(subH,:)), nanstd(deltatc1_ts(subH,:))./sqrt(sum(isfinite(deltatc1_ts(subH,:)))),'bs-')
hold on
e2=errorbar(1:trials, nanmean(deltatc2_ts(subH,:)), nanstd(deltatc2_ts(subH,:))./sqrt(sum(isfinite(deltatc2_ts(subH,:)))),'rs-')
l2=legend([e1,e2],'VP_{1}','VP_{2}')
line([1 trials],[0 0],'col','k')
set(l1,'Location','north')
legend boxoff
xlim([1 trials])
ylim([-2 2])

figure
set(gcf,'pos',[100 100 250 350])
patch([12.5 132.5 132.5 12.5],[-2 -2 2 2],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group VH')
hold on

e1=errorbar(1:trials, nanmean(deltatc1_ts(subVH,:)), nanstd(deltatc1_ts(subVH,:))./sqrt(sum(isfinite(deltatc1_ts(subVH,:)))),'bs-')
hold on
e2=errorbar(1:trials, nanmean(deltatc2_ts(subVH,:)), nanstd(deltatc2_ts(subVH,:))./sqrt(sum(isfinite(deltatc2_ts(subVH,:)))),'rs-')
l2=legend([e1,e2],'VP_{1}','VP_{2}')
line([1 trials],[0 0],'col','k')
set(l1,'Location','north')
legend boxoff

xlim([1 trials])
ylim([-2 2])




%% Plot crossng times pre-post 

% (own via-point)

figure
set(gcf,'pos',[100 100 250 350])


errorbar(2+[-0.15 0.15], mean(allvals11(subVH,[2 11])), ...
                         std(allvals11(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(allvals11(subH,[2 11])), ...
                         std(allvals11(subH,[2 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals11(subH,2)) mean(allvals11(subH,11)); ...
     mean(allvals11(subVH,2)) mean(allvals11(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel(plotlabel)
ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
title('VP1')

figure
set(gcf,'pos',[100 100 250 350])


errorbar(2+[-0.15 0.15], mean(allvals22(subVH,[2 11])), ...
                         std(allvals22(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(allvals22(subH,[2 11])), ...
                         std(allvals22(subH,[2 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals22(subH,2)) mean(allvals22(subH,11)); ...
     mean(allvals22(subVH,2)) mean(allvals22(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel(plotlabel)
ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
title('VP2')



%% Plot average crossing times, pre-post other via-point
figure
set(gcf,'pos',[100 100 250 350])
errorbar(2+[-0.15 0.15], nanmean(allvals21(subVH,[2 11])), ...
                     nanstd(allvals21(subVH,[2 11]))./sqrt(sum(isfinite(allvals21(subVH,[2 11])))),'k.')
hold on
errorbar(1+[-0.15 0.15], nanmean(allvals21(subH,[2 11])), ...
                 nanstd(allvals21(subH,[2 11]))./sqrt(sum(isfinite(allvals21(subH,[2 11])))),'k.')

h=bar([nanmean(allvals21(subH,2)) nanmean(allvals21(subH,11)); ...
     nanmean(allvals21(subVH,2)) nanmean(allvals21(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel(plotlabel)
ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
title('VP1')

figure
set(gcf,'pos',[100 100 250 350])

errorbar(2+[-0.15 0.15], nanmean(allvals12(subVH,[2 11])), ...
                         nanstd(allvals12(subVH,[2 11]))./sqrt(sum(isfinite(allvals12(subVH,[2 11])))),'k.')
hold on
errorbar(1+[-0.15 0.15], nanmean(allvals12(subH,[2 11])), ...
                         nanstd(allvals12(subH,[2 11]))./sqrt(sum(isfinite(allvals12(subH,[2 11])))),'k.')

h=bar([nanmean(allvals12(subH,2)) nanmean(allvals12(subH,11)); ...
     nanmean(allvals12(subVH,2)) nanmean(allvals12(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel(plotlabel)
ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
title('VP2')

%% Plot delta crossing time
% This is 'late' only (in 'early' very often there is no crossing)

figure
set(gcf,'pos',[100 100 250 350])
errorbar(1, nanmean(deltatc1(subH,11)), ...
            nanstd(deltatc1(subH,11))./sqrt(sum(isfinite(deltatc1(subH,11)))),'k.')
hold on
errorbar(2, nanmean(deltatc1(subVH,11)), ...
            nanstd(deltatc1(subVH,11))./sqrt(sum(isfinite(deltatc1(subVH,11)))),'k.')
hold on
h=bar(1:2, [nanmean(deltatc1(subH,11)) nanmean(deltatc1(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('\Delta tc_{1} [s]')
ylim([-0.2 0.2])
%legend(h,{'early','late'})
%legend boxoff
%set(h(1),'facecol','w')
set(h,'facecol',[0.8 0.8 0.8])
title('VP1')

figure
set(gcf,'pos',[100 100 250 350])
errorbar(1, nanmean(deltatc2(subH,11)), ...
            nanstd(deltatc2(subH,11))./sqrt(sum(isfinite(deltatc2(subH,11)))),'k.')
hold on
errorbar(2, nanmean(deltatc2(subVH,11)), ...
            nanstd(deltatc2(subVH,11))./sqrt(sum(isfinite(deltatc2(subVH,11)))),'k.')
hold on
h=bar(1:2, [nanmean(deltatc2(subH,11)) nanmean(deltatc2(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('\Delta tc_{2} [s]')
ylim([-0.2 0.2])
%legend(h,{'early','late'})
%legend boxoff
%set(h(1),'facecol','w')
set(h,'facecol',[0.8 0.8 0.8])
title('VP2')

%% Stats analysis

group = zeros(subjs,1);
StatMatrix = [group, deltatc1(:,[11]), deltatc2(:,[11])]
StatMatrix(6,:)=[]; % take out subject 6
StatMatrix(subH,1)=1;
StatMatrix(6:end,1)=2;
save('crossingtimes.txt','StatMatrix', '-ascii')


%% Plot delta crossing time computed from per-trial 
% This is 'late' only (in 'early' very often there is no crossing)

figure
set(gcf,'pos',[100 100 250 350])
late_tr = (1:12)+12*(11-1);

deltatc1_late = nanmean(deltatc1_ts(:,late_tr),2);

errorbar(1, nanmean(deltatc1_late(subH)), ...
             nanstd(deltatc1_late(subH))./sqrt(length(subH)),'k.')
hold on
errorbar(2, nanmean(deltatc1_late(subVH)), ...
             nanstd(deltatc1_late(subVH))./sqrt(length(subVH)),'k.')
hold on
h=bar(1:2, [nanmean(deltatc1_late(subH)) nanmean(deltatc1_late(subVH))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('\Delta tc_{1} [s]')
ylim([-0.15 0.15])
%legend(h,{'early','late'})
%legend boxoff
%set(h(1),'facecol','w')
set(h,'facecol',[0.8 0.8 0.8])
title('VP1')

figure
set(gcf,'pos',[100 100 250 350])
late_tr = (1:12)+12*(11-1);

% eliminate outliers...
%deltatc2_ts(:,late_tr(find(deltatc2_ts(:,late_tr) <-0.2)))=nan;

deltatc2_late = nanmean(deltatc2_ts(:,late_tr),2);


errorbar(1, nanmean(deltatc2_late(subH)), ...
             nanstd(deltatc2_late(subH))./sqrt(length(subH)),'k.')
hold on
errorbar(2, nanmean(deltatc2_late(subVH)), ...
             nanstd(deltatc2_late(subVH))./sqrt(length(subVH)),'k.')
hold on
h=bar(1:2, [nanmean(deltatc2_late(subH)) nanmean(deltatc2_late(subVH))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('\Delta tc_{2} [s]')
ylim([-0.15 0.15])
%legend(h,{'early','late'})
%legend boxoff
%set(h(1),'facecol','w')
set(h,'facecol',[0.8 0.8 0.8])
title('VP2')

%% Reaction time, movement time per subjects

trials = size(tc11_ts,2);
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([36.5 180.5 180.5 36.5],[0 0 40 40],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('trials')
    ylabel('time[s]')
    title(sprintf('Dyad %d',subj))
    rtime_ts(subj,:) = min([ts1_ts(subj,:); ts2_ts(subj,:)]); 
    mtime_ts(subj,:) = max([te1_ts(subj,:)-rtime_ts(subj,:); te2_ts(subj,:)-rtime_ts(subj,:)]);
   
    deltatc1_ts(subj,:) = tc21_ts(subj,:)-tc11_ts(subj,:);
    deltatc2_ts(subj,:) = tc12_ts(subj,:)-tc22_ts(subj,:);
   
    a(1)=line(1:trials,mtime_ts(subj,:),'col','k','lines','-');
    a(2)=line(1:trials,rtime_ts(subj,:),'col','r','lines','-');
    a(3)=line(1:trials,ts1_ts(subj,:),'col','b','lines','--');
    a(4)=line(1:trials,ts2_ts(subj,:),'col','r','lines','--');
   % a(5)=line(1:trials,(te1_ts(subj,:)-ts1_ts(subj,:)));
   % a(6)=line(1:trials,(te2_ts(subj,:)-ts2_ts(subj,:)));
    a(5)=line(1:trials,tc11_ts(subj,:),'col','g','lines','--');
    a(6)=line(1:trials,tc12_ts(subj,:),'col','m','lines','--');
    
         %line(1:trials,md11_ts(subj,:),'col','b','marker','.','lines',':');
   % a(5)=line(1:trials,md12_ts(subj,:),'col','b','marker','s','lines','-');
   % a(6)=line(1:trials,md21_ts(subj,:),'col','r','marker','s','lines','-');
         %line(1:trials,md22_ts(subj,:),'col','r','marker','.','lines',':');
 
    legend(a,{'movement time','reaction time','ts_1','ts_2','sub1 mtime','subj2 mtime'})
    legend boxoff
    xlim([1 trials])
    
    
   
end
%%
trials = size(tc11_ts,2);
for subj=1:subjs
    
    figure
    [a,b,c]=plotyy(1:trials,rtime_ts(subj,:),1:trials,mdavg_ts(subj,:))
    xlabel('Trials')
    ylabel(a(1),'reaction time')
    ylabel(a(2),'mdavg_{12}')
    title(sprintf('Dyad %d, rtime.mdavg',subj))
    
    figure
    [a,b,c]=plotyy(1:trials,mtime_ts(subj,:),1:trials,mdavg_ts(subj,:))
    xlabel('Trials')
    ylabel(a(1),'movement time')
    ylabel(a(2),'mdavg_{12}')
    title(sprintf('Dyad %d, mtime.mdavg',subj))

   
end

%%

alltsvals_avg = 0.5*(alltsvals12+alltsvals21)
figure
subplot(2,1,1)

[a b c]=plotyy(1:trials,mean(alltsvals_avg(subH,:)),1:trials,mean(mdavg_ts(subH,:)))

xlabel('trials')
title('H')
ylabel(a(1),'avg crossing times(tc12,tc21)')
ylabel(a(2),'avg minimum dist(tc12,tc21)')
ylim(a(1),[0 4])
ylim(a(2),[0 0.06])
%patch([12.5 132.5 132.5 12.5],[0 0 4 4],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
subplot(2,1,2)
[x b c]=plotyy(1:trials,mean(alltsvals_avg(subVH,:)),1:trials,mean(mdavg_ts(subVH,:)))
xlabel('trials')
title('VH')
ylim(x(1),[0 4])
ylim(x(2),[0 0.06])
%patch([12.5 132.5 132.5 12.5],[0 0 4 4],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
%ylabel(a(1),'avg crossing times(tc12,tc21)')
%ylabel(a(2),'avg minimum dist(tc12,tc21)')


%% to check difference in starting time
trials = size(tc11_ts,2);
for subj=1:subjs
    
    figure
    [a,b,c]=plotyy(1:trials,rtime_ts(subj,:),1:trials,mdavg_ts(subj,:))
    xlabel('Trials')
    ylabel(a(1),'reaction time')
    ylabel(a(2),'mdavg_{12}')
    title(sprintf('Dyad %d, rtime.mdavg',subj))
    
    figure
    [a,b,c]=plotyy(1:trials,mtime_ts(subj,:),1:trials,mdavg_ts(subj,:))
    xlabel('Trials')
    ylabel(a(1),'movement time')
    ylabel(a(2),'mdavg_{12}')
    title(sprintf('Dyad %d, mtime.mdavg',subj))
 
end


