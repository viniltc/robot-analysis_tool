%% Plot tstart and tend for the two groups

clear all
close all

% Set location of data
datadir = '../../Data/';
% Details on where are results
resdir = [datadir, 'results/'];

phases = {'baseline','training','aftereffect'}; 
order = {1,2:11,12:13};  % target sets in experiment 

subH = 1:5;
subVH = 6:9;  % 2 removed


load([resdir, 'ts1_statmat.txt']);
load([resdir, 'ts1_statmat_catch.txt']);
load([resdir, 'ts2_statmat.txt']);
load([resdir, 'ts2_statmat_catch.txt']);
load([resdir, 'te1_statmat.txt']);
load([resdir, 'te1_statmat_catch.txt']);
load([resdir, 'te2_statmat.txt']);
load([resdir, 'te2_statmat_catch.txt']);
group = ts1_statmat(:,1);

ts1 = ts1_statmat(:,2:end);
ts1_c = ts1_statmat_catch();
ts2 = ts2_statmat(:,2:end);
ts2_c = ts2_statmat_catch();
te1 = te1_statmat(:,2:end);
te1_c = te1_statmat_catch();
te2 = te2_statmat(:,2:end);
te2_c = te2_statmat_catch();

% load per-trial data
load ([resdir,'ts1.mat']);
ts1_ts = ind_ts;
load ([resdir,'ts2.mat']);
ts2_ts = ind_ts;
load ([resdir,'te1.mat']);
te1_ts = ind_ts;
load ([resdir,'te2.mat']);
te2_ts = ind_ts;

%%Plot per-subject, per-epoch
subjs=size(ts1,1);
tsets=size(ts1,2);
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([1.5 11.5 11.5 1.5],[0 0 6 6],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel('time [s]')
    title(sprintf('S%d',subj))
        line(1:tsets,ts2(subj,:),'col','b','marker','s','lines',':');
        line(1:tsets,te2(subj,:),'col','b','marker','s');
        line(1:tsets,ts1(subj,:),'col','r','marker','s','lines',':');
        line(1:tsets,te1(subj,:),'col','r','marker','s');
    xlim([1 tsets])
    ylim([0 6])
end

% same but subtracting the common reaction time
allvals1 = [];
allvals2 = [];
allvals1_c = [];
allvals2_c = [];
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([1.5 11.5 11.5 1.5],[0 0 6 6],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel('duration [s]')
    title(sprintf('S%d',subj))
    
    rtime = min([ts1(subj,:); ts2(subj,:)]);
    rtime_c = min([ts1_c(subj,:); ts2_c(subj,:)]);
    
    
    allvals1 = [allvals1; te1(subj,:)-rtime];
    allvals2 = [allvals2; te2(subj,:)-rtime];
    
    allvals1_c = [allvals1_c; te1_c(subj,:)-rtime_c];
    allvals2_c = [allvals2_c; te2_c(subj,:)-rtime_c];
    
    line(1:tsets,te1(subj,:)-rtime, 'col','b','marker','s','lines',':');
    line(1:tsets,te2(subj,:)-rtime, 'col','r','marker','s','lines',':');
  xlim([1 tsets])
    ylim([0 4])
end  
    
%%Plot average, per-epoch
figure
set(gcf,'pos',[100 100 400 500])
patch([1.5 11.5 11.5 1.5],[0 0 6 6],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel('duration [s]')
title('all subjects')
hold on

e1=errorbar(1:tsets, mean(allvals1(subVH,:)), std(allvals1(subVH,:))./sqrt(length(subVH)),'rs-')
hold on
e2=errorbar(1:tsets, mean(allvals2(subVH,:)), std(allvals2(subVH,:))./sqrt(length(subVH)),'bs-')
hold on
e3=errorbar(1:tsets, mean(allvals1(subH,:),1), std(allvals1(subH,:),1)./sqrt(length(subH)),'ms-')
e4=errorbar(1:tsets, mean(allvals2(subH,:),1), std(allvals2(subH,:),1)./sqrt(length(subH)),'cs-')
l1=legend([e1,e2,e3,e4],'VH:Subj 1','VH:Subj 2','H: Subj 1','H:Subj 2')
legend boxoff
set(l1,'Location','north')

%% Plot per-subject, per-trial
alltsvals1 = [];
alltsvals2 = [];
trials = size(ts1_ts,2);
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([12.5 132.5 132.5 12.5],[0 0 6 6],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('trials')
    ylabel('duration[s]')
    title(sprintf('S%d',subj))
        line(1:trials,ts1_ts(subj,:),'col','r','marker','s');
        line(1:trials,ts2_ts(subj,:),'col','b','marker','s');
    xlim([1 trials])
    ylim([0 6])
    
    rtime_ts = min([ts1_ts(subj,:); ts2_ts(subj,:)]);
    
    alltsvals1 = [alltsvals1; te1_ts(subj,:)-rtime_ts];
    alltsvals2 = [alltsvals2; te2_ts(subj,:)-rtime_ts];
end



%% Plot average, pre-post
figure
set(gcf,'pos',[100 100 300 400])
errorbar(2+[-0.15 0.15], mean(allvals1(subVH,[2 11])), ...
                         std(allvals1(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(allvals1(subH,[2 11])), ...
                         std(allvals1(subH,[2 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals1(subH,2)) mean(allvals1(subH,11)); ...
     mean(allvals1(subVH,2)) mean(allvals1(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('duration[s]')
%ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])

figure
set(gcf,'pos',[100 100 300 400])
errorbar(2+[-0.15 0.15], mean(allvals2(subVH,[2 11])), ...
                         std(allvals2(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(allvals2(subH,[2 11])), ...
                         std(allvals2(subH,[2 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals2(subH,2)) mean(allvals2(subH,11)); ...
     mean(allvals2(subVH,2)) mean(allvals2(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('duration[s]')
%ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])


%% Plots average per-trial
figure
set(gcf,'pos',[100 100 300 400])
patch([12.5 132.5 132.5 12.5],[0 0 6 6],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel('duration[s]')
title('all subjects')
hold on


errorbar(1:trials, nanmean(alltsvals1(subVH,:)), nanstd(alltsvals1(subVH,:))./sqrt(sum(isfinite(alltsvals1(subVH,:)))),'rs-')
hold on
errorbar(1:trials, nanmean(alltsvals2(subVH,:)), nanstd(alltsvals2(subVH,:))./sqrt(sum(isfinite(alltsvals1(subVH,:)))),'bs-')
hold on
errorbar(1:trials, nanmean(alltsvals1(subH,:),1), nanstd(alltsvals1(subH,:),1)./sqrt(sum(isfinite(alltsvals1(subH,:)))),'ms-')
errorbar(1:trials, nanmean(alltsvals2(subH,:),1), nanstd(alltsvals2(subH,:),1)./sqrt(sum(isfinite(alltsvals1(subH,:)))),'cs-')

xlim([1 trials])
ylim([0 6])


%% Plots average per- catch trial
% figure
% set(gcf,'pos',[100 100 300 400])
% patch([12.5 132.5 132.5 12.5],[0 0 6 6],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% xlabel('trials')
% ylabel('duration[s]')
% title('all subjects')
% hold on
% 
% 
% errorbar(1:trials, nanmean(alltsvals1_c(subVH,:)), nanstd(alltsvals1_c(subVH,:))./sqrt(sum(isfinite(alltsvals1_c(subVH,:)))),'rs-')
% hold on
% errorbar(1:trials, nanmean(alltsvals2_c(subVH,:)), nanstd(alltsvals2_c(subVH,:))./sqrt(sum(isfinite(alltsvals1_c(subVH,:)))),'bs-')
% hold on
% errorbar(1:trials, nanmean(alltsvals1_c(subH,:),1), nanstd(alltsvals1_c(subH,:),1)./sqrt(sum(isfinite(alltsvals1_c(subH,:)))),'ms-')
% errorbar(1:trials, nanmean(alltsvals2_c(subH,:),1), nanstd(alltsvals2_c(subH,:),1)./sqrt(sum(isfinite(alltsvals1_c(subH,:)))),'cs-')
% 
% xlim([1 trials])
% ylim([0 6])
    
 %%for catch trials

% figure
% set(gcf,'pos',[100 100 300 400])
% errorbar(2+[-0.15 0.15], mean(allvals1_c(subVH,[2 11])), ...
%                          std(allvals1_c(subVH,[2 11]))./sqrt(length(subVH)),'k.')
% hold on
% errorbar(1+[-0.15 0.15], mean(allvals1_c(subH,[2 11])), ...
%                          std(allvals1_c(subH,[2 11]))./sqrt(length(subH)),'k.')
% 
% j=bar([mean(allvals1_c(subH,2)) mean(allvals1_c(subH,11)); ...
%      mean(allvals1_c(subVH,2)) mean(allvals1_c(subVH,11))])
% box off
% set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
% xlabel('group')
% ylabel('duration[s]')
% %ylim(yrange)
% legend(j,{'early','late'})
% legend boxoff
% set(j(1),'facecol','w')
% set(j(2),'facecol',[0.8 0.8 0.8])
% 
% 
% figure
% set(gcf,'pos',[100 100 300 400])
% errorbar(2+[-0.15 0.15], mean(allvals2_c(subVH,[2 11])), ...
%                          std(allvals2_c(subVH,[2 11]))./sqrt(length(subVH)),'k.')
% hold on
% errorbar(1+[-0.15 0.15], mean(allvals2_c(subH,[2 11])), ...
%                          std(allvals2_c(subH,[2 11]))./sqrt(length(subH)),'k.')
% 
% j=bar([mean(allvals2_c(subH,2)) mean(allvals2_c(subH,11)); ...
%      mean(allvals2_c(subVH,2)) mean(allvals2_c(subVH,11))])
% box off
% set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
% xlabel('group')
% ylabel('duration[s]')
% %ylim(yrange)
% legend(j,{'early','late'})
% legend boxoff
% set(j(1),'facecol','w')
% set(j(2),'facecol',[0.8 0.8 0.8])
%  