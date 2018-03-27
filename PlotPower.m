%% Plot power figures

clear all
close all

% Set location of data
datadir = '../../Data/';
% Details on where are results
resdir = [datadir, 'results/'];

phases = {'baseline','training','aftereffect'}; 
order = {1,2:11,12:13};  % target sets in experiment 

subH = 1:5;
subVH = 7:11;  % S6 removed

% load per-epoch data
load([resdir, 'pow1_statmat.txt']);
load([resdir, 'pow1_statmat_catch.txt']);
load([resdir, 'pow2_statmat.txt']);
load([resdir, 'pow2_statmat_catch.txt']);

% load per-trial data
load ([resdir,'pow1.mat']);
pow1_ts = (ind_ts);
load ([resdir,'pow2.mat']);
pow2_ts = (ind_ts);
load([resdir, 'forces.mat']);

pow1_ts(find(forces==0))= nan;
pow2_ts(find(forces==0))= nan;

plotlabel = 'power (w)';
yrange = [0 10]; % range of variation of effort

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

rtime = 0*ts1;

rt_ts = min(ts1_ts,ts2_ts);
mt1_ts = te1_ts-rt_ts;
mt2_ts = te2_ts-rt_ts;
mt_ts = max(mt1_ts,mt2_ts);


% convert NaN to 0

pow1_ts(find(isnan(pow1_ts)))=0;
pow2_ts(find(isnan(pow2_ts)))=0;



%% Plot per-subject, per-epoch
pow1 = sqrt(pow1_statmat(:,2:end));
pow2 = sqrt(pow2_statmat(:,2:end));
subjs=size(pow1,1);
tsets=size(pow1,2);
allvals1 = [];
allvals2 = [];
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel(plotlabel)
    title(sprintf('dyad%d',subj))
        line(1:tsets,pow1(subj,:),'col','r','marker','s');
        line(1:tsets,pow2(subj,:),'col','b','marker','s');
    xlim([1 tsets])
    %ylim(yrange)
    
    allvals1 = [allvals1; pow1(subj,:)];
    allvals2 = [allvals2; pow2(subj,:)];
    
  
    
end



%% Plots average per-epoch
figure
set(gcf,'pos',[100 100 300 350])
 patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('all subjects')
hold on

e1=errorbar(1:tsets, mean(allvals1(subVH,:)), std(allvals1(subVH,:))./sqrt(length(subVH)),'rs-')
hold on
e2=errorbar(1:tsets, mean(allvals2(subVH,:)), std(allvals2(subVH,:))./sqrt(length(subVH)),'bs-')
hold on
e3=errorbar(1:tsets, mean(allvals1(subH,:),1), std(allvals1(subH,:),1)./sqrt(length(subH)),'ms-')
e4=errorbar(1:tsets, mean(allvals2(subH,:),1), std(allvals2(subH,:),1)./sqrt(length(subH)),'cs-')



%errorbar(1:tsets, allvals1(subH,:), 0*allvals1(subH,:),'rs-')
%errorbar(1:tsets, allvals2(subH,:), 0*allvals2(subH,:),'bs-')
l1=legend([e1,e2,e3,e4],'VH:Subj 1','VH:Subj 2','H: Subj 1','H: Subj 2')
%rect=[1.25, 1.25, .25, .zeros25]
set(l1,'Location','north')
legend boxoff
% ylim(yrange)


%% Plot average, pre-post
figure
set(gcf,'pos',[100 100 250 350])
title('Subject 1')
hold on
errorbar(2+[-0.225 0 0.225], mean(allvals1(subVH,[2 6 11])), ...
                         std(allvals1(subVH,[2 6 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(allvals1(subH,[2 6 11])), ...
                         std(allvals1(subH,[2 6 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals1(subH,2)) mean(allvals1(subH,6)) mean(allvals1(subH,11)); ...
     mean(allvals1(subVH,2)) mean(allvals1(subVH,6)) mean(allvals1(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel(plotlabel)
%ylim([0 4])
%ylim(yrange)
legend(h,{'early','mid','late'})
 legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])


figure
set(gcf,'pos',[100 100 250 350])
title('Subject 2')
hold on
errorbar(2+[-0.225 0 0.225], mean(allvals2(subVH,[2 6 11])), ...
                         std(allvals2(subVH,[2 6 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(allvals2(subH,[2 6 11])), ...
                         std(allvals2(subH,[2  6 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals2(subH,2)) mean(allvals2(subH,6)) mean(allvals2(subH,11)); ...
     mean(allvals2(subVH,2)) mean(allvals2(subVH,6)) mean(allvals2(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel(plotlabel)
%ylim([0 4])
%ylim(yrange)
% legend(h,{'early','mid','late'})
 legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])




%% Sum power

NP = abs(pow1+ pow2);
allvals_sum = [];
NP_H = NP(subH,11);
NP_VH = NP(subVH,11);
NP_tot = [NP_H; NP_VH];

for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 350])
    patch([1.5 11.5 11.5 1.5],[0 0 30 30],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel('sum power [W]')
    title(sprintf('D%d',subj))
        h(1)=line(1:tsets,NP(subj,:),'col','r','marker','s');
        h(2)=line(1:tsets,NP(subj,:),'col','b','marker','s');
    xlim([1 tsets])
    ylim([0 30])
    
    allvals_sum = [allvals_sum; NP(subj,:)];
   
end

 %% Plots average per-epoch

figure
set(gcf,'pos',[100 100 300 350])
patch([1.5 11.5 11.5 1.5],[0 0 30 30],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel('net power [W]')

hold on
e2=errorbar(1:tsets, mean(allvals_sum(subH,:)), std(allvals_sum(subH,:))./sqrt(length(subH)),'bs-');
%e3=errorbar(1:tsets_d, mean(allvals_sum(dyads_a,:),1), std(allvals_sum(dyads_a,:),1)./sqrt(length(dyads_a)),'ms-');
e4=errorbar(1:tsets, mean(allvals_sum(subVH,:),1), std(allvals_sum(subVH,:),1)./sqrt(length(subVH)),'rs-');
l1=legend([e2,e4],'H','VH');

set(l1,'Location','northeast')
%legend boxoff
ylim([0 30])
xlim([1 tsets])

%% Plot average, pre-post
figure
set(gcf,'pos',[100 100 250 350])
% title('Subject 1')
hold on
errorbar(2+[-0.225 0 0.225], mean(allvals_sum(subVH,[2 6 11])), ...
                         std(allvals_sum(subVH,[2 6 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(allvals_sum(subH,[2 6 11])), ...
                         std(allvals_sum(subH,[2 6 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals_sum(subH,2)) mean(allvals_sum(subH,6)) mean(allvals_sum(subH,11)); ...
     mean(allvals_sum(subVH,2)) mean(allvals_sum(subVH,6)) mean(allvals_sum(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('Sum power [N]')
%ylim([0 4])
%ylim(yrange)
legend(h,{'early','mid','late'})
 legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])


% %% Plot per-subject, per-trial
% % alltsvals1 = [];
% % alltsvals2 = [];
% trials = size(pow1_ts,2);
% for subj=1:subjs
%     figure
%     set(gcf,'pos',[100 100 300 400])
%     patch([12.5 132.5 132.5 12.5],[-30 -30 40 40],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
%     xlabel('trials')
%     ylabel(plotlabel)
%     title(sprintf('S%d',subj))
%         line(1:trials,pow1_ts(subj,:),'col','r','marker','s');
%         line(1:trials,pow2_ts(subj,:),'col','b','marker','s');
%     xlim([1 trials])
%     %ylim(yrange)
%     
% %     alltsvals1 = [alltsvals1; pow1_ts(subj,:)];
% %     alltsvals2 = [alltsvals2; pow2_ts(subj,:)];
% end


% %% Plot leadership index
% trials = size(pow1_ts,2);
% for subj=1:subjs
%     figure
%     set(gcf,'pos',[100 100 300 400])
%     patch([12.5 132.5 132.5 12.5],[-15 -15 15 15],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
%     xlabel('trials')
%     ylabel('duration(s)')
%     title(sprintf('S%d',subj))
%     hold on
% 
% %   line(1:trials,sign(pow1_ts(subj,:)),'col','b')
% %   line(1:trials,sign(pow2_ts(subj,:)),'col','r')
% % 
% %   bar(12:132,sign(pow1_ts(subj,12:132)).*mt_ts(subj,12:132),'g' ), hold on
% %   bar(12:132,sign(pow2_ts(subj,12:132)).*mt_ts(subj,12:132),'m' )
%     
%     bar(1:trials,sign(pow1_ts(subj,:)).*mt_ts(subj,:),'g' ), hold on
%     bar(1:trials,sign(pow2_ts(subj,:)).*mt_ts(subj,:),'m' )
%     xlim([1 trials])
%     %ylim(yrange)
%     
% %     alltsvals1 = [alltsvals1; pow1_ts(subj,:)];
% %     alltsvals2 = [alltsvals2; pow2_ts(subj,:)];
% end



