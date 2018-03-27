%% Plot md11,md12,md21,md22 for the two groups

% NB: P1 is red, P2 is blue!!!

clear all
close all

% Set location of data
datadir = '../../Data/';
% Details on where are results
resdir = [datadir, 'results/'];
epsdir = [datadir, 'eps/'];
figdir = [datadir, 'figs/'];

folder = 'C:\Users\Asus\Desktop\journal_photoshop\illustrator\exp';
folder = 'C:\Users\Asus\Dropbox\pHHI sfn paper\Data\eps\new\panel_sep\';
stat_folder = 'C:\Users\Asus\Dropbox\pHHI sfn paper\Statistics\';

subH = 1:5;
subVH = 7:11;  % 5 removed
subPV = 12:16;  % 

load([resdir, 'md11_statmat.txt']);
load([resdir, 'md11_statmat_nofor.txt']);
load([resdir, 'md12_statmat.txt']);
load([resdir, 'md12_statmat_nofor.txt']);
load([resdir, 'md21_statmat.txt']);
load([resdir, 'md21_statmat_nofor.txt']);
load([resdir, 'md22_statmat.txt']);
load([resdir, 'md22_statmat_nofor.txt']);


load([resdir, 'eff1_statmat.txt']);
load([resdir, 'eff1_statmat_catch.txt']);
load([resdir, 'eff2_statmat.txt']);
load([resdir, 'eff2_statmat_catch.txt']);
load([resdir, 'forcenorm_statmat.txt']);
load([resdir, 'forcenorm_statmat_catch.txt']);


load ([resdir,'forcenorm.mat']);
forcenorm_ts = sqrt(ind_ts);
load ([resdir,'eff1.mat']);
eff1_ts = sqrt(ind_ts);
load ([resdir,'eff2.mat']);
eff2_ts = sqrt(ind_ts);
load([resdir, 'forces.mat']);

load('asd_score_arr.txt'); % ASD score for 22 subjects
asd_overall = asd_score_arr(:,6);
asd_subj1 = [asd_score_arr(1,:);asd_score_arr(3,:);asd_score_arr(5,:);asd_score_arr(7,:);asd_score_arr(9,:);asd_score_arr(11,:);asd_score_arr(13,:);asd_score_arr(15,:);asd_score_arr(17,:);asd_score_arr(19,:);asd_score_arr(21,:)]
asd_subj2 = [asd_score_arr(2,:);asd_score_arr(4,:);asd_score_arr(6,:);asd_score_arr(8,:);asd_score_arr(10,:);asd_score_arr(12,:);asd_score_arr(14,:);asd_score_arr(16,:);asd_score_arr(18,:);asd_score_arr(20,:);asd_score_arr(22,:)]



group = md11_statmat(:,1);

% load per-trial data
load ([resdir,'md11.mat']);
md11_ts = (ind_ts);
load ([resdir,'md12.mat']);
md12_ts = (ind_ts);
load ([resdir,'md21.mat']);
md21_ts = (ind_ts);
load ([resdir,'md22.mat']);
md22_ts = (ind_ts);
load([resdir, 'forces.mat']);
forces(4, [30, 47:48]) = 1;

forces_c = forces(:, [13:120]); %% catch trial in phase training

md11_ts_cc = md11_ts;
md12_ts_cc = md12_ts;
md21_ts_cc = md21_ts;
md22_ts_cc = md22_ts;

% md11_ts(find(forces_c==0))= nan;
% md12_ts(find(forces_c==0))= nan;
% md21_ts(find(forces_c==0))= nan;
% md22_ts(find(forces_c==0))= nan;


%load catch trial
% load ([resdir,'md11_bidir_catch.mat']);
% md11c_ts = sqrt(ind_ts);
% load ([resdir,'md12_bidir_catch.mat']);
% md12c_ts = sqrt(ind_ts);
% load ([resdir,'md21_bidir_catch.mat']);
% md21c_ts = sqrt(ind_ts);
% load ([resdir,'md22_bidir_catch.mat']);
% md22c_ts = sqrt(ind_ts);
% load([resdir, 'forces.mat']);

%md11_ts(find(forces==0))= nan;
%md22_ts(find(forces==0))= na12*3n;

%catch trial



md11 = md11_statmat(:,2:end);
md12 = md12_statmat(:,2:end);
md21 = md21_statmat(:,2:end);
md22 = md22_statmat(:,2:end);

forcenorm = forcenorm_statmat(:,2:end);


md_avg = 0.5*(md12+md21); % minimumdistance average

md = [md12,md21]

plotlabel = 'minimum distance [m]';
yrange = [0 0.06]; 

subjs=size(md11,1);
tsets=size(md22,2);
allvals11 = [];
allvals22 = [];
allvals12 = [];
allvals21 = [];
avgallvalsij = [];
for subj=1:subjs
    figure
    %patch([1:tsets tsets:-1:1],[ts2(subj,:) fliplr(te2(subj,:))],'b');
    %patch([1:tsets tsets:-1:1],[ts1(subj,:) fliplr(te1(subj,:))],'r');
    set(gcf,'pos',[100 100 250 350])
    patch([1.5 11.5 11.5 1.5],[0 0 6 6],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel(plotlabel)
    title(sprintf('S%d',subj))
    h(1)=line(1:tsets,md11(subj,:),'col','b','lines',':');
    h(2)=line(1:tsets,md12(subj,:),'col','b');
    h(3)=line(1:tsets,md21(subj,:),'col','r');
    h(4)=line(1:tsets,md22(subj,:),'col','r','lines',':');
    legend(h,{'md_{11}','md_{12}','md_{21}','md_{22}'})
     
    allvals11 = [allvals11; md11(subj,:)];
    allvals22 = [allvals22; md22(subj,:)];
    allvals12 = [allvals12; md12(subj,:)];
    allvals21 = [allvals21; md21(subj,:)];
    avgallvalsij = [avgallvalsij; md_avg(subj,:)];
    
    xlim([1 13])
    ylim(yrange)
   
    
end


%% stats analysis
group = zeros(subjs,1);
StatMatrix = [group, allvals12(:,[2 11]), allvals21(:,[2 11])];
StatMatrix(6,:)=[]; % take out subject 6
StatMatrix(subH,1)=1;
StatMatrix(6:end,1)=2;
save('mindistances_early_late_epoch.txt','StatMatrix', '-ascii')
    
%% Plot average, per-epoch 
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0 0 0.05 0.05],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
 xlabel('Epochs', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('VH group','FontName', 'Times New Roman')
hold on

e2=errorbar(1:tsets, mean(allvals21(subVH,:)), std(allvals21(subVH,:))./sqrt(length(subVH)),'rs-','LineWidth',1)
hold on
e5=errorbar(1:tsets, mean(allvals12(subVH,:)), std(allvals12(subVH,:))./sqrt(length(subVH)),'bs-','LineWidth',1)


% catch trials:
e3=errorbar(1:tsets,nanmean(md21_statmat_nofor(subVH,2:end)),...
    nanstd(md21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'r:','LineWidth',1)
e4=errorbar(1:tsets,nanmean(md12_statmat_nofor(subVH,2:end)),...
    nanstd(md12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'b:','LineWidth',1)

%l1=legend([e2,e5,e3,e4],'md_{21}','md_{12}','catch md_{21}','catch md_{12}')
%set(l1,'Location','northeast')
legend boxoff
ylim([0.005 0.05])
%%Plot average, per-epoch (12,21)
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0 0 0.05 0.05],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
 xlabel('Epochs', 'fontsize', 14,'FontName', 'Times New Roman')
 ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('H group','FontName', 'Times New Roman')
hold on

e(1)=errorbar(1:tsets, mean(allvals21(subH,:),1), std(allvals21(subH,:),1)./sqrt(length(subH)),'rs-','LineWidth',1)
hold on
e(2)=errorbar(1:tsets, mean(allvals12(subH,:),1), std(allvals12(subH,:),1)./sqrt(length(subH)),'bs-','Linewidth',1)
legend boxoff

% catch trials:
e(3)=errorbar(1:tsets,nanmean(md21_statmat_nofor(subH,2:end)),...
    nanstd(md21_statmat_nofor(subH,2:end))./sqrt(length(subH)),'r:','LineWidth',1)
e(4)=errorbar(1:tsets,nanmean(md12_statmat_nofor(subH,2:end)),...
    nanstd(md12_statmat_nofor(subH,2:end))./sqrt(length(subH)),'b:','LineWidth',1)
% legend(e,{'Partner 2','Partner 1','Catch 2','Catch 1'})
legend([e(2), e(1)],{'Subject 1','Subject 2'},'FontName', 'Times New Roman')
legend boxoff
ylim([0.005 0.05])




%% Plot average, pre-post (11,22)
figure
set(gcf,'pos',[100 100 300 400])
title(plotlabel)
errorbar(2+[-0.15 0.15], mean(allvals11(subVH,[2 11])), ...
                         std(allvals11(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(2+[-0.15 0.15], mean(allvals11(subH,[2 11])), ...
                         std(allvals11(subH,[2 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals11(subH,2)) mean(allvals11(subH,11)); ...
     mean(allvals11(subVH,2)) mean(allvals11(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel(plotlabel)
ylim([0 0.04])
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])


figure
set(gcf,'pos',[100 100 300 400])
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
ylim([0 0.04])
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])


%% Plot average, pre-post (12,21)
figure
set(gcf,'pos',[100 100 250 300])
%title('md_{12}')
title('Partner 1')
hold on
errorbar(2+[-0.225 0 0.225], mean(allvals12(subVH,[2 6 11])), ...
                         std(allvals12(subVH,[2 6 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(allvals12(subH,[2 6 11])), ...
                         std(allvals12(subH,[2 6 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals12(subH,2)) mean(allvals12(subH,6)) mean(allvals12(subH,11)); ...
     mean(allvals12(subVH,2)) mean(allvals12(subVH,6)) mean(allvals12(subVH,11))])
 
sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225]},[0.0005, 0.0005])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group', 'fontsize', 14)
ylabel(plotlabel, 'fontsize', 14)
ylim([0 0.05])
legend(h,{'early','mid','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])


figure
%title('md_{21}')
title('Partner 2')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(allvals21(subVH,[2 6 11])), ...
                         std(allvals21(subVH,[2 6 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(allvals21(subH,[2 6 11])), ...
                         std(allvals21(subH,[2 6 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals21(subH,2)) mean(allvals21(subH,6)) mean(allvals21(subH,11)); ...
     mean(allvals21(subVH,2)) mean(allvals21(subVH,6)) mean(allvals21(subVH,11))])
 
 
sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225]},[0.0005, 0.0005])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group', 'fontsize', 14)
ylabel(plotlabel, 'fontsize', 14)
ylim([0 0.05])
%legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])

%% avg pre post

figure
%title('md_{21}')
% title('Partner 2')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(avgallvalsij(subVH,[2 6 11])), ...
                         std(avgallvalsij(subVH,[2 6 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(avgallvalsij(subH,[2 6 11])), ...
                         std(avgallvalsij(subH,[2 6 11]))./sqrt(length(subH)),'k.')

h=bar([mean(avgallvalsij(subH,2)) mean(avgallvalsij(subH,6)) mean(avgallvalsij(subH,11)); ...
     mean(avgallvalsij(subVH,2)) mean(avgallvalsij(subVH,6)) mean(avgallvalsij(subVH,11))])
 
 
sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225], [1,2]},[0.0005, 0.0005, 0.05])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
 xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.06])
legend(h,{'Early','Middle','Late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])




%% Stats analysis (epoch wise)

group = zeros(subjs,1);
StatMatrix = [group, allvals12(:,[2 6 11]), allvals21(:,[2 6 11]), avgallvalsij(:,[2 6 11])];
StatMatrix(6,:)=[]; % take out subject 6
StatMatrix(subH,1)=1;
StatMatrix(6:end,1)=2;
save('mindist_epoch1.txt','StatMatrix', '-ascii')

%% catch trail 

group = zeros(subjs,1);
StatMatrix3 = [group, md12_statmat_nofor(:,[3 6 10]), md21_statmat_nofor(:,[3 6 10])];
StatMatrix3(6,:)=[]; % take out subject 6
StatMatrix3(subH,1)=1;
StatMatrix3(6:end,1)=2;
save('mindist_catch_epoch1.txt','StatMatrix3', '-ascii')

%% Plot per-subject, per-trial
alltsvals11 = [];
alltsvals22 = [];
alltsvals12 = [];
alltsvals21 = [];
trials = size(md11_ts,2);
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 250 350])
    patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('trials')
    ylabel(plotlabel)
    title(sprintf('S%d',subj))
        line(1:trials,md11_ts(subj,:),'col','b','marker','.','lines',':');
        line(1:trials,md12_ts(subj,:),'col','b','marker','s','lines','-');
        line(1:trials,md21_ts(subj,:),'col','r','marker','s','lines','-');
        line(1:trials,md22_ts(subj,:),'col','r','marker','.','lines',':');
    xlim([1 trials])
    ylim(yrange)
    
    alltsvals11 = [alltsvals11; md11_ts(subj,:)];
    alltsvals22 = [alltsvals22; md22_ts(subj,:)];
    alltsvals12 = [alltsvals12; md12_ts(subj,:)];
    alltsvals21 = [alltsvals21; md21_ts(subj,:)];
end


%% Stats analysis (trial wise)

group = zeros(subjs,1);
StatMatrix1 = [group, alltsvals12(:,[(13:24) (121:132)]), alltsvals21(:,[(13:24) (121:132)])];
StatMatrix1(6,:)=[]; % take out subject 6
StatMatrix1(subH,1)=1;
StatMatrix1(6:end,1)=2;
save('mindistances_early_late_trial.txt','StatMatrix1', '-ascii')

%% Plots average per-trial
figure
set(gcf,'pos',[100 100 250 350])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel(plotlabel)
title('VH group')
hold on
%s1=errorbar(1:trials, nanmean(alltsvals11(subVH,:)), nanstd(alltsvals11(subVH,:))./sqrt(sum(isfinite(alltsvals11(subVH,:)))),'bs-')
hold on
s2=errorbar(1:trials, nanmean(md21_ts(subVH,:)), nanstd(md21_ts(subVH,:))./sqrt(sum(isfinite(alltsvals21(subVH,:)))),'rs:')
hold on
s5=errorbar(1:trials, nanmean(md12_ts(subVH,:)), nanstd(md12_ts(subVH,:))./sqrt(sum(isfinite(alltsvals12(subVH,:)))),'bs:')
%s6=errorbar(1:trials, nanmean(alltsvals22(subVH,:)), nanstd(alltsvals22(subVH,:))./sqrt(sum(isfinite(alltsvals22(subVH,:)))),'rs-')
%t1=legend([s1,s2,s5,s6],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')
t1=legend([s2,s5],'md_{21}','md_{12}')
set(t1,'Location','north')
legend boxoff
xlim([1 trials])
ylim(yrange)

% catch trials:
center_trials = 6+12*((1:tsets)-1);
hold on
errorbar(center_trials,nanmean(md21_statmat_nofor(subVH,2:end)),...
    nanstd(md21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'k')
errorbar(center_trials,nanmean(md12_statmat_nofor(subVH,2:end)),...
    nanstd(md12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),'k')



figure
set(gcf,'pos',[100 100 250 350])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel(plotlabel)
title('H group')
hold on

s4=errorbar(1:trials, nanmean(md21_ts(subH,:),1), nanstd(md21_ts(subH,:),1)./sqrt(sum(isfinite(alltsvals21(subH,:)))),'rs:')



s7=errorbar(1:trials, nanmean(md12_ts(subH,:),1), nanstd(md12_ts(subH,:),1)./sqrt(sum(isfinite(alltsvals12(subH,:)))),'bs:')

t2=legend([s4,s7],'md_{21}','md_{12}')
set(t2,'Location','north')
legend boxoff
xlim([1 trials])
ylim(yrange)

%%%%catch trial try!
% catch trials:
errorbar(center_trials,nanmean(md21_statmat_nofor(subH,2:end)),...
    nanstd(md21_statmat_nofor(subH,2:end))./sqrt(length(subH)),'k')
errorbar(center_trials,nanmean(md12_statmat_nofor(subH,2:end)),...
    nanstd(md12_statmat_nofor(subH,2:end))./sqrt(length(subH)),'k')


%% patch try per trial

figure
set(gcf,'pos',[100 100 250 300])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel(plotlabel)
title('H group')
hold on
patch([1:trials trials:-1:1],[nanmean(md21_ts(subH,:),1)-nanstd(md21_ts(subH,:),1) ...
                        fliplr(nanmean(md21_ts(subH,:),1)+nanstd(md21_ts(subH,:),1))],0.5+[0.5 0 0.5],'edgecol',0.5+[0.5 0 0.5]);

patch([1:trials trials:-1:1],[nanmean(md12_ts(subH,:),1)-nanstd(md12_ts(subH,:),1) ...
                        fliplr(nanmean(md12_ts(subH,:),1)+nanstd(md12_ts(subH,:),1))],0.3+[0.3 0 0.3],'edgecol',0.3+[0.3 0 0.3]);
                    legend boxoff
xlim([1 trials])
ylim(yrange)

figure
set(gcf,'pos',[100 100 250 350])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel(plotlabel)
title('VH group')
hold on
patch([1:trials trials:-1:1],[nanmean(md21_ts(subVH,:),1)-nanstd(md21_ts(subVH,:),1) ...
                        fliplr(nanmean(md21_ts(subVH,:),1)+nanstd(md21_ts(subVH,:),1))],0.5+[0.5 0 0.5],'edgecol',0.5+[0.5 0 0.5]);

patch([1:trials trials:-1:1],[nanmean(md12_ts(subVH,:),1)-nanstd(md12_ts(subVH,:),1) ...
                        fliplr(nanmean(md12_ts(subVH,:),1)+nanstd(md12_ts(subVH,:),1))],0.3+[0.3 0 0.3],'edgecol',0.3+[0.3 0 0.3]);
                    legend boxoff
xlim([1 trials])
ylim(yrange)

%% shaded error bar with trasperancy using fcn (see shadedErrorBar function in the directory)

name = sprintf('md_H');
fullFileName = fullfile(folder, name);
figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.0005 0.0005 0.075 0.075],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('H', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
e1=shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subH,:),1)),smooth(nanstd(alltsvals12(subH,:),1))./sqrt(length(subH)),{'-','Color',[.1 .4 .9], 'LineWidth',1},1);
hold on
e2=shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subH,:),1)),smooth(nanstd(alltsvals21(subH,:),1))./sqrt(length(subH)),{'-','Color',[.9 .1 .4], 'LineWidth',1},1);


% catch trials..
% errorbar(center_trials,nanmean(md21_statmat_nofor(subVH,2:end)),...
%     nanstd(md21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)*0),'Color',[.9 .4 .4],'LineStyle',':', 'LineWidth',2)
% errorbar(center_trials,nanmean(md12_statmat_nofor(subVH,2:end)),...
%     nanstd(md12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)*0),'Color',[.4 .4 .9],'LineStyle',':', 'LineWidth',2)

xlim([1 trials])
ylim([0 0.075])
legend([e1.mainLine,e2.mainLine],'MD_{12}','MD_{21}','FontName', 'Times New Roman','Location','northwest')
legend boxoff
print(gcf,fullFileName,'-depsc', '-r300')

name = sprintf('md_VH');
fullFileName = fullfile(folder, name);
figure
set(gcf,'pos',[100 100 230 300])
patch([12 132 132 12],[0.0005 0.0005 0.075 0.075],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
%  ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman','Color', 'w')
% ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('VH', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
% shadedErrorBar([1:trials],nanmean(alltsvals12(subVH,:),1),nanstd(alltsvals12(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
% hold on
% shadedErrorBar([1:trials],nanmean(alltsvals21(subVH,:),1),nanstd(alltsvals21(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subVH,:),1)),smooth(nanstd(alltsvals12(subVH,:),1))./sqrt(length(subVH)),{'-','Color',[.1 .4 .9],'LineWidth',1},1);
hold on
shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subVH,:),1)),smooth(nanstd(alltsvals21(subVH,:),1))./sqrt(length(subVH)),{'-','Color',[.9 .1 .4] 'LineWidth',1},1);
xlim([1 trials])
ylim([0 0.075])

% pos = get(gca, 'Position'); % space adustment without ylabel
% pos(1) = 0.15;
% pos(3) = 0.78;
% set(gca, 'Position', pos)
print(gcf,fullFileName,'-depsc', '-r300')

name = sprintf('md_PV');
fullFileName = fullfile(folder, name);
figure
set(gcf,'pos',[100 100 230 300])
patch([12 132 132 12],[0.0005 0.0005 0.075 0.075],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
%  ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman','Color', 'w')
% ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('PV', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
% shadedErrorBar([1:trials],nanmean(alltsvals12(subVH,:),1),nanstd(alltsvals12(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
% hold on
% shadedErrorBar([1:trials],nanmean(alltsvals21(subVH,:),1),nanstd(alltsvals21(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subPV,:),1)),smooth(nanstd(alltsvals12(subPV,:)),1)./sqrt(length(subPV)),{'-','Color',[.1 .4 .9],'LineWidth',1},1);
hold on
shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subPV,:),1)),smooth(nanstd(alltsvals21(subPV,:)),1)./sqrt(length(subPV)),{'-','Color',[.9 .1 .4] 'LineWidth',1},1);
xlim([1 trials])
ylim([0 0.075])

% pos = get(gca, 'Position'); % space adustment without ylabel
% pos(1) = 0.15;
% pos(3) = 0.78;
% set(gca, 'Position', pos)
print(gcf,fullFileName,'-depsc', '-r300')


%% 3 groups same plot figure per subjects
%fullFileName = fullfile(folder, name);
c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.0005 0.0005 0.075 0.075],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('MD_{12} [m]', 'fontsize', 14,'FontName', 'Times New Roman')
% ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('Subject 1', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
% shadedErrorBar([1:trials],nanmean(alltsvals12(subVH,:),1),nanstd(alltsvals12(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
% hold on
% shadedErrorBar([1:trials],nanmean(alltsvals21(subVH,:),1),nanstd(alltsvals21(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
 e1=shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subH,:),1),10),smooth(nanstd(alltsvals12(subH,:),1),10)./sqrt(length(subH)),{'-','Color',c1,'LineWidth',1.5},1);
% e1=shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subH,:),1)),smooth(nanstd(alltsvals12(subH,:),1))./sqrt(length(subH)),{'LineStyle','none','Color',c1,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);

hold on
 %e2=shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subVH,:),1)),smooth(nanstd(alltsvals12(subVH,:),1))./sqrt(length(subVH)),{'LineStyle','none','Color',c2,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
  e2=shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subVH,:),1),10),smooth(nanstd(alltsvals12(subVH,:),1),10)./sqrt(length(subVH)),{'-','Color',c2, 'LineWidth',1.5},1);
hold on
 e3=shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subPV,:),1),10),smooth(nanstd(alltsvals12(subPV,:),1),10)./sqrt(length(subPV)),{'-','Color',c3, 'LineWidth',1.5},1);
 %e3=shadedErrorBar([1:trials],smooth(nanmean(alltsvals12(subPV,:),1)),smooth(nanstd(alltsvals12(subPV,:),1))./sqrt(length(subPV)),{'LineStyle','none','Color',c3,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
h=legend([e1.mainLine,e2.mainLine,e3.mainLine],'H','VH','PV','FontName', 'Times New Roman')

rect = [0.63, 0.75, .1, .1]; % to control position of legend 
set(h, 'Position', rect)
legend boxoff

xlim([1 trials])
ylim([0 0.06])
name = sprintf('md_all_sub1')
fullFileName = fullfile(folder, name)
print(gcf,fullFileName,'-depsc', '-r300')

figure
set(gcf,'pos',[100 100 230 300])
patch([12 132 132 12],[0.0005 0.0005 0.075 0.075],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
%ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
 ylabel('MD_{21} [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('Subject 2', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
% shadedErrorBar([1:trials],nanmean(alltsvals12(subVH,:),1),nanstd(alltsvals12(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
% hold on
% shadedErrorBar([1:trials],nanmean(alltsvals21(subVH,:),1),nanstd(alltsvals21(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
   e1=shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subH,:),1),10),smooth(nanstd(alltsvals21(subH,:),1),10)./sqrt(length(subH)),{'-','Color',c1,'LineWidth',1.5},1);
% e1=shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subH,:),1)),smooth(nanstd(alltsvals21(subH,:),1))./sqrt(length(subH)),{'LineStyle','none','Color',c1,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
hold on
 e2=shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subVH,:),1),10),smooth(nanstd(alltsvals21(subVH,:),1),10)./sqrt(length(subVH)),{'-','Color',c2, 'LineWidth',1.5},1);
% e2=shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subVH,:),1)),smooth(nanstd(alltsvals21(subVH,:),1))./sqrt(length(subVH)),{'LineStyle','none','Color',c2,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
hold on
  e3=shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subPV,:),1),10),smooth(nanstd(alltsvals21(subPV,:),1),10)./sqrt(length(subPV)),{'-','Color',c3, 'LineWidth',1.5},1);
% e3=shadedErrorBar([1:trials],smooth(nanmean(alltsvals21(subPV,:),1)),smooth(nanstd(alltsvals21(subPV,:),1))./sqrt(length(subPV)),{'LineStyle','none','Color',c3,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
% legend([e1.mainLine,e2.mainLine,e3.mainLine],'H','VH','PV','FontName', 'Times New Roman','Location','northwest')
% legend boxoff
set(gca, 'ytick', [])
xlim([1 trials])
ylim([0 0.06])
name = sprintf('md_all_sub2')
fullFileName = fullfile(folder, name)
print(gcf,fullFileName,'-depsc', '-r300')

%% Catch trial plot
md11_ts_cc(find(forces_c==1))=nan;
md22_ts_cc(find(forces_c==1))=nan;
md12_ts_cc(find(forces_c==1))=nan;
md21_ts_cc(find(forces_c==1))=nan;


md11_ts_cc(find(forces(:,[1:12, 133:156])==0))=nan;
md22_ts_cc(find(forces(:,[1:12, 133:156])==0))=nan;
md12_ts_cc(find(forces(:,[1:12, 133:156])==0))=nan;
md21_ts_cc(find(forces(:,[1:12, 133:156])==0))=nan;
% md11_ts_c(617:619,:)=[];
% md12_ts_c(617:619,:)=[];
% md21_ts_c(617:619,:)=[];
% md22_ts_c(617:619,:)=[];
% md11_ts_c = reshape(md11_ts_c,11,156);
% md22_ts_c = reshape(md22_ts_c,11,156);
% md12_ts_c = reshape(md12_ts_c,11,156);
% md21_ts_c = reshape(md21_ts_c,11,156);
% c1 = smooth(md12_ts_cc,20);
% c2 = smooth(md21_ts_cc,20);
% % c11 = smooth(md12_ts_c);
% % c22 = smooth(md21_ts_c);
% md12_ts_ccc = reshape(c1,11,156);
% md21_ts_ccc = reshape(c2,11,156);
% md12_ts_c = reshape(c11,11,56);
% md21_ts_c = reshape(c22,11,56);


%% plot
figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.0005 0.0005 0.075 0.075],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
hold on
% e1=shadedErrorBar([1:trials],smooth(nanmean(md12_ts_cc(subH,:),1),5),nanstd(md12_ts_cc(subH,:),1)./sqrt(length(subH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
e1=shadedErrorBar([1:trials],nanmean(md12_ts_cc(subH,:),1),nanstd(md12_ts_cc(subH,:),1)./sqrt(length(subH)),{'-','Color',[.1 .4 .9], 'LineWidth',1},1);
hold on
% e2=shadedErrorBar([1:trials],smooth(nanmean(md21_ts_cc(subH,:),1),5),nanstd(md21_ts_cc(subH,:),1)./sqrt(length(subH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
e2=shadedErrorBar([1:trials],nanmean(md21_ts_cc(subH,:),1),nanstd(md21_ts_cc(subH,:),1)./sqrt(length(subH)),{'-','Color',[.9 .1 .4],'LineWidth',1},1);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
%  ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman','Color', 'w')
ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('H', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
legend([e1.mainLine,e2.mainLine],'MD_{12}','MD_{21}','FontName', 'Times New Roman','Location','northwest')
legend boxoff
 xlim([1 trials])
ylim([0 0.075])
name = sprintf('catch_H');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

figure
set(gcf,'pos',[100 100 230 300])
patch([12 132 132 12],[0.0005 0.0005 0.075 0.075],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
hold on
%shadedErrorBar([1:trials],smooth(nanmean(md12_ts_cc(subVH,:),1),5),nanstd(md12_ts_cc(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
shadedErrorBar([1:trials],nanmean(md12_ts_cc(subVH,:),1),nanstd(md12_ts_cc(subVH,:),1)./sqrt(length(subVH)),{'-','Color',[.1 .4 .9], 'LineWidth',1},1);
hold on
%shadedErrorBar([1:trials],smooth(nanmean(md21_ts_cc(subVH,:),1),5),nanstd(md21_ts_cc(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
shadedErrorBar([1:trials],nanmean(md21_ts_cc(subVH,:),1),nanstd(md21_ts_cc(subVH,:),1)./sqrt(length(subVH)),{'-','Color',[.9 .1 .4],'LineWidth',1},1);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
%  ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman','Color', 'w')
% ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('VH', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
xlim([1 trials])
ylim([0 0.075])
% pos = get(gca, 'Position'); % space adustment without ylabel
% pos(1) = 0.15;
% pos(3) = 0.78;
% set(gca, 'Position', pos)
% set(gca, 'YTick', [])
name = sprintf('catch_VH');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')
%% avg pre post baseline and aftereffects
name = sprintf('md_stat');
fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
% title('Partner 2')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.3 -0.15 0 0.15 0.3], mean(avgallvalsij(subVH,[1 2 6 11 12])), ...
                         std(avgallvalsij(subVH,[1 2 6 11 12]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.3 -0.15 0 0.15 0.3], mean(avgallvalsij(subH,[1 2 6 11 12])), ...
                         std(avgallvalsij(subH,[1 2 6 11 12]))./sqrt(length(subH)),'k.')

h=bar([mean(avgallvalsij(subH,1)) mean(avgallvalsij(subH,2)) mean(avgallvalsij(subH,6)) mean(avgallvalsij(subH,11)) mean(avgallvalsij(subH,12)); ...
     mean(avgallvalsij(subVH,1)) mean(avgallvalsij(subVH,2)) mean(avgallvalsij(subVH,6)) mean(avgallvalsij(subVH,11)) mean(avgallvalsij(subVH,12))])
 
 
sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1 , 2]},[0.0005, 0.0005, 0.05])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Mean of MD_{12} and MD_{21}[m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.075])
 legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
set(h(2),'facecol','w')
set(h(3),'facecol',[0.9 0.9 0.9])
set(h(4),'facecol',[0.4 0.4 0.4])
set(h(1),'facecol',[0.4 0.9 0.4])
set(h(5),'facecol',[.7 .7 .1])
print(gcf,fullFileName,'-depsc', '-r300')

%% avg pre post baseline and aftereffects md12
% name = sprintf('md12_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
% title('Partner 2')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.3 -0.15 0 0.15 0.3], mean(allvals12(subVH,[1 2 6 11 12])), ...
                         std(allvals12(subVH,[1 2 6 11 12]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.3 -0.15 0 0.15 0.3], mean(allvals12(subH,[1 2 6 11 12])), ...
                         std(allvals12(subH,[1 2 6 11 12]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals12(subH,1)) mean(allvals12(subH,2)) mean(allvals12(subH,6)) mean(allvals12(subH,11)) mean(allvals12(subH,12)); ...
     mean(allvals12(subVH,1)) mean(allvals12(subVH,2)) mean(allvals12(subVH,6)) mean(allvals12(subVH,11)) mean(allvals12(subVH,12))])
 
 
sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15]},[0.0005, 0.0005])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
title('MD_{12}', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight', 'Normal')
ylabel('Distance[m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.075])
 legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
set(h(2),'facecol','w')
set(h(3),'facecol',[0.9 0.9 0.9])
set(h(4),'facecol',[0.4 0.4 0.4])
set(h(1),'facecol',[0.4 0.9 0.4])
set(h(5),'facecol',[.7 .7 .1])
% print(gcf,fullFileName,'-depsc', '-r300')

% avg pre post baseline and aftereffects md12
% name = sprintf('md21_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
% title('Partner 2')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.3 -0.15 0 0.15 0.3], mean(allvals21(subVH,[1 2 6 11 12])), ...
                         std(allvals21(subVH,[1 2 6 11 12]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.3 -0.15 0 0.15 0.3], mean(allvals21(subH,[1 2 6 11 12])), ...
                         std(allvals21(subH,[1 2 6 11 12]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals21(subH,1)) mean(allvals21(subH,2)) mean(allvals21(subH,6)) mean(allvals21(subH,11)) mean(allvals21(subH,12)); ...
     mean(allvals21(subVH,1)) mean(allvals21(subVH,2)) mean(allvals21(subVH,6)) mean(allvals21(subVH,11)) mean(allvals21(subVH,12))])
 
 
sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15]},[0.0005, 0.0005])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
title('MD_{21}', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight', 'Normal')
% ylabel('Distance[m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.075])
% legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
set(h(2),'facecol','w')
set(h(3),'facecol',[0.9 0.9 0.9])
set(h(4),'facecol',[0.4 0.4 0.4])
set(h(1),'facecol',[0.4 0.9 0.4])
set(h(5),'facecol',[.7 .7 .1])

pos = get(gca, 'Position'); % space adustment without ylabel
pos(1) = 0.15;
pos(3) = 0.78;
set(gca, 'Position', pos)
% print(gcf,fullFileName,'-depsc', '-r300')

%% epoch-wise shaded error bar

figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0 0 0.05 0.05],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
 xlabel('Epochs', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('VH group','FontName', 'Times New Roman')
hold on

shadedErrorBar([1:tsets], mean(allvals21(subVH,:)), std(allvals21(subVH,:))./sqrt(length(subVH)),{'-','Color',[.9 .1 .4],'LineWidth',1},1)
hold on
shadedErrorBar([1:tsets], mean(allvals12(subVH,:)), std(allvals12(subVH,:))./sqrt(length(subVH)),{'-','Color',[.1 .4 .9],'LineWidth',1},1)


% catch trials:
e3=shadedErrorBar(1:tsets,nanmean(md21_statmat_nofor(subVH,2:end)),...
    nanstd(md21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),{':','Color',[.9 .1 .4],'LineWidth',1},1)
e4=shadedErrorBar(1:tsets,nanmean(md12_statmat_nofor(subVH,2:end)),...
    nanstd(md12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)),{':','Color',[.1 .4 .9],'LineWidth',1},1)

%l1=legend([e2,e5,e3,e4],'md_{21}','md_{12}','catch md_{21}','catch md_{12}')
%set(l1,'Location','northeast')
legend boxoff
ylim([0.005 0.05])
%%Plot average, per-epoch (12,21)
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0 0 0.05 0.05],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
 xlabel('Epochs', 'fontsize', 14,'FontName', 'Times New Roman')
 ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('H group','FontName', 'Times New Roman')
hold on

shadedErrorBar(1:tsets, mean(allvals21(subH,:),1), std(allvals21(subH,:),1)./sqrt(length(subH)),{'-','Color',[.9 .1 .4],'LineWidth',1},1)
hold on
shadedErrorBar(1:tsets, mean(allvals12(subH,:),1), std(allvals12(subH,:),1)./sqrt(length(subH)),{'-','Color',[.1 .4 .9],'Linewidth',1},1)
legend boxoff

% catch trials:
% e(3)=shadedErrorBar(1:tsets,nanmean(md21_statmat_nofor(subH,2:end)),...
%     nanstd(md21_statmat_nofor(subH,2:end))./sqrt(length(subH)),'r:','LineWidth',1)
% e(4)=shadedErrorBar(1:tsets,nanmean(md12_statmat_nofor(subH,2:end)),...
%     nanstd(md12_statmat_nofor(subH,2:end))./sqrt(length(subH)),'b:','LineWidth',1)
% legend(e,{'Partner 2','Partner 1','Catch 2','Catch 1'})
% legend([e(2), e(1)],{'Subject 1','Subject 2'},'FontName', 'Times New Roman')
% legend boxoff
ylim([0.005 0.05])
%% subject-wise fitting

% alltsvals_corrected12 = alltsvals12;
% alltsvals_corrected21 = alltsvals21;
% alltsvals_corrected12(isnan(alltsvals_corrected12))=0;
% alltsvals_corrected21(isnan(alltsvals_corrected21))=0;
% coeffs_exp2 = zeros(11,4);
% coeffs_exp1 = zeros(11,2);
% coeffs_expfit12 = zeros(11,1);
% coeffs_expfit21 = zeros(11,1);
% for subj=1:subjs
% 
%     fit_per_sub = fit((13:132)', (alltsvals_corrected12(subj,13:132))','exp2');
%     coeffs_exp2(subj,:) = coeffvalues(fit_per_sub);
%     %coeffs_exp1(subj,:) = coeffvalues(fit_per_sub);
%     coeffs_expfit12(subj,:) = expfit(alltsvals_corrected12(subj,:));
%     coeffs_expfit21(subj,:) = expfit(alltsvals_corrected21(subj,:));
%     figure
%     set(gcf,'pos',[100 100 300 400])
% %     patch([12.5 132.5 132.5 12.5],[0 0 4 4],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% %     hold on;
%     xlabel('trials')
%     ylabel(plotlabel)
%     title(sprintf('Dyad%d',subj))
% %   plot(alltsvals_corrected(subj,:),'bo')
%     plot(alltsvals12(subj,:),'bo')
%     hold on
%     plot(fit_per_sub,'b-')
%     xlabel('trials')
%     ylabel('effort[N]')
%     xlim([1 trials])
%     ylim([0 0.07])
% 
% end
% 
% 
% %% error bar plot for fitting coefficients
% 
% figname = 'md12_fitting'
% figure
% set(gcf,'pos',[100 100 250 350])
% %title('fitting coefficents')
% hold on
% errorbar(2, mean(coeffs_expfit12(subVH,:)), ...
%                          std(coeffs_expfit12(subVH,:))./sqrt(length(subVH)),'k.')
% hold on
% errorbar(1, mean(coeffs_expfit12(subH,:)), ...
%                          std(coeffs_expfit12(subH,:))./sqrt(length(subH)),'k.')
%  
% h=bar([mean(coeffs_expfit12(subH,:)) ; ...
% mean(coeffs_expfit12(subVH,:))])
% % sigstar({[1,2]},[0.05])
% box off
% set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
% xlabel('group')
% ylim(yrange)
% ylabel('mean of fitting coefficient')
% % set(h(1),'facecol','w')
% set(h,'facecol',[0.8 0.8 0.8])
% saveas(gcf,[figdir,figname],'fig');
% eval(['print -depsc ', epsdir,figname])
% 
% 
% figname = 'md21_fitting'
% figure
% set(gcf,'pos',[100 100 250 350])
% %title('fitting coefficents')
% hold on
% errorbar(2, mean(coeffs_expfit21(subVH,:)), ...
%                          std(coeffs_expfit21(subVH,:))./sqrt(length(subVH)),'k.')
% hold on
% errorbar(1, mean(coeffs_expfit21(subH,:)), ...
%                          std(coeffs_expfit21(subH,:))./sqrt(length(subH)),'k.')
%  
% h=bar([mean(coeffs_expfit21(subH,:)) ; ...
% mean(coeffs_expfit21(subVH,:))]);
% % sigstar({[1,2]},[0.05])
% box off
% set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
% xlabel('group')
% ylim(yrange)
% ylabel('mean of fitting coefficient')
% % set(h(1),'facecol','w')
% set(h,'facecol',[0.8 0.8 0.8])
%  saveas(gcf,[figdir,figname],'fig');
%  eval(['print -depsc ', epsdir,figname])
% 
% %% Stats analysis 
% group = zeros(subjs,1);
% StatMatrix3 = [group, coeffs_expfit12,coeffs_expfit21];
% StatMatrix3(6,:)=[]; % take out subject 6
% StatMatrix3(subH,1)=1;
% StatMatrix3(6:end,1)=2;
% save('md_fitting_coeffs.txt','StatMatrix3', '-ascii')
% 
% 


%% trial ststistics...
alltsvals12_1 = nanmean(alltsvals12(:, 1:12),2)
alltsvals12_2 = nanmean(alltsvals12(:, 13:24),2)
alltsvals12_6 = nanmean(alltsvals12(:, 60:72),2)
alltsvals12_11 = nanmean(alltsvals12(:, 120:132),2)
alltsvals12_12 = nanmean(alltsvals12(:, 133:156),2)
sub1stat = [ alltsvals12_2 alltsvals12_6 alltsvals12_11]
alltsvals21_1 = nanmean(alltsvals21(:, 1:12),2)
alltsvals21_2 = nanmean(alltsvals21(:, 13:24),2)
alltsvals21_6 = nanmean(alltsvals21(:, 60:72),2)
alltsvals21_11 = nanmean(alltsvals21(:, 120:132),2)
alltsvals21_12 = nanmean(alltsvals21(:, 133:156),2)
sub2stat = [ alltsvals21_2 alltsvals21_6 alltsvals21_11 ]

%% figure

% name = sprintf('score1_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
title('MD_{12}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(sub1stat(subVH,:)), ...
                         std(sub1stat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(sub1stat(subH,:)), ...
                         std(sub1stat(subH,:))./sqrt(length(subH)),'k.')

h=bar([ mean(sub1stat(subH,1)) mean(sub1stat(subH,2)) mean(sub1stat(subH,3)) ; ...
      mean(sub1stat(subVH,1)) mean(sub1stat(subVH,2)) mean(sub1stat(subVH,3)) ])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225]},[0.0005, 0.0005])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.06])
%  legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])

name = sprintf('md12_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')


figure
%title('md_{21}')
% title('Partner 2')
set(gcf,'pos',[100 100 230 300])
errorbar(2+[-0.225 0 0.225], mean(sub2stat(subVH,:)), ...
                         std(sub2stat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(sub2stat(subH,:)), ...
                         std(sub2stat(subH,:))./sqrt(length(subH)),'k.')

h=bar([ mean(sub2stat(subH,1)) mean(sub2stat(subH,2)) mean(sub2stat(subH,3)) ; ...
      mean(sub2stat(subVH,1)) mean(sub2stat(subVH,2)) mean(sub2stat(subVH,3)) ])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1+0.225,2+0.225]},[0.0005, 0.0005, 0.05])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
set(gca, 'ytick', [])
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
% ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
title('MD_{21}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
 ylim([0 0.06])
%  legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])

set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])

% pos = get(gca, 'Position'); % space adustment without ylabel
% pos(1) = 0.15;
% pos(3) = 0.78;
% set(gca, 'Position', pos)

name = sprintf('md21_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%% figure 2


% name = sprintf('score1_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
title('MD_{12}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(sub1stat(subVH,:)), ...
                         std(sub1stat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(sub1stat(subH,:)), ...
                         std(sub1stat(subH,:))./sqrt(length(subH)),'k.')
hold on
errorbar(3+[-0.225 0 0.225], mean(sub1stat(subPV,:)), ...
                         std(sub1stat(subPV,:))./sqrt(length(subPV)),'k.')

h=bar([mean(sub1stat(subH,1)) mean(sub1stat(subH,2)) mean(sub1stat(subH,3)) ; ...
      mean(sub1stat(subVH,1)) mean(sub1stat(subVH,2)) mean(sub1stat(subVH,3));
      mean(sub1stat(subPV,1)) mean(sub1stat(subPV,2)) mean(sub1stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225], [3-0.225,3+0.225]},[0.0005, 0.0005, 0.0005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.06])
%  legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])

name = sprintf('md12_all_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')


figure
%title('md_{21}')
% title('Partner 2')
set(gcf,'pos',[100 100 230 300])
errorbar(2+[-0.225 0 0.225], mean(sub2stat(subVH,:)), ...
                         std(sub2stat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(sub2stat(subH,:)), ...
                         std(sub2stat(subH,:))./sqrt(length(subH)),'k.')
hold on
errorbar(3+[-0.225 0 0.225], mean(sub2stat(subPV,:)), ...
                         std(sub2stat(subPV,:))./sqrt(length(subPV)),'k.')

h=bar([ mean(sub2stat(subH,1)) mean(sub2stat(subH,2)) mean(sub2stat(subH,3)) ; ...
      mean(sub2stat(subVH,1)) mean(sub2stat(subVH,2)) mean(sub2stat(subVH,3));
       mean(sub2stat(subPV,1)) mean(sub2stat(subPV,2)) mean(sub2stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1+0.225,2+0.225],[3-0.225,3+0.225], [1+0.225,3+0.225]},[0.0005, 0.0005, 0.05, 0.0005, 0.005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'},'FontName', 'Times New Roman')
set(gca, 'ytick', [])
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
% ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
title('MD_{21}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
title('Subject 2','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
 ylim([0 0.06])
%  legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])

set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])

% pos = get(gca, 'Position'); % space adustment without ylabel
% pos(1) = 0.15;
% pos(3) = 0.78;
% set(gca, 'Position', pos)

name = sprintf('md21_all_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')
%% figure 3 

c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
% name = sprintf('score1_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
title('MD_{12}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(1+[-0.225 0 0.225], [mean(sub1stat(subH,1)) mean(sub1stat(subVH,1)) mean(sub1stat(subPV,1))], ...
                         [std(sub1stat(subH,1))./sqrt(length(subH)) std(sub1stat(subVH,1))./sqrt(length(subVH)) std(sub1stat(subPV,1))./sqrt(length(subPV))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [mean(sub1stat(subH,2)) mean(sub1stat(subVH,2)) mean(sub1stat(subPV,2))], ...
                         [std(sub1stat(subH,2))./sqrt(length(subH)) std(sub1stat(subVH,2))./sqrt(length(subVH)) std(sub1stat(subPV,2))./sqrt(length(subPV))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [mean(sub1stat(subH,3)) mean(sub1stat(subVH,3)) mean(sub1stat(subPV,3))], ...
                         [std(sub1stat(subH,3))./sqrt(length(subH)) std(sub1stat(subVH,3))./sqrt(length(subVH)) std(sub1stat(subPV,3))./sqrt(length(subPV))],'k.')
 hold on
h=bar([mean(sub1stat(subH,1)) mean(sub1stat(subVH,1)) mean(sub1stat(subPV,1)) ; ...
      mean(sub1stat(subH,2)) mean(sub1stat(subVH,2)) mean(sub1stat(subPV,2));
      mean(sub1stat(subH,3)) mean(sub1stat(subVH,3)) mean(sub1stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
  sigstar({[2,2+0.225],[2-0.225,2+0.225], [3-0.225,3+0.225]},[0.05, 0.0005, 0.05])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'},'FontName', 'Times New Roman')
xlabel('Epoch', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('MD_{12} [m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.06])
%  legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)

name = sprintf('md12_all_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')


figure
%title('md_{21}')
title('MD_{12}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
title('Subject 2','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(1+[-0.225 0 0.225], [mean(sub2stat(subH,1)) mean(sub2stat(subVH,1)) mean(sub2stat(subPV,1))], ...
                         [std(sub2stat(subH,1))./sqrt(length(subH)) std(sub2stat(subVH,1))./sqrt(length(subVH)) std(sub2stat(subPV,1))./sqrt(length(subPV))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [mean(sub2stat(subH,2)) mean(sub2stat(subVH,2)) mean(sub2stat(subPV,2))], ...
                         [std(sub2stat(subH,2))./sqrt(length(subH)) std(sub2stat(subVH,2))./sqrt(length(subVH)) std(sub2stat(subPV,2))./sqrt(length(subPV))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [mean(sub2stat(subH,3)) mean(sub2stat(subVH,3)) mean(sub2stat(subPV,3))], ...
                         [std(sub2stat(subH,3))./sqrt(length(subH)) std(sub2stat(subVH,3))./sqrt(length(subVH)) std(sub2stat(subPV,3))./sqrt(length(subPV))],'k.')
 hold on
h=bar([mean(sub2stat(subH,1)) mean(sub2stat(subVH,1)) mean(sub2stat(subPV,1)) ; ...
      mean(sub2stat(subH,2)) mean(sub2stat(subVH,2)) mean(sub2stat(subPV,2));
      mean(sub2stat(subH,3)) mean(sub2stat(subVH,3)) mean(sub2stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
  sigstar({[2-0.225,2+0.225],[3-0.225,3+0.225], [3,3-0.225]},[0.05, 0.0005, 0.0005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'},'FontName', 'Times New Roman')
xlabel('Epoch', 'fontsize', 14,'FontName', 'Times New Roman')
 ylabel('MD_{21} [m]', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 0.06])
%  legend(h,{'Baseline','Early training','Middle training','Late training', 'Washout'}, 'Location', 'northwest')
legend boxoff
set(gca, 'ytick', [])
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)

name = sprintf('md21_all_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%% Stat values
S_1 = [sub1stat]; 
S_1(6,:) = [];
group_S = [md22_statmat(:,1)];
group_S(6,:)=[];
subjects_S_1 = (1:length(S_1(:,1)))';
pre_S_1 = [ S_1(:,1) group_S ones(length(group_S),1) subjects_S_1];
mid_S_1 = [ S_1(:,2) group_S 2*ones(length(group_S),1) subjects_S_1];
pos_S_1 = [ S_1(:,3) group_S 3*ones(length(group_S),1) subjects_S_1];
A_S_1 = [pre_S_1; mid_S_1; pos_S_1];

name = sprintf('md12_three_group_fourcolumn.txt');
fullFileName = fullfile(stat_folder, name);
save(fullFileName,'A_S_1', '-ascii')
BWAOV2(A_S_1);
%
S_2 = [sub2stat]; 
S_2(6,:) = [];
subjects_S_2 = (1:length(S_2(:,1)))';
pre_S_2 = [ S_2(:,1) group_S ones(length(group_S),1) subjects_S_2];
mid_S_2 = [ S_2(:,2) group_S 2*ones(length(group_S),1) subjects_S_2];
pos_S_2 = [ S_2(:,3) group_S 3*ones(length(group_S),1) subjects_S_2];
A_S_2 = [pre_S_2; mid_S_2; pos_S_2];
name = sprintf('md21_three_group_fourcolumn.txt');
fullFileName = fullfile(stat_folder, name);
save(fullFileName,'A_S_2', '-ascii')
BWAOV2(A_S_2);
