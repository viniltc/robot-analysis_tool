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
subPV = 12:16;  

load([resdir, 'score1_statmat.txt']);
load([resdir, 'score1_statmat_nofor.txt']);
load([resdir, 'score2_statmat.txt']);
load([resdir, 'score2_statmat_nofor.txt']);

load([resdir, 'md11_statmat.txt']);
load([resdir, 'md11_statmat_nofor.txt']);
load([resdir, 'md12_statmat.txt']);
load([resdir, 'md12_statmat_nofor.txt']);
load([resdir, 'md21_statmat.txt']);
load([resdir, 'md21_statmat_nofor.txt']);
load([resdir, 'md22_statmat.txt']);
load([resdir, 'md22_statmat_nofor.txt']);

group = md11_statmat(:,1);

% load per-trial data
load ([resdir,'score1.mat']);
score1_ts = (ind_ts);
load ([resdir,'score2.mat']);
score2_ts = (ind_ts);
load([resdir, 'forces.mat']);
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


%catch trial
md11_ts_c=md11_ts(find(forces==0));
md22_ts_c=md22_ts(find(forces==0));
md12_ts_c=md12_ts(find(forces==0));
md21_ts_c=md21_ts(find(forces==0));
score1_ts_c=md11_ts(find(forces==0));
score2_ts_c=md22_ts(find(forces==0));

md11 = md11_statmat(:,2:end);
md12 = md12_statmat(:,2:end);
md21 = md21_statmat(:,2:end);
md22 = md22_statmat(:,2:end);

score1 = score1_statmat(:,2:end);
score2 = score2_statmat(:,2:end);

% load per-trial score data

load([resdir, 'scores.mat']);
n_score = size(scores);
score_1d = (scores(2:3:n_score(1),120:132))'; %% all dyads, sub1, epoch 11
score_2d = (scores(3:3:n_score(1),120:132))';

score_1d = (scores(2:3:n_score(1),:)); %% all dyads , subject 1
score_2d = (scores(3:3:n_score(1),:)); %% all dyads , subject 2

% score_1dd = score_1d;
% score_2dd = score_2d;
% score_1dd(find(forces==0))=[];
% score_2dd(find(forces==0))=[];

score_1d(find(forces_c==0))=nan;
score_2d(find(forces_c==0))=nan;

% score_1dd = reshape(score_1dd, 11, 156);
% score_2dd = reshape(score_2dd, 11, 156);


% score_1d_c=score_1d(find(forces==0));
% score_2d_c=score_2d(find(forces==0));
% score_1d_c = reshape(score_1d_c, 11, 56);
% score_2d_c = reshape(score_2d_c, 11, 56);


% md11_ts_c = reshape(md11_ts_c,11,56);
% md22_ts_c = reshape(md22_ts_c,11,56);
% md12_ts_c = reshape(md12_ts_c,11,56);
% md21_ts_c = reshape(md21_ts_c,11,56);




% score_1d = mean(score_1d);
% score_2d = mean(score_2d);
score_dyads = [score_1d; score_2d];
% score_dyads = mean(score_dyads);

md_avg = 0.5*(md12+md21); % minimumdistance average
subjs=size(md11,1);
tsets=size(md22,2);
allvals11 = [];
allvals22 = [];
allvals12 = [];
allvals21 = [];
avgallvalsij = [];
scoreallvals1 = [];
scoreallvals2 =[];
for subj=1:subjs

     
    allvals11 = [allvals11; md11(subj,:)];
    allvals22 = [allvals22; md22(subj,:)];
    allvals12 = [allvals12; md12(subj,:)];
    allvals21 = [allvals21; md21(subj,:)];
    avgallvalsij = [avgallvalsij; md_avg(subj,:)];
    scoreallvals1 = [scoreallvals1; score1(subj,:)];
    scoreallvals2 = [scoreallvals2; score2(subj,:)];
    
    
     
end

%% Plot per-subject, per-trial
scoretsallvals1 = [];
scoretsallvals2 = [];
alltsvals12 = [];
alltsvals21 = [];
alltsvals11 = [];
alltsvals22 = [];
trials = size(md11_ts,2);
for subj=1:subjs
    
%    figure
%   set(gcf,'pos',[100 100 250 350])
%   patch([12 132 132 12],[0.0005 0.0005 110 110],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
%     xlabel('trials')
%      ylabel('Score')
%     title(sprintf('Dyad %d',subj))
%         line(1:trials,score_1d(subj,:),'col','b','marker','.','lines','-');
%         line(1:trials,score_2d(subj,:),'col','r','marker','s','lines','-');
%         
%     xlim([1 trials])
%     ylim([0 110])

    alltsvals11 = [alltsvals11; md11_ts(subj,:)];
    alltsvals22 = [alltsvals22; md22_ts(subj,:)];
    alltsvals12 = [alltsvals12; md12_ts(subj,:)];
    alltsvals21 = [alltsvals21; md21_ts(subj,:)];
%     scoretsallvals1 = [scoretsallvals1; score_1d(subj,:)];
%     scoretsallvals2 = [scoretsallvals2; score_2d(subj,:)];
    scoretsallvals1 = [scoretsallvals1; score_1d(subj,:)];
    scoretsallvals2 = [scoretsallvals2; score_2d(subj,:)];
    


%   c1 = smooth(scoretsallvals1, 20);
%   c2 = smooth(scoretsallvals2, 20);

%   c11 = smooth(scoretsallvals1);
%   c22 = smooth(scoretsallvals2);
end

%   c1 = smooth(scoretsallvals1,10);
%   c2 = smooth(scoretsallvals2,10);
%   c1 = smooth(scoretsallvals1,'sgolay');
%   c2 = smooth(scoretsallvals2, 'sgolay');
%   c1 = smooth(scoretsallvals1,'loess');
%   c2 = smooth(scoretsallvals2, 'loess');
  
%   
%   scoretsallvals11 = reshape(c1,11,156);
%   scoretsallvals22 = reshape(c2,11,156);
%   scoretsallvals11 = reshape(c1,11,100);
%   scoretsallvals22 = reshape(c2,11,100);
%   scoretsallvals111 = reshape(c11,11,156);
%   scoretsallvals222 = reshape(c22,11,156);
  
  
%% Score over trial shaded error bar
%trials = 100;

figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.0005 0.0005 110 110],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% patch([1 132 132 1],[0.0005 0.0005 110 110],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
title('H', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
e1=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subH,:),1),10,'sgolay'),smooth(nanstd(scoretsallvals1(subH,:),1),10,'sgolay')./sqrt(length(subH)),{'-','Color',[.1 .4 .9], 'LineWidth',1},1);
%e1=shadedErrorBar([1:trials],nanmean(scoretsallvals1(subVH,:),1),nanstd(scoretsallvals1(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
hold on
e2=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals2(subH,:),1),10,'sgolay'),smooth(nanstd(scoretsallvals2(subH,:),1),10,'sgolay')./sqrt(length(subH)),{'-','Color',[.9 .1 .4], 'LineWidth',1},1);
%e2=shadedErrorBar([1:trials],nanmean(scoretsallvals2(subVH,:),1),nanstd(alltsvals21(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);


xlim([1 trials])
ylim([0 110])
legend([e1.mainLine,e2.mainLine],'Subject 1','Subject 2','FontName', 'Times New Roman','Location','northwest')
legend boxoff
name = sprintf('score_H');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')


figure
set(gcf,'pos',[100 100 230 300])
 patch([12 132 132 12],[0.0005 0.0005 110 110],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% patch([1 132 132 1],[0.0005 0.0005 110 110],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
% ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman','Color', [1 1 1])
% ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('VH', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
% shadedErrorBar([1:trials],nanmean(alltsvals12(subVH,:),1),nanstd(alltsvals12(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
% hold on
% shadedErrorBar([1:trials],nanmean(alltsvals21(subVH,:),1),nanstd(alltsvals21(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subVH,:),1),10,'sgolay'),smooth(nanstd(scoretsallvals1(subVH,:),1),10,'sgolay')./sqrt(length(subVH)),{'-','Color',[.1 .4 .9],'LineWidth',1},1);
hold on
shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals2(subVH,:),1),10,'sgolay'),smooth(nanstd(scoretsallvals2(subVH,:),1),10,'sgolay')./sqrt(length(subVH)),{'-','Color',[.9 .1 .4] 'LineWidth',1},1);
xlim([1 trials])
ylim([0 110])
 
name = sprintf('score_VH');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

figure
set(gcf,'pos',[100 100 230 300])
 patch([12 132 132 12],[0.0005 0.0005 110 110],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% patch([1 132 132 1],[0.0005 0.0005 110 110],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
% ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman','Color', [1 1 1])
% ylabel('Minimum distance [m]', 'fontsize', 14,'FontName', 'Times New Roman')
title('PV', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
% shadedErrorBar([1:trials],nanmean(alltsvals12(subVH,:),1),nanstd(alltsvals12(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.1 .4 .9],'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
% hold on
% shadedErrorBar([1:trials],nanmean(alltsvals21(subVH,:),1),nanstd(alltsvals21(subVH,:),1)./sqrt(length(subVH)),{'LineStyle','none','Color',[.9 .1 .4],'Marker','.', 'MarkerSize',15,'LineWidth',1},1);
shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subPV,:),1),10,'sgolay'),smooth(nanstd(scoretsallvals1(subPV,:),1),10,'sgolay')./sqrt(length(subPV)),{'-','Color',[.1 .4 .9],'LineWidth',1},1);
hold on
shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals2(subPV,:),1),10,'sgolay'),smooth(nanstd(scoretsallvals2(subPV,:),1),10,'sgolay')./sqrt(length(subPV)),{'-','Color',[.9 .1 .4] 'LineWidth',1},1);
xlim([1 trials])
ylim([0 110])
 
name = sprintf('score_PV');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%% for stats

score1vals_1 = nanmean(scoretsallvals1(:, 1:12),2)
score1vals_2 = nanmean(scoretsallvals1(:, 13:24),2)
score1vals_6 = nanmean(scoretsallvals1(:, 60:72),2)
score1vals_11 = nanmean(scoretsallvals1(:, 120:132),2)
score1vals_12 = nanmean(scoretsallvals1(:, 133:156),2)
%score1stat = [score1vals_1 score1vals_2 score1vals_6 score1vals_11 score1vals_12]
score1stat = [score1vals_2 score1vals_6 score1vals_11]
score2vals_1 = nanmean(scoretsallvals2(:, 1:12),2)
score2vals_2 = nanmean(scoretsallvals2(:, 13:24),2)
score2vals_6 = nanmean(scoretsallvals2(:, 60:72),2)
score2vals_11 = nanmean(scoretsallvals2(:, 120:132),2)
score2vals_12 = nanmean(scoretsallvals2(:, 133:156),2)
%score2stat = [score2vals_1 score2vals_2 score2vals_6 score2vals_11 score2vals_12]
score2stat = [score2vals_2 score2vals_6 score2vals_11]

name = sprintf('score1_stat');
fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(score1stat(subVH,:)), ...
                         std(score1stat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(score1stat(subH,:)), ...
                         std(score1stat(subH,:))./sqrt(length(subH)),'k.')
 hold on
errorbar(3+[-0.225 0 0.225], mean(score1stat(subPV,:)), ...
                         std(score1stat(subPV,:))./sqrt(length(subPV)),'k.')

h=bar([ mean(score1stat(subH,1)) mean(score1stat(subH,2)) mean(score1stat(subH,3)) ; ...
      mean(score1stat(subVH,1)) mean(score1stat(subVH,2)) mean(score1stat(subVH,3));
      mean(score1stat(subPV,1)) mean(score1stat(subPV,2)) mean(score1stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1+0.225,2+0.225],[3+0.225, 2+0.225], [3+0.225, 1+0.225],[3+0.225, 3-0.225] },[0.0005, 0.0005, 0.005, 0.05, 0.005,0.005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH' , 'PV'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 130])
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
%print(gcf,fullFileName,'-depsc', '-r300')

name = sprintf('score2_stat');
fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
% title('Partner 2')
hold on
set(gcf,'pos',[100 100 230 300])
errorbar(2+[-0.225 0 0.225], mean(score2stat(subVH,:)), ...
                         std(score2stat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(score2stat(subH,:)), ...
                         std(score2stat(subH,:))./sqrt(length(subH)),'k.')
hold on
                     errorbar(3+[-0.225 0 0.225], mean(score2stat(subPV,:)), ...
                         std(score2stat(subPV,:))./sqrt(length(subPV)),'k.')

h=bar([ mean(score2stat(subH,1)) mean(score2stat(subH,2)) mean(score2stat(subH,3)) ; ...
      mean(score2stat(subVH,1)) mean(score2stat(subVH,2)) mean(score2stat(subVH,3));
       mean(score2stat(subPV,1)) mean(score2stat(subPV,2)) mean(score2stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1+0.225,2+0.225],[3+0.225, 2+0.225], [3+0.225, 1+0.225], [3-0.225,3+0.225]},[0.0005, 0.0005, 0.005, 0.05, 0.005, 0.005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH' 'PV'},'FontName', 'Times New Roman')
set(gca, 'ytick', []);
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
% ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
title('Subject 2','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
ylim([0 130])
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

%print(gcf,fullFileName,'-depsc', '-r300')



%% Three groups in 1
c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
window =5;
figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.0005 0.0005 100 100],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
% title('PV', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
%e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subPV,:),1),nanstd(alltsvalsint(subPV,:)./sqrt(length(subPV)),1),{'-k', 'LineWidth',1},1);
%e1=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subH,:),1),20,'sgolay'),smooth(nanstd(scoretsallvals1(subH,:),1),20,'sgolay')./sqrt(length(subH)),{'-','Color',c1,'LineWidth',1},1);
e1=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subH,:),1),window),smooth(nanstd(scoretsallvals1(subH,:),1),window)./sqrt(length(subH)),{'LineStyle','none','Color',c1,'Marker','.', 'MarkerSize',15 ,'LineWidth',1},1);
hold on
%e2=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subVH,:),1),20,'sgolay'),smooth(nanstd(scoretsallvals1(subVH,:),1),20,'sgolay')./sqrt(length(subVH)),{'-','Color',c2,'LineWidth',1},1);
e2=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subVH,:),1),window),smooth(nanstd(scoretsallvals1(subVH,:),1),window)./sqrt(length(subVH)),{'LineStyle','none','Color',c2,'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
hold on
%e3=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subPV,:),1),20,'sgolay'),smooth(nanstd(scoretsallvals1(subPV,:),1),20,'sgolay')./sqrt(length(subPV)),{'-','Color',c3,'LineWidth',1},1);
e3=shadedErrorBar([1:trials],smooth(nanmean(scoretsallvals1(subPV,:),1),window),smooth(nanstd(scoretsallvals1(subPV,:),1),window)./sqrt(length(subPV)),{'LineStyle','none','Color',c3,'Marker','.', 'MarkerSize',15, 'LineWidth',1},1);
plot([1 160],[100 100],'--','Color',[.01 .01 .01])
h=legend([e1.mainLine,e2.mainLine,e3.mainLine],'H','VH','PV','FontName', 'Times New Roman')
hold on
rect = [0.62, 0.25, .5, .1]; % to control position of legend 
set(h, 'Position', rect)
 legend boxoff
% set(h,'position',[0.62 0.25 0.2 0.2])
xlim([1 156])
ylim([0 100])
name = sprintf('score_all');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%%
figure
%title('md_{21}')
% title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(score1stat(subVH,:)), ...
                         std(score1stat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(score1stat(subH,:)), ...
                         std(score1stat(subH,:))./sqrt(length(subH)),'k.')
 hold on
errorbar(3+[-0.225 0 0.225], mean(score1stat(subPV,:)), ...
                         std(score1stat(subPV,:))./sqrt(length(subPV)),'k.')

h=bar([ mean(score1stat(subH,1)) mean(score1stat(subH,2)) mean(score1stat(subH,3)) ; ...
      mean(score1stat(subVH,1)) mean(score1stat(subVH,2)) mean(score1stat(subVH,3));
      mean(score1stat(subPV,1)) mean(score1stat(subPV,2)) mean(score1stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1+0.225,2+0.225],[3+0.225, 2+0.225], [3+0.225, 1+0.225],[3+0.225, 3-0.225] },[0.0005, 0.0005, 0.005, 0.05, 0.005,0.005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH' , 'PV'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 140])
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
% name = sprintf('score_all_stat');
% fullFileName = fullfile(folder, name);
% print(gcf,fullFileName,'-depsc', '-r300')

%% figure 3 

c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
% name = sprintf('score1_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
% title('MD_{12}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
% title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(1+[-0.225 0 0.225], [mean(score1stat(subH,1)) mean(score1stat(subVH,1)) mean(score1stat(subPV,1))], ...
                         [std(score1stat(subH,1))./sqrt(length(subH)) std(score1stat(subVH,1))./sqrt(length(subVH)) std(score1stat(subPV,1))./sqrt(length(subPV))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [mean(score1stat(subH,2)) mean(score1stat(subVH,2)) mean(score1stat(subPV,2))], ...
                         [std(score1stat(subH,2))./sqrt(length(subH)) std(score1stat(subVH,2))./sqrt(length(subVH)) std(score1stat(subPV,2))./sqrt(length(subPV))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [mean(score1stat(subH,3)) mean(score1stat(subVH,3)) mean(score1stat(subPV,3))], ...
                         [std(score1stat(subH,3))./sqrt(length(subH)) std(score1stat(subVH,3))./sqrt(length(subVH)) std(score1stat(subPV,3))./sqrt(length(subPV))],'k.')
 hold on
h=bar([mean(score1stat(subH,1)) mean(score1stat(subVH,1)) mean(score1stat(subPV,1)) ; ...
      mean(score1stat(subH,2)) mean(score1stat(subVH,2)) mean(score1stat(subPV,2));
      mean(score1stat(subH,3)) mean(score1stat(subVH,3)) mean(score1stat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
   sigstar({[1-0.225,1+0.225],[1,1+0.225],[2,2+0.225],[2-0.225,2+0.225],[3-0.225,3+0.225],[3-0.225,3]},[0.0005, 0.0005, 0.0005, 0.0005, 0.0005,0.005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'},'FontName', 'Times New Roman')
xlabel('Epoch', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Score', 'fontsize', 14,'FontName', 'Times New Roman')
ylim([0 140])
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

name = sprintf('score_all_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')





%% Stat values
S_1 = [score1stat]; 
S_1(6,:) = [];
group_S = [md22_statmat(:,1)];
group_S(6,:)=[];
subjects_S_1 = (1:length(S_1(:,1)))';
pre_S_1 = [ S_1(:,1) group_S ones(length(group_S),1) subjects_S_1];
mid_S_1 = [ S_1(:,2) group_S 2*ones(length(group_S),1) subjects_S_1];
pos_S_1 = [ S_1(:,3) group_S 3*ones(length(group_S),1) subjects_S_1];
A_S_1 = [pre_S_1; mid_S_1; pos_S_1];
name = sprintf('score_1_three_group_fourcolumn.txt');
fullFileName = fullfile(stat_folder, name);
save(fullFileName,'A_S_1', '-ascii')
BWAOV2(A_S_1);

%
S_2 = [score2stat]; 
S_2(6,:) = [];
subjects_S_2 = (1:length(S_2(:,1)))';
pre_S_2 = [ S_2(:,1) group_S ones(length(group_S),1) subjects_S_2];
mid_S_2 = [ S_2(:,2) group_S 2*ones(length(group_S),1) subjects_S_2];
pos_S_2 = [ S_2(:,3) group_S 3*ones(length(group_S),1) subjects_S_2];
A_S_2 = [pre_S_2; mid_S_2; pos_S_2];
name = sprintf('score_2_three_group_fourcolumn.txt');
fullFileName = fullfile(stat_folder, name);
save(fullFileName,'A_S_2', '-ascii')
BWAOV2(A_S_2);
