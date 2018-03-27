%% Plot effort figures

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
phases = {'baseline','training','aftereffect'}; 
order = {1,2:11,12:13};  % target sets in experiment 

subH = 1:5;
subVH = 7:11;  % S5 removed
subPV = 12:16;  % 


% load per-epoch effort data
load([resdir, 'eff1_statmat.txt']);
load([resdir, 'eff1_statmat_catch.txt']);
load([resdir, 'eff2_statmat.txt']);
load([resdir, 'eff2_statmat_catch.txt']);
load([resdir, 'R_statmat.txt']);
load([resdir, 'R_statmat_catch.txt']);

% load per-epoch speed data
load([resdir, 'speed1_statmat.txt']);
load([resdir, 'speed1_statmat_catch.txt']);
load([resdir, 'speed2_statmat.txt']);
load([resdir, 'speed2_statmat_catch.txt']);

group = eff1_statmat(:,1);

% load per-trial data
load ([resdir,'eff1.mat']);
eff1_ts = (ind_ts);
load ([resdir,'eff2.mat']);
eff2_ts = (ind_ts);
load([resdir, 'forces.mat']);
load ([resdir,'R.mat']);
R_ts = sqrt(ind_ts);

intforce_ts = 0.5*(eff1_ts+eff2_ts);

% load per-trial data
load ([resdir,'speed1.mat']);
speed1_ts = sqrt(ind_ts);
load ([resdir,'speed2.mat']);
speed2_ts = sqrt(ind_ts);
load([resdir, 'forces.mat']);

eff1_ts(find(forces==0))= nan;
eff2_ts(find(forces==0))= nan;
% forces_c = forces(:, [13:120]); %% catch trial in phase training
% 
% eff1_ts(find(forces_c==0))= mean(eff1_ts(12:13,:),2);
% eff2_ts(find(forces_c==0))= mean(eff2_ts(12:13,:),2);
%speed1_ts(find(forces==0))= nan;
%speed2_ts(find(forces==0))= nan;

plotlabel = 'RMS effort [N]';
yrange = [0 6]; % range of variation of effort


% %% Plot per-subject, per-epoch speed....
% speed1 = sqrt(speed1_statmat(:,2:end));
% speed2 = sqrt(speed2_statmat(:,2:end));
% R = sqrt(R_statmat(:,2:end));
% 
% 
% subjs=size(speed1,1);
% tsets=size(speed1,2);
% allvalsv1 = [];
% allvalsv2 = [];
% for subj=1:subjs
% %     figure
% %     set(gcf,'pos',[100 100 300 400])
% %     patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% %     xlabel('epochs')
% %     ylabel(plotlabel)
% %     title(sprintf('dyad%d',subj))
% %         line(1:tsets,speed1(subj,:),'col','r','marker','s');
% %         line(1:tsets,speed2(subj,:),'col','b','marker','s');
% %        % line(1:tsets,xcorr(speed1(subj,:),speed2(subj,:)),'col','g','marker','s');
% %     xlim([1 tsets])
%     %ylim(yrange)
%     
%     allvalsv1 = [allvalsv1; speed1(subj,:)];
%     allvalsv2 = [allvalsv2; speed2(subj,:)];
% 
% end


%% Plot per-subject, per-epoch
eff1 = (eff1_statmat(:,2:end));
eff2 = (eff2_statmat(:,2:end));
intforce=0.5*(eff1+eff2);

subjs=size(eff1,1);
tsets=size(eff1,2);
allvals1 = [];
allvals2 = [];
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel(plotlabel)
    title(sprintf('S%d',subj))
        line(1:tsets,eff1(subj,:),'col','r','marker','s');
        line(1:tsets,eff2(subj,:),'col','b','marker','s');
    xlim([1 tsets])
    ylim(yrange)
    
    allvals1 = [allvals1; eff1(subj,:)];
    allvals2 = [allvals2; eff2(subj,:)];
    allvalsint = 0.5*(allvals1+allvals2);
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
ylim(yrange)

%% Stats analysis (epoch wise)

group = zeros(subjs,1);
StatMatrix = [group, allvals1(:,[2 11]), allvals2(:,[2 11])];
StatMatrix(6,:)=[]; % take out subject 6
StatMatrix(subH,1)=1;
StatMatrix(6:end,1)=2;
save('effort_early_late_epoch.txt','StatMatrix', '-ascii')


%% Plot average, pre-post
figure
set(gcf,'pos',[100 100 250 350])
title('Subject 1')
hold on
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
ylabel(plotlabel)
ylim([0 4])
%ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])


figure
set(gcf,'pos',[100 100 250 350])
title('Subject 2')
hold on
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
ylabel(plotlabel)
ylim([0 4])
%ylim(yrange)
%legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])

 
%% Plot per-subject, per-trial
alltsvals1 = [];
alltsvals2 = [];
trials = size(eff1_ts,2);
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('trials')
    ylabel(plotlabel)
    title(sprintf('S%d',subj))
        line(1:trials,eff1_ts(subj,:),'col','r','marker','s');
        line(1:trials,eff2_ts(subj,:),'col','b','marker','s');
    xlim([1 trials])
    ylim(yrange)
    
    alltsvals1 = [alltsvals1; eff1_ts(subj,:)];
    alltsvals2 = [alltsvals2; eff2_ts(subj,:)];
    alltsvalsint = 0.5*(alltsvals1+alltsvals2);
end


%% Stats analysis (trial wise)

group = zeros(subjs,1);
StatMatrix1 = [group, alltsvals1(:,[(13:24) (121:132)]), alltsvals2(:,[(13:24) (121:132)])];
StatMatrix1(6,:)=[]; % take out subject 6
StatMatrix1(subH,1)=1;
StatMatrix1(6:end,1)=2;
save('effort_early_late_trial.txt','StatMatrix1', '-ascii')

 %% Plot average per trial 
figure
set(gcf,'pos',[100 100 300 400])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel('RMS interaction force [N]')
title('group VH')
hold on
% errorbar(1:trials, nanmean(alltsvals1(subVH,:)), nanstd(alltsvals1(subVH,:))./sqrt(sum(isfinite(alltsvals1(subVH,:)))),'rs-')
% hold on
errorbar(1:trials, nanmean(alltsvals2(subVH,:)), nanstd(alltsvals2(subVH,:))./sqrt(sum(isfinite(alltsvals1(subVH,:)))),'bs-')
xlim([1 trials])
ylim(yrange)

figure
set(gcf,'pos',[100 100 300 400])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel('RMS interaction force [N]')
title('group H')
hold on
% errorbar(1:trials, nanmean(alltsvals1(subH,:)), nanstd(alltsvals1(subH,:))./sqrt(sum(isfinite(alltsvals1(subH,:)))),'rs-')
% hold on
errorbar(1:trials, nanmean(alltsvals2(subH,:)), nanstd(alltsvals2(subH,:))./sqrt(sum(isfinite(alltsvals1(subH,:)))),'bs-')
xlim([1 trials])
ylim(yrange)


%% plot average, per epoch , VH and H
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0 0 5 5],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group VH')
hold on
%%e=11,22,11,22 , e=12,21,12,21
%e1=errorbar(1:tsets, mean(allvals11(subVH,:)), std(allvals11(subVH,:))./sqrt(length(subVH)),'rs-')
e2=errorbar(1:tsets, mean(allvals1(subVH,:)), std(allvals1(subVH,:))./sqrt(length(subVH)),'bs-')
hold on
e5=errorbar(1:tsets, mean(allvals2(subVH,:)), std(allvals2(subVH,:))./sqrt(length(subVH)),'rs-')
%e6=errorbar(1:tsets, mean(allvals22(subVH,:)), std(allvals22(subVH,:))./sqrt(length(subVH)),'bs-')



l1=legend([e2,e5],'Subj_{1}','Subj_{2}')
set(l1,'Location','north')
legend boxoff
ylim([0 4])
%%Plot average, per-epoch (12,21)
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0 0 5 5],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
ylabel(plotlabel)
title('group H')
hold on

e4=errorbar(1:tsets, mean(allvals1(subH,:),1), std(allvals1(subH,:),1)./sqrt(length(subH)),'bs-')
hold on
e7=errorbar(1:tsets, mean(allvals2(subH,:),1), std(allvals2(subH,:),1)./sqrt(length(subH)),'rs-')
%l1=legend([e5,e6,e7,e8],'VH:Subj 1','VH:Subj 2','H: Subj 1','H: Subj 2')
%l2=legend([e4,e7],'md_{21}','md_{12}')
%set(l2,'Location','north')
ylim([0 4])
legend boxoff


%% plot average, per epoch 'interaction effort' , VH and H
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0.02 0.02 10 10],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% xlabel('epochs','fontsize',14)
% ylabel('Interaction force[N]','fontsize',14)
title('Group VH','fontsize',14,'FontName', 'Times New Roman')
hold on
%%e=11,22,11,22 , e=12,21,12,21
%e1=errorbar(1:tsets, mean(allvals11(subVH,:)), std(allvals11(subVH,:))./sqrt(length(subVH)),'rs-')
e(2)=errorbar(1:tsets, mean(allvalsint(subVH,:)), std(allvalsint(subVH,:))./sqrt(length(subVH)),'ks-')

%e6=errorbar(1:tsets, mean(allvals22(subVH,:)), std(allvals22(subVH,:))./sqrt(length(subVH)),'bs-')
set(e(2),'Linewidth',1)
 xlabel('epochs','fontsize',14,'FontName', 'Times New Roman')
 ylabel('Interaction force [N]','fontsize',14,'FontName', 'Times New Roman')

% l1=legend([e2,e5],'Subj_{1}','Subj_{2}')
% set(l1,'Location','north')
% legend boxoff
ylim([0 10])
%%Plot average, per-epoch (12,21)
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0.02 0.02 10 10],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
 xlabel('epochs','fontsize',14,'FontName', 'Times New Roman')
 ylabel('Interaction force [N]','fontsize',14,'FontName', 'Times New Roman')
% ylabel('Interaction force[N]','fontsize',14)
title('Group H','fontsize',14,'FontName', 'Times New Roman')
hold on

e(4)=errorbar(1:tsets, mean(allvalsint(subH,:),1), std(allvalsint(subH,:),1)./sqrt(length(subH)),'ks-')
set(e(4),'Linewidth',1)
% hold on
% e7=errorbar(1:tsets, mean(allvals2(subH,:),1), std(allvals2(subH,:),1)./sqrt(length(subH)),'rs-')
%l1=legend([e5,e6,e7,e8],'VH:Subj 1','VH:Subj 2','H: Subj 1','H: Subj 2')
%l2=legend([e4,e7],'md_{21}','md_{12}')
%set(l2,'Location','north')
ylim([0 10])
% legend boxoff

%% plot pre-post 'interaction effort' , VH and H
figname = 'eff_prepost_all_star'
figure
set(gcf,'pos',[100 100 250 300])
%title('Subject 1')
hold on
errorbar(2+[-0.15 0.15], mean(allvalsint(subVH,[2 11])), ...
                         std(allvalsint(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(allvalsint(subH,[2 11])), ...
                         std(allvalsint(subH,[2 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvals1(subH,2)) mean(allvalsint(subH,11)); ...
     mean(allvals1(subVH,2)) mean(allvalsint(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group','fontsize',14)
ylabel('Interaction force','fontsize',14)
sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15]},[0.005, 0.005])
ylim([0 10])
%ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
saveas(gcf,[figdir,figname],'fig');
eval(['print -depsc ', epsdir,figname])

%% shaded error bar try
figure
set(gcf,'pos',[100 100 250 300])
patch([1.5 11.5 11.5 1.5],[0 0 10 10],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
myeb( allvals1(subVH,:))
ylim([0 10])
ylabel(plotlabel)
title('group VH')
figure
patch([1.5 11.5 11.5 1.5],[0 0 5 5],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('epochs')
set(gcf,'pos',[100 100 250 350])
myeb(allvals1(subH,:))
ylim([0 4])
title('group H')
ylabel(plotlabel)

%% Patch try
figure
set(gcf,'pos',[100 100 250 300])
patch([12.5 132.5 132.5 12.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel(plotlabel)
title('H group')
hold on


%s3=errorbar(1:trials, nanmean(alltsvals11(subH,:),1), nanstd(alltsvals11(subH,:),1)./sqrt(sum(isfinite(alltsvals11(subH,:)))),'cs-')


%hold on
%s4=errorbar(1:trials, nanmean(eff1_ts(subH,:),1), nanstd(eff1_ts(subH,:),1)./sqrt(sum(isfinite(alltsvals1(subH,:)))),'rs:')
hold on
patch([1:trials trials:-1:1],[nanmean(eff1_ts(subH,:),1)-nanstd(eff1_ts(subH,:),1) ...
                        fliplr(nanmean(eff1_ts(subH,:),1)+nanstd(eff1_ts(subH,:),1))],0.5+[0.5 0 0.5],'edgecol',0.5+[0.5 0 0.5]);
%s4=line(1:trials, nanmean(md21_ts(subH,:),1), 'col','r','marker','s','lines','-')


%s7=errorbar(1:trials, nanmean(eff2_ts(subH,:),1), nanstd(eff2_ts(subH,:),1)./sqrt(sum(isfinite(alltsvals2(subH,:)))),'bs:')
% hold on
% patch([1:trials trials:-1:1],[nanmean(eff2_ts(subH,:),1)-nanstd(eff2_ts(subH,:),1) ...
%                         fliplr(nanmean(eff2_ts(subH,:),1)+nanstd(eff2_ts(subH,:),1))],0.3+[0.3 0 0.3],'edgecol',0.3+[0.3 0 0.3]);
%t2=legend([s3,s4,s7,s8],'S1-VP1','S2-VP1','S1-VP2','S2-VP2')
%t2=legend([s4,s7],'md_{21}','md_{12}')
%set(t2,'Location','north')
legend boxoff
xlim([1 trials])
ylim(yrange)

% %%
% figure
% for i=1:5
% % a1 = plot(intforce(i,[2 11]),R(i,[2 11]),'b-','LineWidth',2),hold on
% % plot(intforce(i,[11]),R(i,[11]),'b*','LineWidth',2),hold on
% a1 = scatter(R(i,[11]),intforce(i,[11]),'b'),hold on
% end
% for i=7:11
% % a2 = plot(intforce(i,[2 11]),R(i,[2 11]),'r-','LineWidth',2),hold on
% % plot(intforce(i,[11]),R(i,[11]),'r*','LineWidth',2),hold on
% a2 = scatter(R(i,[11]),intforce(i,[11]),'r'),hold on
% end
% % t1=legend([a1,a2],'H','VH')
% % legend boxoff
% xlabel('corr coeff','fontsize',14)
% ylabel('int force','fontsize',14)
% title('new')
% %xlim([0 0.05])
% %ylim([0 0.05])
% axis square


%% fitting try
meanVH_epoch = mean(allvalsint(subVH,:));
meanH_epoch = mean(allvalsint(subH,:));
meanVH_trial = nanmean(alltsvalsint(subVH,:));
meanH_trial = nanmean(alltsvalsint(subH,:));

[VH_hat, VH_c] = expfit(meanVH_epoch);
[H_hat, H_c] = expfit(meanH_epoch);

% fit_VH = fit((2:11)', (meanVH_epoch(:,2:11))','exp1', 'StartPoint',[2.5,1.5]);
% fit_H = fit((2:11)', (meanH_epoch(:,2:11))','exp1', 'StartPoint',[2.5,1.5]);
fit_VH = fit((2:11)', (meanVH_epoch(:,2:11))','exp1');
fit_H = fit((2:11)', (meanH_epoch(:,2:11))','exp1');
[fitt_VH, gof_VH] = fit((13:132)', (meanVH_trial(:,13:132))','exp2'); % considering goodness of fit
[fitt_H, gof_H] = fit((13:132)', (meanH_trial(:,13:132))','exp2');
%% epoch wise
figure
a(1) = plot(meanVH_epoch,'ro')
hold on
a(2) = plot(meanH_epoch,'bo')
plot(fit_VH,'r-')
plot(fit_H,'b-')
legend(a,{'VH','H'})
xlabel('epochs')
ylabel('mean effort[N]')
% title('')
legend boxoff 
axis square

%% trial wise
figure
set(gcf,'pos',[100 100 250 350])
title('group VH')
patch([12.5 132.5 132.5 12.5],[0 0 10 10],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
hold on
b1= plot(meanVH_trial,'ro')
hold on
b2= plot(fitt_VH,'r-')
legend(b2,{'fitted curve for VH'})
xlabel('trials')
ylabel('mean effort[N]')
ylim([0 10])
legend boxoff 
%axis square
str3 = sprintf('r^2 = %f',gof_VH.rsquare);
text(20,1.5,str3 ,'FontSize',10)


figure
set(gcf,'pos',[100 100 250 350])
title('group H')
patch([12.5 132.5 132.5 12.5],[0 0 10 10],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
hold on
b3 = plot(meanH_trial,'bo')
hold on
b4=plot(fitt_H,'b-')
legend([b4],{ 'fitted curve for H'})
xlabel('trials')
ylabel('mean effort[N]')
legend boxoff 
ylim([0 10])
%axis square
str2 = sprintf('r^2 = %f',gof_H.rsquare);
text(20,1.5,str2 ,'FontSize',10)


%% subject-wise fitting




coeffs_exp2_b = zeros(11,1);
coeffs_exp2_d = zeros(11,1);
coeffs_exp2_a = zeros(11,1);
coeffs_exp2_c = zeros(11,1);
coeffs_exp1 = zeros(11,1);
coeffs_expfit = zeros(11,1);

for subj=1:subjs
    
    alltsvals = alltsvalsint(subj,13:132);
    alltsvals_corrected = alltsvals;
    alltsvals_corrected(isnan(alltsvals))=[];
    
    trials = 1:length(alltsvals);
    trials_corrected = trials;
    trials_corrected(isnan(alltsvals))=[];
    
    

    [fit_per_sub_exp2, ggf2] = fit(trials_corrected', alltsvals_corrected','exp2');
    [fit_per_sub_exp1, ggf1] = fit(trials_corrected', alltsvals_corrected','exp1');
%   coeffs_exp2(subj,:) = coeffvalues(fit_per_sub);
    coeffs_exp2_b(subj,:) = fit_per_sub_exp2.b;
    coeffs_exp2_d(subj,:) = fit_per_sub_exp2.d;
    coeffs_exp2_a(subj,:) = fit_per_sub_exp2.a;
    coeffs_exp2_c(subj,:) = fit_per_sub_exp2.c;
    coeffs_exp1(subj,:) = fit_per_sub_exp1.b;
    %coeffs_exp1(subj,:) = coeffvalues(fit_per_sub);
    coeffs_expfit(subj,:) = expfit(alltsvals_corrected);
    
    figure
    set(gcf,'pos',[100 100 300 400])
%     patch([12.5 132.5 132.5 12.5],[0 0 10 10],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
%     hold on;
    xlabel('trials')
    ylabel(plotlabel)
    title(sprintf('Dyad%d',subj))
%   plot(alltsvals_corrected,'bo')
    plot(alltsvals_corrected,'bo')
    hold on
    a=plot(fit_per_sub_exp2,'r-')
    b=plot(fit_per_sub_exp1,'b-')
    xlabel('trials','fontsize',14)
    ylabel('effort[N]','fontsize',14)
    legend([a,b],{ 'fitted with exp2', 'fitted with exp1'})
    %xlim([1 trials])
    ylim([0 10])
    str = sprintf('r^2 = %f',ggf1.rsquare);
    text(20,1.5,str ,'Color','b','FontSize',10)
    str = sprintf('r^2 = %f',ggf2.rsquare);
    text(20,2.5,str ,'Color','r','FontSize',10)


end


%% coeffs
one= ones(11,1);
two= 2*ones(11,1);
three= 3*ones(11,1);
four = 4*ones(11,1);

figure
set(gcf,'pos',[100 100 250 300]);
plot(one, coeffs_exp2_a,'ro');
hold on
plot(two, coeffs_exp2_b, 'bo');
plot(three, coeffs_exp2_c, 'go');
plot(four, coeffs_exp2_d, 'yo');
set(gca,'xtick',[1 2 3 4],'xticklabel',{'a', 'b', 'c', 'd'})
%set(gca,'xtick',[1 2],'xticklabel',{'a', 'b'})
xlim([0 5])



%% coeffs values
figure
% for i=1:5
% c(1) = plot(coeffs_expfit(i,:),'ro')
% hold on
% end 
% 
% for i= 7:11
% c(2) = plot(coeffs_expfit(i,:),'bo')
% end
c(1) = plot(coeffs_expfit(subH,:),'ro')
hold on
c(2) = plot(coeffs_expfit(subVH,:),'bo')
legend(c,{ 'fitting coeff for H','fitting coeff for VH'})
xlim([0 6])
ylim([0 2])
axis square

%% coeffs values a, b,c d
figure
% for i=1:5
% c(1) = plot(coeffs_expfit(i,:),'ro')
% hold on
% end 
% 
% for i= 7:11
% c(2) = plot(coeffs_expfit(i,:),'bo')
% end
c(1) = plot(coeffs_exp2_a(subH,:))
hold on
c(2) = plot(coeffs_exp2_a(subVH,:))
% legend(c,{ 'fitting coeff for H','fitting coeff for VH'})
% xlim([0 6])
% ylim([0 2])


% set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
% hold on;
% 
% plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(ydata, 1), 2, 1), 'k-')

% ylim([0 max(ydata(:)+1)])
% ylim([-5 5])



%% Stats analysis 
group = zeros(subjs,1);
tt = [coeffs_exp2_b,coeffs_exp2_d];
tt_exp = sort(tt,2);
StatMatrix2 = [group, coeffs_expfit,  coeffs_exp1, tt_exp];
StatMatrix2(6,:)=[]; %take out subject 6
StatMatrix2(subH,1)=1;
StatMatrix2(6:end,1)=2;
save('fitting_coeffs.txt','StatMatrix2', '-ascii')

%% error bar plot for fitting coefficients

figname = 'eff_fitting'
figure
set(gcf,'pos',[100 100 250 300])
%title('fitting coefficents')
hold on
errorbar(2, mean(coeffs_expfit(subVH,:)), ...
                         std(coeffs_expfit(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1, mean(coeffs_expfit(subH,:)), ...
                         std(coeffs_expfit(subH,:))./sqrt(length(subH)),'k.')
 
h=bar([mean(coeffs_expfit(subH,:)) ; ...
mean(coeffs_expfit(subVH,:))])
sigstar({[1,2]},[0.05])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group','fontsize',14)
ylim([0 2])
ylabel('Mean of fitting coefficient','fontsize',14)
% set(h(1),'facecol','w')
set(h,'facecol',[0.8 0.8 0.8])
saveas(gcf,[figdir,figname],'fig');
eval(['print -depsc ', epsdir,figname])





%% plot pre-post 'interaction effort' , VH and H
% name = sprintf('intforce_stat');
fullFileName = fullfile(folder, name);

figure
set(gcf,'pos',[100 100 250 300])
%title('Subject 1')
hold on
errorbar(2+[-0.225 0 0.225], mean(allvalsint(subVH,[2 6 11])), ...
                         std(allvalsint(subVH,[2 6 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(allvalsint(subH,[2 6 11])), ...
                         std(allvalsint(subH,[2 6 11]))./sqrt(length(subH)),'k.')

h=bar([mean(allvalsint(subH,2)) mean(allvalsint(subH,6)) mean(allvalsint(subH,11)); ...
     mean(allvalsint(subVH,2))  mean(allvalsint(subVH,6)) mean(allvalsint(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
 xlabel('Group','fontsize',14,'FontName', 'Times New Roman')
 ylabel('Interaction force [N]','fontsize',14,'FontName', 'Times New Roman')
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
ylim([0 11.5])
%ylim(yrange)
 legend(h,{'Early training','Middle training','Late training'}, 'Location', 'northwest')
 legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])
print(gcf,fullFileName,'-depsc', '-r300')




%% shaded error bar with trasperancy using fcn (see shadedErrorBar function in the directory)
name = sprintf('intforce_H');
fullFileName = fullfile(folder, name);
colB = [.2 .2 .2];
figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.02 0.02 12 12],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Interaction force [N]', 'fontsize', 14,'FontName', 'Times New Roman')
title('H', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
%  e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subH,:),1),nanstd(alltsvalsint(subH,:)./sqrt(length(subH)),1),{'-k', 'LineWidth',1},1);
e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subH,:),1),nanstd(alltsvalsint(subH,:)./sqrt(length(subH)),1),{'LineStyle','none','Color',colB,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
% catch trials..
% errorbar(center_trials,nanmean(md21_statmat_nofor(subVH,2:end)),...
%     nanstd(md21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)*0),'Color',[.9 .4 .4],'LineStyle',':', 'LineWidth',2)
% errorbar(center_trials,nanmean(md12_statmat_nofor(subVH,2:end)),...
%     nanstd(md12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)*0),'Color',[.4 .4 .9],'LineStyle',':', 'LineWidth',2)

xlim([1 156])
ylim([0 11.5])
% legend([e1.mainLine,e2.mainLine],'Subject 1','Subject 2','FontName', 'Times New Roman')
% legend boxoff
print(gcf,fullFileName,'-depsc', '-r300')

name = sprintf('intforce_VH');
fullFileName = fullfile(folder, name);

figure
set(gcf,'pos',[100 100 230 300])
patch([12 132 132 12],[0.02 0.02 12 12],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
% ylabel('Interaction force [N]', 'fontsize', 14,'FontName', 'Times New Roman')
title('VH', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
% shadedErrorBar([1:156],nanmean(alltsvalsint(subVH,:),1),nanstd(alltsvalsint(subVH,:)./sqrt(length(subVH)),1),{'-k', 'LineWidth',1},1);
e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subVH,:),1),nanstd(alltsvalsint(subVH,:)./sqrt(length(subVH)),1),{'LineStyle','none','Color',colB,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
xlim([1 156])
ylim([0 11.5])
pos = get(gca, 'Position'); % space adustment without ylabel
pos(1) = 0.15;
pos(3) = 0.78;
set(gca, 'Position', pos)
print(gcf,fullFileName,'-depsc', '-r300')

name = sprintf('intforce_PV');
fullFileName = fullfile(folder, name);

figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.02 0.02 12 12],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Interaction force [N]', 'fontsize', 14,'FontName', 'Times New Roman')
title('PV', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
%  e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subPV,:),1),nanstd(alltsvalsint(subPV,:)./sqrt(length(subPV)),1),{'-k', 'LineWidth',1},1);
 e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subPV,:),1),nanstd(alltsvalsint(subPV,:)./sqrt(length(subPV)),1),{'LineStyle','none','Color',colB,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);

%          
%                         
% errorbar(center_trials,nanmean(md21_statmat_nofor(subVH,2:end)),...
%     nanstd(md21_statmat_nofor(subVH,2:end))./sqrt(length(subVH)*0),'Color',[.9 .4 .4],'LineStyle',':', 'LineWidth',2)
% errorbar(center_trials,nanmean(md12_statmat_nofor(subVH,2:end)),...
%     nanstd(md12_statmat_nofor(subVH,2:end))./sqrt(length(subVH)*0),'Color',[.4 .4 .9],'LineStyle',':', 'LineWidth',2)

xlim([1 156])
ylim([0 11.5])
% legend([e1.mainLine,e2.mainLine],'Subject 1','Subject 2','FontName', 'Times New Roman')
% legend boxoff


%% Three groups in 1
c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
training= 12:132;
figure
set(gcf,'pos',[100 100 250 300])
patch([12 132 132 12],[0.02 0.02 12 12],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('Trials', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Interaction force [N]', 'fontsize', 14,'FontName', 'Times New Roman')
% title('PV', 'fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
%e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subPV,:),1),nanstd(alltsvalsint(subPV,:)./sqrt(length(subPV)),1),{'-k', 'LineWidth',1},1);
e1=shadedErrorBar([1:156],nanmean(alltsvalsint(subH,:),1),nanstd(alltsvalsint(subH,:)./sqrt(length(subH)),1),{'LineStyle','none','Color',c1,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
hold on
e2=shadedErrorBar([1:156],nanmean(alltsvalsint(subVH,:),1),nanstd(alltsvalsint(subVH,:)./sqrt(length(subVH)),1),{'LineStyle','none','Color',c2,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
hold on
e3=shadedErrorBar([1:156],nanmean(alltsvalsint(subPV,:),1),nanstd(alltsvalsint(subPV,:)./sqrt(length(subPV)),1),{'LineStyle','none','Color',c3,'Marker','.', 'MarkerSize',15 'LineWidth',1},1);
% h=legend([e1.mainLine,e2.mainLine,e3.mainLine],'H','VH','PV','FontName', 'Times New Roman')

rect = [0.63, 0.75, .1, .1]; % to control position of legend 
set(h, 'Position', rect)
legend boxoff
xlim([1 156])
ylim([0 11.5])
name = sprintf('intforce_all');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%% Stats analysis (epoch wise)

group = zeros(subjs,1);
StatMatrix = [group, allvalsint(:,[2 6 11])];
StatMatrix(6,:)=[]; % take out subject 6
StatMatrix(subH,1)=1;
StatMatrix(6:end,1)=2;
save('intforce_epoch1.txt','StatMatrix', '-ascii')

%% trial ststistics...
alltsvals1_1 = nanmean(alltsvalsint(:, 1:12),2)
alltsvals1_2 = nanmean(alltsvalsint(:, 13:24),2)
alltsvals1_6 = nanmean(alltsvalsint(:, 60:72),2)
alltsvals1_11 = nanmean(alltsvalsint(:, 120:132),2)
alltsvals1_12 = nanmean(alltsvalsint(:, 133:156),2)
intstat = [alltsvals1_2 alltsvals1_6 alltsvals1_11]

    

%% figure

% name = sprintf('score1_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
% title('Interaction force','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(intstat(subVH,:)), ...
                         std(intstat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(intstat(subH,:)), ...
                         std(intstat(subH,:))./sqrt(length(subH)),'k.')

h=bar([ mean(intstat(subH,1)) mean(intstat(subH,2)) mean(intstat(subH,3)) ; ...
      mean(intstat(subVH,1)) mean(intstat(subVH,2)) mean(intstat(subVH,3)) ])
 
 
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1+0.225,2+0.225]},[0.0005, 0.0005, 0.05])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Interaction force [N]', 'fontsize', 14,'FontName', 'Times New Roman')
 ylim([0 11.5])
 legend(h,{'Early training','Middle training','Late training'}, 'Location', 'northwest')
 legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])

set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
set(h(3),'facecol',[0.4 0.4 0.4])

name = sprintf('intforce_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%% three
% name = sprintf('intforce_all_stat');
% fullFileName = fullfile(folder, name);
figure
%title('md_{21}')
% title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(2+[-0.225 0 0.225], mean(intstat(subVH,:)), ...
                         std(intstat(subVH,:))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.225 0 0.225], mean(intstat(subH,:)), ...
                         std(intstat(subH,:))./sqrt(length(subH)),'k.')
 hold on
errorbar(3+[-0.225 0 0.225], mean(intstat(subPV,:)), ...
                         std(intstat(subPV,:))./sqrt(length(subPV)),'k.')

h=bar([ mean(intstat(subH,1)) mean(intstat(subH,2)) mean(intstat(subH,3)) ; ...
      mean(intstat(subVH,1)) mean(intstat(subVH,2)) mean(intstat(subVH,3));
      mean(intstat(subPV,1)) mean(intstat(subPV,2)) mean(intstat(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
 sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1+0.225,2+0.225],[3+0.225, 2+0.225], [3+0.225, 1+0.225],[3+0.225, 3-0.225] },[0.0005, 0.0005, 0.05, 0.05, 0.05,0.0005])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH' , 'PV'},'FontName', 'Times New Roman')
xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Interaction force [N]', 'fontsize', 14,'FontName', 'Times New Roman')
 ylim([0 11.5])
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
print(gcf,fullFileName,'-depsc', '-r300')



%% figure 4

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
errorbar(1+[-0.225 0 0.225], [mean(intstat(subH,1)) mean(intstat(subVH,1)) mean(intstat(subPV,1))], ...
                         [std(intstat(subH,1))./sqrt(length(subH)) std(intstat(subVH,1))./sqrt(length(subVH)) std(intstat(subPV,1))./sqrt(length(subPV))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [mean(intstat(subH,2)) mean(intstat(subVH,2)) mean(intstat(subPV,2))], ...
                         [std(intstat(subH,2))./sqrt(length(subH)) std(intstat(subVH,2))./sqrt(length(subVH)) std(intstat(subPV,2))./sqrt(length(subPV))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [mean(intstat(subH,3)) mean(intstat(subVH,3)) mean(intstat(subPV,3))], ...
                         [std(intstat(subH,3))./sqrt(length(subH)) std(intstat(subVH,3))./sqrt(length(subVH)) std(intstat(subPV,3))./sqrt(length(subPV))],'k.')
 hold on
h=bar([mean(intstat(subH,1)) mean(intstat(subVH,1)) mean(intstat(subPV,1)) ; ...
      mean(intstat(subH,2)) mean(intstat(subVH,2)) mean(intstat(subPV,2));
      mean(intstat(subH,3)) mean(intstat(subVH,3)) mean(intstat(subPV,3))])
 
 
%   sigstar({[1-0.225,1+0.225],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.05, 0.0005, 0.005])
   sigstar({[2-0.225,2+0.225], [3,3-0.225], [3-0.225,3+0.225]},[0.05, 0.05, 0.05])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'},'FontName', 'Times New Roman')
xlabel('Epoch', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Interaction force [N]', 'fontsize', 14,'FontName', 'Times New Roman')
 ylim([0 11.5])
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
% set(h(1),'facecol',c1)
% set(h(2),'facecol',c2)
% set(h(3),'facecol',c3)

name = sprintf('intforce_all_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')




%% Stat values
S_1 = [intstat]; 
S_1(6,:) = [];
group_S = [eff1_statmat(:,1)];
group_S(6,:)=[];
subjects_S_1 = (1:length(S_1(:,1)))';
pre_S_1 = [ S_1(:,1) group_S ones(length(group_S),1) subjects_S_1];
mid_S_1 = [ S_1(:,2) group_S 2*ones(length(group_S),1) subjects_S_1];
pos_S_1 = [ S_1(:,3) group_S 3*ones(length(group_S),1) subjects_S_1];
A_S_1 = [pre_S_1; mid_S_1; pos_S_1];
name = sprintf('intf_three_group_fourcolumn.txt');
fullFileName = fullfile(stat_folder, name);
save(fullFileName,'A_S_1', '-ascii')
BWAOV2(A_S_1);
