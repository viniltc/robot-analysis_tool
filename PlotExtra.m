

%% Extra plots..

clear all
close all

% Set location of data
datadir = '../../Data/';
% Details on where are results
resdir = [datadir, 'results/'];
epsdir = [datadir, 'eps/'];
folder = 'C:\Users\Asus\Desktop\journal_photoshop\illustrator\exp';
folder = 'C:\Users\Asus\Dropbox\pHHI sfn paper\Data\eps\new\panel_sep\';

subH = 1:5;
subVH = 7:11;  % 5 removed
subPV = 12:16;  

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

load([resdir, 'pow1_statmat.txt']);
load([resdir, 'pow1_statmat_catch.txt']);
load([resdir, 'pow2_statmat.txt']);
load([resdir, 'pow2_statmat_catch.txt']);
load([resdir, 'pow_fac_statmat.txt']);
load([resdir, 'pow_fac_statmat_catch.txt']);

load([resdir, 'lead1_statmat.txt']);
load([resdir, 'lead1_statmat_catch.txt']);
load([resdir, 'lead2_statmat.txt']);
load([resdir, 'lead2_statmat_catch.txt']);

load ([resdir,'forcenorm.mat']);
forcenorm_ts = sqrt(ind_ts);
load ([resdir,'eff1.mat']);
eff1_ts = sqrt(ind_ts);
load ([resdir,'eff2.mat']);
eff2_ts = sqrt(ind_ts);
load([resdir, 'forces.mat']);

load ([resdir,'pow1.mat']);
pow1_ts = (ind_ts);
load ([resdir,'pow2.mat']);
pow2_ts = (ind_ts);
load ([resdir,'pow_fac.mat']);
pow_fac_ts = (ind_ts);


load ([resdir,'lead1.mat']);
lead1_ts = (ind_ts);
load ([resdir,'lead2.mat']);
lead2_ts = (ind_ts);

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

%catch trial
md11_ts_c=md11_ts(find(forces==0));
md22_ts_c=md22_ts(find(forces==0));
md12_ts_c=md12_ts(find(forces==0));
md21_ts_c=md21_ts(find(forces==0));

md11 = md11_statmat(:,2:end);
md12 = md12_statmat(:,2:end);
md21 = md21_statmat(:,2:end);
md22 = md22_statmat(:,2:end);
pow1 = pow1_statmat(:,2:end);
pow2 = pow2_statmat(:,2:end);
pow_fac = pow_fac_statmat(:,2:end);
lead1 = lead1_statmat(:,2:end);
lead2 = lead2_statmat(:,2:end);

lead1(isnan(lead1))=0;
lead2(isnan(lead2))=0;

forcenorm = forcenorm_statmat(:,2:end);
md_avg = md12+md21; % minimumdistance average

md = [md12,md21];
subjs=size(pow1,1);
tsets=size(pow1,2);
trials = size(pow1_ts,2);


%% learning traces


col1 =[.1 .6 .3]
col2 =[.6 .1 .4]
col1 =[.1 .6 .3];
col2 =[.6 .1 .4];

c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
figure
set(gcf,'pos',[100 100 300 300])
for i=1:5
a1 = plot(md21(i,[2 11]),md12(i,[2 11]),'-','LineWidth',1, 'Color',c1),hold on
plot(md21(i,[11]),md12(i,[11]),'.','MarkerSize',20,'LineWidth',1, 'Color', c1),hold on
% plot(md21(i,[11]),md12(i,[11]),'.','MarkerSize',20,'LineWidth',1, 'Color', w),hold on
labelpoints(md21(i,11), md12(i,11), i, 'Color', c1, 'FontSize', 8)
end
for i=7:11
a2 = plot(md21(i,[2 11]),md12(i,[2 11]),'-','LineWidth',1, 'Color', c2),hold on
plot(md21(i,[11]),md12(i,[11]),'.','MarkerSize',20,'LineWidth',1, 'Color',c2),hold on
labelpoints(md21(i,11), md12(i,11), i, 'Color', c2, 'FontSize', 8)
end
for i=12:16
a3 = plot(md21(i,[2 11]),md12(i,[2 11]),'-','LineWidth',1, 'Color', c3),hold on
plot(md21(i,[11]),md12(i,[11]),'.','MarkerSize',20,'LineWidth',1, 'Color',c3),hold on
labelpoints(md21(i,11), md12(i,11), i, 'Color', c3, 'FontSize', 8)
end

t1=legend([a1,a2, a3],'H','VH','PV','Late training', 'Location', 'northwest')
legend boxoff
xlabel('MD_{21}[m]','fontsize',14,'FontName', 'Times New Roman')
ylabel('MD_{12}[m]','fontsize',14,'FontName', 'Times New Roman')
% title('Learning traces')
xlim([0 0.06])
ylim([0 0.06])
axis square
box off
name = sprintf('learning');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')
%% asd score vs md (md_21 vs ASD)

figure
for i=1:5
a1 = plot(md21(i,[11]),asd_subj2(i,[6]),'b*'), hold on
end
for i=7:11
a2 = plot(md21(i,[11]),asd_subj2(i,[6]),'r*'), hold on
end
%t1=legend([a1,a2],'H','VH')
%legend boxoff
xlabel('md_{21}')
ylabel('ASD score')
xlim([0 0.05])
ylim([5 30])
title('ASD score vs minimum distance:subject 2')
axis square

%% asd score vs md (md_12 vs ASD)


figure
for i=1:5
b1 = plot(md12(i,[11]),asd_subj1(i,[6]),'b*'), hold on
end
for i=7:11
b2 = plot(md12(i,[11]),asd_subj1(i,[6]),'r*'), hold on
end
%t1=legend([b1,b2],'H','VH')
%legend boxoff
xlabel('md_{12}')
ylabel('ASD score')
ylim([5 30])
xlim([0 0.05])
title('ASD score vs minimum distance:subject 1')
axis square

%% intforce vs avg minimum distance (collaboration vs compliance with task)


figure
for i=1:5
%b1 = plot(forcenorm(i,[11]),md_avg(i,[11]),'b*'), hold on
a1 = scatter(forcenorm(i,1),md_avg(i,1),'MarkerEdgeColor',[.7 0 .7],...
              'MarkerFaceColor',[.9 0 .7],...
              'LineWidth',1), hold on
a2 = scatter(forcenorm(i,2),md_avg(i,2),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
a3 = scatter(forcenorm(i,11),md_avg(i,11),'MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
a4 = scatter(forcenorm(i,12),md_avg(i,12),'MarkerEdgeColor',[.5 .5 0],...
              'MarkerFaceColor',[.9 .7 0],...
              'LineWidth',1), hold on
a5 = scatter(-1,-1,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1), hold on
end

for i=7:11
%b2 = plot(forcenorm(i,[11]),md_avg(i,[11]),'r*'), hold on
b1 = scatter(forcenorm(i,1),md_avg(i,1),'d','MarkerEdgeColor',[.5 .5 0],...
              'MarkerFaceColor',[.9 .7 0],...
              'LineWidth',1), hold on
b2 = scatter(forcenorm(i,2),md_avg(i,2),'d','MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
b3 = scatter(forcenorm(i,11),md_avg(i,11),'d','MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
b4 = scatter(forcenorm(i,12),md_avg(i,12),'d','MarkerEdgeColor',[.7 0 .7],...
              'MarkerFaceColor',[.9 0 .7],...
              'LineWidth',1), hold on
b5 = scatter(-1,-1,'d','MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1), hold on

end
%t1=legend([a5,b5],'H','VH')
t2=legend([a1,a2,a3,a4,a5,b5],'Baseline','Training inital','Training final','Washout','H group','VH group')
legend boxoff
xlabel('interaction force [N]')
ylabel('average minimum distance [m]')
ylim([0 0.1])
xlim([0 4])
title('interaction force vs avg. minimum distance ')
axis square

%% Calculated score (min distance and interaction_force)

% trials = size(md11_ts,2);
% strials = 24;
% % forcenorm_ts_c=size(forcenorm_ts(find(forces==0));
% 
% K = 2*log(100*2-1)/(2-0.25);
% score1_ts = md11_ts+0.5*forcenorm_ts;
% score2_ts = md22_ts+0.5*forcenorm_ts;
% score1_ep = md11+0.5*forcenorm;
% score2_ep = md22+0.5*forcenorm;
% 
% cal_score1_ts = (100./(1+exp(K*(score1_ts-(2+0.25)/2))));
% cal_score2_ts = (100./(1+exp(K*(score2_ts-(2+0.25)/2))));
% cal_score1_ep = (100./(1+exp(K*(score1_ep-(2+0.25)/2))));
% cal_score2_ep = (100./(1+exp(K*(score2_ep-(2+0.25)/2))));
% 
% %% Calculated score per trial
% 
% 
% for subj=1:subjs
%        figure
%        plot(1:trials,cal_score1_ts(subj,:),'-rs'), hold on
%        plot(1:trials,cal_score2_ts(subj,:),'-bs');
%        xlabel('trials')
%        ylabel('calculated score')
%        title(sprintf('Dyad %d',subj))
%        xlim([1 trials])
%        ylim([1 100])
%        axis square
% end
% 
% %% Calculated score per epoch
% for subj=1:subjs
%        
%     figure
%    
%        plot(1:tsets,cal_score1_ep(subj,:),'-rs'), hold on
%        plot(1:tsets,cal_score2_ep(subj,:),'-bs');
%        xlabel('epoch')
%        ylabel('calculated score')
%        title(sprintf('Dyad %d',subj))
%        xlim([1 tsets])
%        ylim([1 100])
%        axis square
% end
% 
% 
% %% From strored score (*.pro file)
% foldername ={'../../Data/K_S/', '../../Data/L_J/', '../../Data/P_O/','../../Data/P_At/','../../Data/A_J/','../../Data/N_F/','../../Data/Z_A/','../../Data/P_C/','../../Data/M_F/','../../Data/P_A/','../../Data/S_N/'};
% extension='*.pro';
% score_epoch1 = zeros(13,3);
% score_epoch2 = zeros(13,3);
% 
% for k = 1:length(foldername)
% 
%     name=strcat(foldername{k}, extension);
%     fileset=dir(name);
%     a = zeros(12,3);                               
%   
%    for i = 1:length(fileset)
%         a=load(strcat(foldername{k}, fileset(i).name)); 
%      
% %         figure
% %         
% %         sub1=plot(a(:,1),a(:,2),'-bs'), hold on
% %         sub2=plot(a(:,1),a(:,3),'-rs')
% %         ylabel('score')
% %         xlabel('trial')
% %         legend([sub1, sub2],'sub1','sub2')
% %         legend boxoff
% %         title(sprintf(['Dyad' strrep(foldername{k}, '_', ' ') ': %d'],i));
% %         axis square
%         score_epoch1(i,k) = mean(a(:,2));
%         score_epoch2(i,k) = mean(a(:,3));
%     end
%     
%     figure
%     ep1 = plot (score_epoch1(:,k),'--bs'), hold on
%     ep2 = plot (score_epoch2(:,k),'--rs')
%     ylabel('score per epoch')
%     xlabel('epoch')
%     legend([ep1, ep2],'sub1','sub2')
%     legend boxoff
%     title(sprintf(foldername{k}));
%     axis square
%     
% end


%% intforce vs avg minimum distance (collaboration vs compliance with task)


figure
for i=1:5
%b1 = plot(forcenorm(i,[11]),md_avg(i,[11]),'b*'), hold on
a1 = scatter(forcenorm(i,1),md12(i,1),'MarkerEdgeColor',[.7 0 .7],...
              'MarkerFaceColor',[.9 0 .7],...
              'LineWidth',1), hold on
a2 = scatter(forcenorm(i,2),md12(i,2),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
a3 = scatter(forcenorm(i,11),md12(i,11),'MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
a4 = scatter(forcenorm(i,12),md12(i,12),'MarkerEdgeColor',[.5 .5 0],...
              'MarkerFaceColor',[.9 .7 0],...
              'LineWidth',1), hold on
a5 = scatter(-1,-1,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1), hold on
end

for i=7:11
%b2 = plot(forcenorm(i,[11]),md_avg(i,[11]),'r*'), hold on
b1 = scatter(forcenorm(i,1),md12(i,1),'d','MarkerEdgeColor',[.5 .5 0],...
              'MarkerFaceColor',[.9 .7 0],...
              'LineWidth',1), hold on
b2 = scatter(forcenorm(i,2),md12(i,2),'d','MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
b3 = scatter(forcenorm(i,11),md12(i,11),'d','MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
b4 = scatter(forcenorm(i,12),md12(i,12),'d','MarkerEdgeColor',[.7 0 .7],...
              'MarkerFaceColor',[.9 0 .7],...
              'LineWidth',1), hold on
b5 = scatter(-1,-1,'d','MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1), hold on

end
%t1=legend([a5,b5],'H','VH')
t2=legend([a1,a2,a3,a4,a5,b5],'Baseline','Training inital','Training final','Washout','H group','VH group')
legend boxoff
xlabel('interaction force [N]')
ylabel('minimum distance [m]')
ylim([0 0.1])
xlim([0 4])
title('Partner 1: Interaction force vs minimum distance ')
axis square
% Partner 2: intforce vs minimum distance (collaboration vs compliance with task)


figure
for i=1:5
%b1 = plot(forcenorm(i,[11]),md_avg(i,[11]),'b*'), hold on
a1 = scatter(forcenorm(i,1),md21(i,1),'MarkerEdgeColor',[.7 0 .7],...
              'MarkerFaceColor',[.9 0 .7],...
              'LineWidth',1), hold on
a2 = scatter(forcenorm(i,2),md21(i,2),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
a3 = scatter(forcenorm(i,11),md21(i,11),'MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
a4 = scatter(forcenorm(i,12),md21(i,12),'MarkerEdgeColor',[.5 .5 0],...
              'MarkerFaceColor',[.9 .7 0],...
              'LineWidth',1), hold on
a5 = scatter(-1,-1,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1), hold on
end

for i=7:11
%b2 = plot(forcenorm(i,[11]),md_avg(i,[11]),'r*'), hold on
b1 = scatter(forcenorm(i,1),md21(i,1),'d','MarkerEdgeColor',[.5 .5 0],...
              'MarkerFaceColor',[.9 .7 0],...
              'LineWidth',1), hold on
b2 = scatter(forcenorm(i,2),md21(i,2),'d','MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
b3 = scatter(forcenorm(i,11),md21(i,11),'d','MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
b4 = scatter(forcenorm(i,12),md21(i,12),'d','MarkerEdgeColor',[.7 0 .7],...
              'MarkerFaceColor',[.9 0 .7],...
              'LineWidth',1), hold on
b5 = scatter(-1,-1,'d','MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1), hold on

end
%t1=legend([a5,b5],'H','VH')
% t2=legend([a1,a2,a3,a4,a5,b5],'Baseline','Training inital','Training final','Washout','H group','VH group')
% legend boxoff
xlabel('interaction force [N]')
ylabel('minimum distance [m]')
ylim([0 0.1])
xlim([0 4])
title('Partner 2: Interaction force vs minimum distance ')
axis square

%% another representation


figure
for i=1:5
%b1 = plot(forcenorm(i,[11]),md_avg(i,[11]),'b*'), hold on

a1 = scatter(forcenorm(i,2),md12(i,2),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
a3 = scatter(forcenorm(i,2),md21(i,2),'d','MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
% a5 = scatter(-1,-1,'MarkerEdgeColor',[0 0 0],...
%               'MarkerFaceColor',[1 1 1],...
%               'LineWidth',1), hold on
end

for i=7:11
%b2 = plot(forcenorm(i,[11]),md_avg(i,[11]),'r*'), hold on

b2 = scatter(forcenorm(i,2),md12(i,2),'MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
b3 = scatter(forcenorm(i,2),md21(i,2),'d','MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on

% b5 = scatter(-1,-1,'d','MarkerEdgeColor',[0 0 0],...
%               'MarkerFaceColor',[1 1 1],...
%               'LineWidth',1), hold on

end
%t1=legend([a5,b5],'H','VH')
t2=legend([a1,a3,b2,b3],'Partner 1, H','Partner 2, H','Partner 1, VH','Partner 2, VH')
legend boxoff
xlabel('interaction force [N]')
ylabel('minimum distance [m]')
ylim([0 0.1])
xlim([0 4])
title('Early: Interaction force vs minimum distance ')
axis square
% Partner 2: intforce vs minimum distance (collaboration vs compliance with task)


figure
for i=1:5
%b1 = plot(forcenorm(i,[11]),md_avg(i,[11]),'b*'), hold on

a2 = scatter(forcenorm(i,11),md12(i,11),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
a3 = scatter(forcenorm(i,11),md21(i,11),'d','MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .9 .7],...
              'LineWidth',1), hold on
% a5 = scatter(-1,-1,'MarkerEdgeColor',[0 0 0],...
%               'MarkerFaceColor',[1 1 1],...
%               'LineWidth',1), hold on
end

for i=7:11
%b2 = plot(forcenorm(i,[11]),md_avg(i,[11]),'r*'), hold on

b2 = scatter(forcenorm(i,11),md12(i,11),'MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on
b3 = scatter(forcenorm(i,11),md21(i,11),'d','MarkerEdgeColor',[.5 0 .5],...
              'MarkerFaceColor',[.7 0 .9],...
              'LineWidth',1), hold on

% b5 = scatter(-1,-1,'d','MarkerEdgeColor',[0 0 0],...
%               'MarkerFaceColor',[1 1 1],...
%               'LineWidth',1), hold on

end
%t1=legend([a5,b5],'H','VH')
% t2=legend([a2,a3,b2,b3],'Partner 1, H','Partner 2, H','Partner 1, VH','Partner 2, VH')
% legend boxoff
xlabel('interaction force [N]')
ylabel('minimum distance [m]')
ylim([0 0.1])
xlim([0 4])
title('Late: Interaction force vs minimum distance ')
axis square

%% Power calculations..
% 
% pow_vp1 = pow2-pow;
% pow_vp2 = pow2-pow1;
figure
for i=1:5
a1 = plot(pow1(i,[2 11]),pow2(i,[2 11]),'b-','LineWidth',1),hold on
plot(pow1(i,[11]),pow2(i,[11]),'b*','LineWidth',1),hold on
end
for i=7:11
a2 = plot(pow1(i,[2 11]),pow2(i,[2 11]),'r-','LineWidth',1),hold on
plot(pow1(i,[11]),pow2(i,[11]),'r*','LineWidth',1),hold on
end
legend([a1,a2],'H','VH')
legend boxoff
xlabel('Subject 1 ')
ylabel('Subject 2 [N]')
title('Power')
% xlim([0 0.05])
% ylim([0 0.05])
axis square
box off

%%
%pow_fac= pow1-pow2;
figure
for i=1:5
a1 = plot(pow_fac(i,:),'b','LineWidth',1),hold on
%  plot(pow_fac(i,[11]),'b*','LineWidth',1),hold on
end
for i=7:11
a2 = plot(pow_fac(i,:),'r','LineWidth',1),hold on
%  plot(pow_fac(i,[11]),'r*','LineWidth',1),hold on
end
t1=legend([a1,a2],'H','VH')
legend boxoff
xlabel('Dyad no')
ylabel('Difference in power between two subjects [N]')
title('Power')
% xlim([0 0.05])
% ylim([0 0.05])
axis square
box off



%% net power  represented in error bar
%%%%%%%%%%%%%%%%% p1+p2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'pos',[100 100 250 350])
title('Net power(p1+p2)')
hold on
errorbar(2+[-0.15 0.15], mean(pow_fac(subVH,[2 11])), ...
                         std(pow_fac(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(pow_fac(subH,[2 11])), ...
                         std(pow_fac(subH,[2 11]))./sqrt(length(subH)),'k.')
 
 h=bar([mean(pow_fac(subH,2)) mean(pow_fac(subH,11)); ...
 mean(pow_fac(subVH,2)) mean(pow_fac(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('p1+p2 [w]')
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
print(gcf,'-depsc','D:\eps\netpower.eps')
%%%%%%%%%%%%%%%%% Power 1%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'pos',[100 100 250 350])
title('Power 1')
hold on
errorbar(2+[-0.15 0.15], mean(pow1(subVH,[2 11])), ...
                         std(pow1(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(pow1(subH,[2 11])), ...
                         std(pow1(subH,[2 11]))./sqrt(length(subH)),'k.')
 
 h=bar([mean(pow1(subH,2)) mean(pow1(subH,11)); ...
 mean(pow1(subVH,2)) mean(pow1(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('power 1 [w]')
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
print(gcf,'-depsc','D:\eps\pow1.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Power 2 %%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'pos',[100 100 250 350])
title('Power 2')
hold on
errorbar(2+[-0.15 0.15], mean(pow2(subVH,[2 11])), ...
                         std(pow2(subVH,[2 11]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(pow2(subH,[2 11])), ...
                         std(pow2(subH,[2 11]))./sqrt(length(subH)),'k.')
 
 h=bar([mean(pow2(subH,2)) mean(pow2(subH,11)); ...
 mean(pow2(subVH,2)) mean(pow2(subVH,11))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('group')
ylabel('power 2 [w]')
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])
print(gcf,'-depsc','D:\eps\pow2.eps')
%%
for subj=1:subjs

    figure
    set(gcf,'pos',[100 100 300 400])
    patch([1.5 11.5 11.5 1.5],[0 0 .3 .3],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    hold on;
    xlabel('epochs')
    title(sprintf('Dyad%d',subj))
    plot(pow_fac(subj,:),'b')
    xlabel('trials')
    ylabel('net power(p1+p2)[w]')
%     xlim([1 trials])
     %ylim([0 0.3])

end
%%
% for subj=1:subjs
% 
%     figure
% %     set(gcf,'pos',[100 100 300 400])
% %     patch([12.5 132.5 132.5 12.5],[0 0 4 4],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
% %     hold on;
%     xlabel('trials')
%     title(sprintf('Dyad%d',subj))
% %   plot(alltsvals_corrected(subj,:),'bo')
%     plot(pow_fac_ts(subj,:),'bo')
%     xlabel('trials')
%     ylabel('p1-p2')
% %     xlim([1 trials])
% %     ylim([0 4])
% end

%% Power calculations: subject-wise per epoch
for subj=1:subjs

    figure
    set(gcf,'pos',[100 100 300 400])
    patch([1.5 11.5 11.5 1.5],[-15 -15 20 20],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    hold on;
    xlabel('epochs')
    title(sprintf('Dyad%d',subj))
    p1=plot(pow1(subj,:),'r')
    hold on
    p2=plot(pow2(subj,:),'b')
    xlabel('epochs')
    ylabel('power [w]')
    legend([p1,p2],'sub 1','sub 2')
    legend boxoff
%     xlim([1 trials])
     %ylim([0 0.3])

end

%% Percentage of leadership
for subj=1:subjs

    figure
    set(gcf,'pos',[100 100 300 400])
%     patch([1.5 11.5 11.5 1.5],[-15 -15 20 20],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
%     hold on;
 
    
    p1=  plot(lead1(subj,:),'r')
    hold on
    p2 = plot(lead2(subj,:),'b')
   
    xlabel('epochs')
    ylabel('power [w]')
    title(sprintf('Dyad%d',subj))
%     xlim([1 trials])
     %ylim([0 0.3])

end

%% subject-wise per trial
for subj=1:subjs

    figure
    set(gcf,'pos',[100 100 300 400])
    patch([12.5 132.5 132.5 12.5],[-30 -30 40 40],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    hold on;
    title(sprintf('Dyad%d',subj))
    plot(pow1_ts(subj,:),'r')
    hold on
    plot(pow2_ts(subj,:),'b')
    hold on
    plot(pow_fac_ts(subj,:),'g')
    xlabel('trials')
    ylabel('power [w]')
    xlim([1 trials])
     %ylim([0 0.3])

end

%% Plots average per-epoch

subjs=size(pow1,1);
tsets=size(pow1,2);
allvals1 = [];
allvals2 = [];
allvals3 = [];
alltsvals1 = [];
alltsvals2 = [];
for subj=1:subjs
    allvals1 = [allvals1; pow1(subj,:)];
    allvals2 = [allvals2; pow2(subj,:)];
    allvals3 = [allvals3; pow_fac(subj,:)];
    alltsvals1 = [alltsvals1; pow1_ts(subj,:)];
    alltsvals2 = [alltsvals2; pow2_ts(subj,:)];
end

%%% all subjects per epoch
figure
set(gcf,'pos',[100 100 400 550])
 patch([2 11 11 2],[-10 -10 15 15],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
 xlabel('epochs')
 ylabel('power [w]')
 title('all subjects')
 hold on

e1=errorbar(1:tsets, mean(allvals1(subVH,:)), std(allvals1(subVH,:))./sqrt(length(subVH)),'rs-','LineWidth',1)
hold on
e2=errorbar(1:tsets, mean(allvals2(subVH,:)), std(allvals2(subVH,:))./sqrt(length(subVH)),'bs-','LineWidth',1)
e5=errorbar(1:tsets, mean(allvals3(subVH,:),1), std(allvals3(subH,:),1)./sqrt(length(subH)),'ks-','LineWidth',1)
hold on
e3=errorbar(1:tsets, mean(allvals1(subH,:),1), std(allvals1(subH,:),1)./sqrt(length(subH)),'ms-','LineWidth',1)
e4=errorbar(1:tsets, mean(allvals2(subH,:),1), std(allvals2(subH,:),1)./sqrt(length(subH)),'cs-','LineWidth',1)
e6=errorbar(1:tsets, mean(allvals3(subH,:),1), std(allvals3(subH,:),1)./sqrt(length(subH)),'gs-','LineWidth',1)



%errorbar(1:tsets, allvals1(subH,:), 0*allvals1(subH,:),'rs-')
%errorbar(1:tsets, allvals2(subH,:), 0*allvals2(subH,:),'bs-')
l1=legend([e1,e2,e3,e4,e5,e6],'VH:Subj 1','VH:Subj 2','H: Subj 1','H: Subj 2','VH: P1+P2','H: P1+P2')
%rect=[1.25, 1.25, .25, .zeros25]
set(l1,'Location','north')
legend boxoff
box off
%ylim(yrange)
print(gcf,'-depsc','D:\eps\total_per_epoch.eps')

%%% all subjects per trial
figure
set(gcf,'pos',[100 100 400 550])
patch([12.5 132.5 132.5 12.5],[-15 -15 20 20],[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
xlabel('trials')
ylabel('power [w]')
title('all subjects')
hold on

e1=errorbar(1:trials, mean(alltsvals1(subVH,:)), std(alltsvals1(subVH,:))./sqrt(length(subVH)),'rs-')
hold on
e2=errorbar(1:trials, mean(alltsvals2(subVH,:)), std(alltsvals2(subVH,:))./sqrt(length(subVH)),'bs-')
hold on
e3=errorbar(1:trials, mean(alltsvals1(subH,:),1), std(alltsvals1(subH,:),1)./sqrt(length(subH)),'ms-')
e4=errorbar(1:trials, mean(alltsvals2(subH,:),1), std(alltsvals2(subH,:),1)./sqrt(length(subH)),'cs-')



%errorbar(1:tsets, allvals1(subH,:), 0*allvals1(subH,:),'rs-')
%erorbar(1:tsets, allvals2(subH,:), 0*allvals2(subH,:),'bs-')
l1=legend([e1,e2,e3,e4],'VH:Subj 1','VH:Subj 2','H: Subj 1','H: Subj 2')
%rect=[1.25, 1.25, .25, .zeros25]
set(l1,'Location','north')
legend boxoff
box off
%ylim(yrange)
 xlim([1 trials])

print(gcf,'-depsc','D:\eps\total_per_trial.eps')


%% model analysis
model_delta_lqg=load('H_model_delta_all_LQG.txt');
model_delta_nash=load('H_model_delta_all_Nash.txt');
model_md_lqg=load('H_model_md_all_LQG.txt');
model_md_nash=load('H_model_md_all_Nash.txt');
model_intf_lqg=load('H_model_intforce_LQG.txt');
model_intf_nash=load('H_model_intforce_Nash.txt');
model_intp_lqg=load('H_model_intpower_LQG.txt');
model_intp_nash=load('H_model_intpower_Nash.txt');


figure
set(gcf,'pos',[100 100 250 300])
%title('Subject 1')
hold on
errorbar(2, mean(model_delta_lqg(:,3)), ...
                         std(model_delta_lqg(:,3))./sqrt(20),'k.')
hold on
errorbar(1, mean(model_delta_nash(:,3)), ...
                         std(model_delta_nash(:,3))./sqrt(20),'k.')

h=bar([model_delta_lqg(:,3); ...
     mean(model_delta_nash(:,3))])
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


%% Leadership index plots

% lead1 = load('Lead_index_1.txt'); % at tc1 0 shift
% lead2 = load('Lead_index_2.txt');
% 
% lead1 = load('Lead_index_11.txt'); % 300ms
% lead2 = load('Lead_index_22.txt');
lead1 = load('Lead_index_111.txt'); % 300ms
lead2 = load('Lead_index_222.txt');
% 
% lead1 = load('Lead_index_111.txt'); % 300ms to tc
% lead2 = load('Lead_index_222.txt');
% 
% lead1 = load('Lead_index_1111.txt'); % at 300ms
% lead2 = load('Lead_index_2222.txt');

% 11 11 11    12 12 12 -> lead 1
% 21 21 21    22 22 22 -> lead 2
figure
set(gcf,'pos',[100 100 250 300])
title('Via point 1')
hold on
errorbar(1+[-0.15  0.15], mean(lead1(:,[4 1])), ...
                         std(lead1(:,[4 1]))./sqrt(length(lead1(:,1))),'k.')
hold on
errorbar(2+[-0.15  0.15], mean(lead1(:,[5 2])), ...
                         std(lead1(:,[5 2]))./sqrt(length(lead1(:,1))),'k.')
hold on
errorbar(3+[-0.15  0.15], mean(lead1(:,[6 3])), ...
                         std(lead1(:,[6 3]))./sqrt(length(lead1(:,1))),'k.')

h=bar([mean(lead1(:,4)) mean(lead1(:,1)); ...
       mean(lead1(:,5)) mean(lead1(:,2)); ...
       mean(lead1(:,6)) mean(lead1(:,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'})
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14)
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
 ylim([-.2  .15])
%ylim(yrange)
legend(h,{'S_1','S_2'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])

figure
set(gcf,'pos',[100 100 250 300])
title('Via point 2')
hold on
errorbar(1+[-0.15  0.15], mean(lead2(:,[4 1])), ...
                         std(lead2(:,[4 1]))./sqrt(length(lead2(:,1))),'k.')
hold on
errorbar(2+[-0.15  0.15], mean(lead2(:,[5 2])), ...
                         std(lead2(:,[5 2]))./sqrt(length(lead2(:,1))),'k.')
hold on
errorbar(3+[-0.15  0.15], mean(lead2(:,[6 3])), ...
                         std(lead2(:,[6 3]))./sqrt(length(lead2(:,1))),'k.')

h=bar([mean(lead2(:,4)) mean(lead2(:,1)); ...
       mean(lead2(:,5)) mean(lead2(:,2)); ...
       mean(lead2(:,6)) mean(lead2(:,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'})
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14)
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
 ylim([-.2  .15])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])
box off

%% LI alternative plots
% 11 11 11    12 12 12 -> lead 1
% 21 21 21    22 22 22 -> lead 2
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{11}')
hold on
errorbar(1, nanmean(lead1(:,1)), ...
                         std(lead1(:,1))./sqrt(length(lead1(:,1))),'k.')
hold on
errorbar(2, nanmean(lead1(:,2)), ...
                         std(lead1(:,2))./sqrt(length(lead1(:,1))),'k.')
hold on
errorbar(3, nanmean(lead1(:,3)), ...
                         std(lead1(:,3))./sqrt(length(lead1(:,1))),'k.')

h=bar([nanmean(lead1(:,1)); ...
       nanmean(lead1(:,2)); ...
       nanmean(lead1(:,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'})
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14)
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
  ylim([-.3  .3])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
 set(h,'facecol',[0.9 0.9 0.9])

figure
set(gcf,'pos',[100 100 250 300])
title('LI_{12}')
hold on
errorbar(1, nanmean(lead1(:,4)), ...
                         std(lead1(:,4))./sqrt(length(lead1(:,1))),'k.')
hold on
errorbar(2, nanmean(lead1(:,5)), ...
                         std(lead1(:,5))./sqrt(length(lead1(:,1))),'k.')
hold on
errorbar(3, nanmean(lead1(:,6)), ...
                         std(lead1(:,6))./sqrt(length(lead1(:,1))),'k.')

h=bar([nanmean(lead1(:,4)) ; ...
       nanmean(lead1(:,5)) ; ...
       nanmean(lead1(:,6))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'})
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14)
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
ylim([-.3  .3])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
 set(h,'facecol',[0.9 0.9 0.9])
box off

%
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{21}')
hold on
errorbar(1, nanmean(lead2(:,1)), ...
                         std(lead2(:,1))./sqrt(length(lead2(:,1))),'k.')
hold on
errorbar(2, nanmean(lead2(:,2)), ...
                         std(lead2(:,2))./sqrt(length(lead2(:,1))),'k.')
hold on
errorbar(3, nanmean(lead2(:,3)), ...
                         std(lead2(:,3))./sqrt(length(lead2(:,1))),'k.')

h=bar([nanmean(lead2(:,1)); ...
       nanmean(lead2(:,2)); ...
       nanmean(lead2(:,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'})
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14)
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
  ylim([-.3  .3])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
 set(h,'facecol',[0.9 0.9 0.9])

figure
set(gcf,'pos',[100 100 250 300])
title('LI_{22}')
hold on
errorbar(1, nanmean(lead2(:,4)), ...
                         std(lead2(:,4))./sqrt(length(lead2(:,1))),'k.')
hold on
errorbar(2, nanmean(lead2(:,5)), ...
                         std(lead2(:,5))./sqrt(length(lead2(:,1))),'k.')
hold on
errorbar(3, nanmean(lead2(:,6)), ...
                         std(lead2(:,6))./sqrt(length(lead2(:,1))),'k.')

h=bar([nanmean(lead2(:,4)) ; ...
       nanmean(lead2(:,5)) ; ...
       nanmean(lead2(:,6))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'})
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14)
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
ylim([-.3  .3])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
 set(h,'facecol',[0.9 0.9 0.9])
box off


%% three groups
c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];

figure
set(gcf,'pos',[100 100 250 300])
title('LI_{11}')
hold on
errorbar(1+[-0.225 0 0.225], [nanmean(lead1(subH,1)) nanmean(lead1(subVH,1)) nanmean(lead1(subPV,1))], ...
                         [std(lead1(subH,1))./sqrt(length(lead1(subH,1))) std(lead1(subVH,1))./sqrt(length(lead1(subVH,1))) std(lead1(subPV,1))./sqrt(length(lead1(subPV,1)))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [nanmean(lead1(subH,2)) nanmean(lead1(subVH,2)) nanmean(lead1(subPV,2))], ...
                         [std(lead1(subH,2))./sqrt(length(lead1(subH,1))) std(lead1(subVH,2))./sqrt(length(lead1(subVH,1))) std(lead1(subPV,2))./sqrt(length(lead1(subPV,1)))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [nanmean(lead1(subH,3)) nanmean(lead1(subVH,3)) nanmean(lead1(subPV,3))], ...
                         [std(lead1(subH,3))./sqrt(length(lead1(subH,1))) std(lead1(subVH,3))./sqrt(length(lead1(subVH,1))) std(lead1(subPV,3))./sqrt(length(lead1(subPV,1)))],'k.')

h=bar([nanmean(lead1(subH,1)) nanmean(lead1(subVH,1)) nanmean(lead1(subPV,1)); ...
       nanmean(lead1(subH,2)) nanmean(lead1(subVH,2)) nanmean(lead1(subPV,2)) ; ...
       nanmean(lead1(subH,3)) nanmean(lead1(subVH,3)) nanmean(lead1(subPV,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14, 'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
  ylim([-.3  .3])
%ylim(yrange)
 legend(h,{'H','VH', 'PV'})
 legend boxoff
% set(h(1),'facecol','w')
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)
name = sprintf('LI11_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

figure
set(gcf,'pos',[100 100 250 300])
title('LI_{12}')
hold on
errorbar(1+[-0.225 0 0.225], [nanmean(lead1(subH,4)) nanmean(lead1(subVH,4)) nanmean(lead1(subPV,4))], ...
                         [std(lead1(subH,4))./sqrt(length(lead1(subH,4))) std(lead1(subVH,4))./sqrt(length(lead1(subVH,4))) std(lead1(subPV,4))./sqrt(length(lead1(subPV,4)))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [nanmean(lead1(subH,5)) nanmean(lead1(subVH,5)) nanmean(lead1(subPV,5))], ...
                         [std(lead1(subH,5))./sqrt(length(lead1(subH,5))) std(lead1(subVH,5))./sqrt(length(lead1(subVH,5))) std(lead1(subPV,5))./sqrt(length(lead1(subPV,5)))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [nanmean(lead1(subH,6)) nanmean(lead1(subVH,6)) nanmean(lead1(subPV,6))], ...
                         [std(lead1(subH,6))./sqrt(length(lead1(subH,6))) std(lead1(subVH,6))./sqrt(length(lead1(subVH,6))) std(lead1(subPV,6))./sqrt(length(lead1(subPV,6)))],'k.')

h=bar([nanmean(lead1(subH,4)) nanmean(lead1(subVH,4)) nanmean(lead1(subPV,4)); ...
       nanmean(lead1(subH,5)) nanmean(lead1(subVH,5)) nanmean(lead1(subPV,5)) ; ...
       nanmean(lead1(subH,6)) nanmean(lead1(subVH,6)) nanmean(lead1(subPV,6))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14, 'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
  ylim([-.3  .3])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)
name = sprintf('LI12_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{21}')
hold on
errorbar(1+[-0.225 0 0.225], [nanmean(lead2(subH,1)) nanmean(lead2(subVH,1)) nanmean(lead2(subPV,1))], ...
                         [std(lead2(subH,1))./sqrt(length(lead2(subH,1))) std(lead2(subVH,1))./sqrt(length(lead2(subVH,1))) std(lead2(subPV,1))./sqrt(length(lead2(subPV,1)))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [nanmean(lead2(subH,2)) nanmean(lead2(subVH,2)) nanmean(lead2(subPV,2))], ...
                         [std(lead2(subH,2))./sqrt(length(lead2(subH,1))) std(lead2(subVH,2))./sqrt(length(lead2(subVH,1))) std(lead2(subPV,1))./sqrt(length(lead2(subPV,1)))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [nanmean(lead2(subH,3)) nanmean(lead2(subVH,3)) nanmean(lead2(subPV,3))], ...
                         [std(lead2(subH,3))./sqrt(length(lead2(subH,3))) std(lead2(subVH,3))./sqrt(length(lead2(subVH,3))) std(lead2(subPV,3))./sqrt(length(lead2(subPV,3)))],'k.')

h=bar([nanmean(lead2(subH,1)) nanmean(lead2(subVH,1)) nanmean(lead2(subPV,1)); ...
       nanmean(lead2(subH,2)) nanmean(lead2(subVH,2)) nanmean(lead2(subPV,2)) ; ...
       nanmean(lead2(subH,3)) nanmean(lead2(subVH,3)) nanmean(lead2(subPV,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
  ylim([-.3  .3])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)
name = sprintf('LI21_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

figure
set(gcf,'pos',[100 100 250 300])
title('LI_{22}')
hold on
errorbar(1+[-0.225 0 0.225], [nanmean(lead2(subH,4)) nanmean(lead2(subVH,4)) nanmean(lead2(subPV,4))], ...
                         [std(lead2(subH,4))./sqrt(length(lead1(subH,4))) std(lead2(subVH,4))./sqrt(length(lead2(subVH,4))) std(lead2(subPV,4))./sqrt(length(lead2(subPV,4)))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [nanmean(lead2(subH,5)) nanmean(lead2(subVH,5)) nanmean(lead2(subPV,5))], ...
                         [std(lead2(subH,5))./sqrt(length(lead2(subH,5))) std(lead2(subVH,5))./sqrt(length(lead2(subVH,5))) std(lead2(subPV,5))./sqrt(length(lead2(subPV,5)))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [nanmean(lead2(subH,6)) nanmean(lead2(subVH,6)) nanmean(lead2(subPV,6))], ...
                         [std(lead2(subH,6))./sqrt(length(lead2(subH,6))) std(lead2(subVH,6))./sqrt(length(lead2(subVH,6))) std(lead2(subPV,6))./sqrt(length(lead2(subPV,6)))],'k.')

h=bar([nanmean(lead2(subH,4)) nanmean(lead2(subVH,4)) nanmean(lead2(subPV,4)); ...
       nanmean(lead2(subH,5)) nanmean(lead2(subVH,5)) nanmean(lead2(subPV,5)) ; ...
       nanmean(lead2(subH,6)) nanmean(lead2(subVH,6)) nanmean(lead2(subPV,6))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
  ylim([-.3  .3])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)

name = sprintf('LI22_stat');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')


%% three groups only late time
c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
col = [0.2 0.5 0.5 ; 0.5 0.2 0.5; 0.5 0.5 0.2]
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{11}')
hold on
errorbar(1, [nanmean(lead1(subH,3)) ], ...
                         [std(lead1(subH,3))./sqrt(length(lead1(subH,3))) ],'k.')
hold on
errorbar(2, [nanmean(lead1(subVH,3))], ...
                         [std(lead1(subVH,3))./sqrt(length(lead1(subVH,3)))],'k.')
hold on
errorbar(3, [ nanmean(lead1(subPV,3))], ...
                         [ std(lead1(subPV,3))./sqrt(length(lead1(subPV,3)))],'k.')

h=bar([nanmean(lead1(subH,3)); ...
        nanmean(lead1(subVH,3)); ...
        nanmean(lead1(subPV,3))])
box off
  set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14, 'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
sigstar({[1,2],[1,3],[2,3]},[0.02, 0.003,0.04])
  ylim([-.4  .4])
%ylim(yrange)
%  legend(h,{'H','VH', 'PV'})
%  legend boxoff
% set(h(1),'facecol','w')
set(h(1),'facecol','w', 'edgecolor', 'k')
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)
% set(,'facecol',col, 'edgecolor', col)
name = sprintf('LI11_stat_late');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{12}')
hold on

errorbar([1 2 3], [nanmean(lead1(subH,6)) nanmean(lead1(subVH,6)) nanmean(lead1(subPV,6))], ...
                         [std(lead1(subH,6))./sqrt(length(lead1(subH,6))) std(lead1(subVH,6))./sqrt(length(lead1(subVH,6))) std(lead1(subPV,6))./sqrt(length(lead1(subPV,6)))],'k.')

h=bar([
       nanmean(lead1(subH,6)) nanmean(lead1(subVH,6)) nanmean(lead1(subPV,6))])
box off
%  set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
 ylabel('LI (w)','fontsize',14, 'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
sigstar({[1,2],[1,3]},[0.04, 0.02])
  ylim([-.4  .4])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
% set(h(1),'facecol',c1, 'edgecolor', c1)
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)
set(h(1),'facecol','w', 'edgecolor', 'k')

name = sprintf('LI12_stat_late');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{21}')
hold on

errorbar([1 2 3], [nanmean(lead2(subH,3)) nanmean(lead2(subVH,3)) nanmean(lead2(subPV,3))], ...
                         [std(lead2(subH,3))./sqrt(length(lead2(subH,3))) std(lead2(subVH,3))./sqrt(length(lead2(subVH,3))) std(lead2(subPV,3))./sqrt(length(lead2(subPV,3)))],'k.')

h=bar([nanmean(lead2(subH,3)) nanmean(lead2(subVH,3)) nanmean(lead2(subPV,3))])
box off
  set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
% sigstar({[1,3],[2,3],[1,2]},[0.005, 0.005,0.05])
  ylim([-.4  .4])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
set(h(1),'facecol','w', 'edgecolor', 'k')
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)
name = sprintf('LI21_stat_late');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

figure
set(gcf,'pos',[100 100 250 300])
title('LI_{22}')
hold on
errorbar([1 2 3], [nanmean(lead2(subH,6)) nanmean(lead2(subVH,6)) nanmean(lead2(subPV,6))], ...
                         [std(lead2(subH,6))./sqrt(length(lead2(subH,6))) std(lead2(subVH,6))./sqrt(length(lead2(subVH,6))) std(lead2(subPV,6))./sqrt(length(lead2(subPV,6)))],'k.')

h=bar([nanmean(lead2(subH,6)) nanmean(lead2(subVH,6)) nanmean(lead2(subPV,6))])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1,3],[2,3],[1,2]},[0.005, 0.005,0.05])
  ylim([-.4  .4])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
set(h(1),'facecol','w', 'edgecolor', 'k')
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)

name = sprintf('LI22_stat_late');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')




%% three groups only late time diff colors
c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
col = [0.2 0.5 0.5 ; 0.5 0.2 0.5; 0.5 0.5 0.2]
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{11}')
hold on
errorbar(1, [nanmean(lead1(subH,3)) ], ...
                         [std(lead1(subH,3))./sqrt(length(lead1(subH,3))) ],'k.')
hold on
errorbar(2, [nanmean(lead1(subVH,3))], ...
                         [std(lead1(subVH,3))./sqrt(length(lead1(subVH,3)))],'k.')
hold on
errorbar(3, [ nanmean(lead1(subPV,3))], ...
                         [ std(lead1(subPV,3))./sqrt(length(lead1(subPV,3)))],'k.')

% h=bar([nanmean(lead1(subH,3)); ...
%         nanmean(lead1(subVH,3)); ...
%         nanmean(lead1(subPV,3))])

hb1=bar(1, [nanmean(lead1(subH,3))]);
hold on
hb2=bar(2, [nanmean(lead1(subVH,3))]);
hb3=bar(3, [nanmean(lead1(subPV,3))]);
set(gca, 'xtick',1:3,'xticklabel',{'H','VH','PV'}, 'FontName', 'Times New Roman')      
set(hb1,'facecol',c1);
set(hb2,'facecol',c2);
set(hb3,'facecol',c3);

box off
%   set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14, 'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
sigstar({[1,2],[1,3],[2,3]},[0.02, 0.003,0.04])
  ylim([-.4  .4])
%ylim(yrange)
%  legend(h,{'H','VH', 'PV'})
%  legend boxoff
% set(h(1),'facecol','w')
%set(h(1),'facecol','w', 'edgecolor', 'k')
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)
% set(,'facecol',col, 'edgecolor', col)
name = sprintf('LI11_stat_late_c');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{12}')
hold on

errorbar([1 2 3], [nanmean(lead1(subH,6)) nanmean(lead1(subVH,6)) nanmean(lead1(subPV,6))], ...
                         [std(lead1(subH,6))./sqrt(length(lead1(subH,6))) std(lead1(subVH,6))./sqrt(length(lead1(subVH,6))) std(lead1(subPV,6))./sqrt(length(lead1(subPV,6)))],'k.')

hb1=bar(1, [nanmean(lead1(subH,6))]);
hold on
hb2=bar(2, [nanmean(lead1(subVH,6))]);
hb3=bar(3, [nanmean(lead1(subPV,6))]);
set(gca, 'xtick',1:3,'xticklabel',{'H','VH','PV'}, 'FontName', 'Times New Roman')      
set(hb1,'facecol',c1);
set(hb2,'facecol',c2);
set(hb3,'facecol',c3);
box off
%  set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
 ylabel('LI (w)','fontsize',14, 'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
sigstar({[1,2],[1,3]},[0.04, 0.02])
  ylim([-.4  .4])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
% set(h(1),'facecol',c1, 'edgecolor', c1)
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)
%set(h(1),'facecol','w', 'edgecolor', 'k')

name = sprintf('LI12_stat_late_c');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%
figure
set(gcf,'pos',[100 100 250 300])
title('LI_{21}')
hold on

errorbar([1 2 3], [nanmean(lead2(subH,3)) nanmean(lead2(subVH,3)) nanmean(lead2(subPV,3))], ...
                         [std(lead2(subH,3))./sqrt(length(lead2(subH,3))) std(lead2(subVH,3))./sqrt(length(lead2(subVH,3))) std(lead2(subPV,3))./sqrt(length(lead2(subPV,3)))],'k.')

hb1=bar(1, [nanmean(lead2(subH,3))]);
hold on
hb2=bar(2, [nanmean(lead2(subVH,3))]);
hb3=bar(3, [nanmean(lead2(subPV,3))]);
set(gca, 'xtick',1:3,'xticklabel',{'H','VH','PV'}, 'FontName', 'Times New Roman')      
set(hb1,'facecol',c1);
set(hb2,'facecol',c2);
set(hb3,'facecol',c3);
box off
  set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
% sigstar({[1,3],[2,3],[1,2]},[0.005, 0.005,0.05])
  ylim([-.4  .4])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
%set(h(1),'facecol','w', 'edgecolor', 'k')
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)
name = sprintf('LI21_stat_late_c');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

figure
set(gcf,'pos',[100 100 250 300])
title('LI_{22}')
hold on
errorbar([1 2 3], [nanmean(lead2(subH,6)) nanmean(lead2(subVH,6)) nanmean(lead2(subPV,6))], ...
                         [std(lead2(subH,6))./sqrt(length(lead2(subH,6))) std(lead2(subVH,6))./sqrt(length(lead2(subVH,6))) std(lead2(subPV,6))./sqrt(length(lead2(subPV,6)))],'k.')

hb1=bar(1, [nanmean(lead2(subH,6))]);
hold on
hb2=bar(2, [nanmean(lead2(subVH,6))]);
hb3=bar(3, [nanmean(lead2(subPV,6))]);
set(gca, 'xtick',1:3,'xticklabel',{'H','VH','PV'}, 'FontName', 'Times New Roman')      
set(hb1,'facecol',c1);
set(hb2,'facecol',c2);
set(hb3,'facecol',c3);
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'H', 'VH', 'PV'}, 'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('LI (w)','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1,3],[2,3],[1,2]},[0.005, 0.005,0.05])
  ylim([-.4  .4])
%ylim(yrange)
% legend(h,{'S_1','S_2'})
% legend boxoff
% set(h(1),'facecol','w')
%set(h(1),'facecol','w', 'edgecolor', 'k')
% set(h(2),'facecol',c2, 'edgecolor', c2)
% set(h(3),'facecol',c3, 'edgecolor', c3)

name = sprintf('LI22_stat_late_c');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')

%% learning traces like plot

lead_vp1 = -lead1(:,1:3) + lead2(:, 1:3);
lead_vp2 = -lead2(:,4:6) + lead1(:, 4:6);
lead_vp1(isnan(lead_vp1))=0;
lead_vp2(isnan(lead_vp2))=0;

% model parameters
nash_lead_vp1=load('leadindex_vp1_1.txt');
nash_lead_vp2=load('leadindex_vp2_1.txt');
lqg_lead_vp1=load('leadindex_vp1_7.txt');
lqg_lead_vp2=load('leadindex_vp2_7.txt');

 lead_vp1(6,:) = nan;
 lead_vp2(6,:) = nan;

col1 =[.4 .5 .7];
col2 =[.7 .3 .2];
c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];

figure
set(gcf,'pos',[100 100 300 300])
for i=1:5
% a1 = plot(lead_vp1(i,[1 3]),lead_vp2(i,[1 3]),'-','LineWidth',1, 'Color',col1),hold on
a1=plot(lead_vp1(i,[3]),lead_vp2(i,[3]),'.','MarkerSize',20,'MarkerEdgeColor', c1, 'MarkerFaceColor', 'w')
 hold on
 plot([0 0],[-0.2 0.7],'k--')
 
 labelpoints(lead_vp1(i,[3]), lead_vp2(i,[3]), i, 'Color', c1, 'FontSize', 8)
end
for i=6:11

% a2 = plot(lead_vp1(i,[1 3]),lead_vp2(i,[1 3]),'-','LineWidth',1, 'Color', col2),hold on
a2=plot(lead_vp1(i,[3]),lead_vp2(i,[3]),'.','MarkerSize',20,'MarkerEdgeColor', c2, 'MarkerFaceColor', 'w')
hold on
plot([-0.2 0.7],[0 0],'k--')

labelpoints(lead_vp1(i,[3]), lead_vp2(i,[3]), i, 'Color', c2,'FontSize', 8)
   
end

for i=12:16

% a2 = plot(lead_vp1(i,[1 3]),lead_vp2(i,[1 3]),'-','LineWidth',1, 'Color', col2),hold on
a3=plot(lead_vp1(i,[3]),lead_vp2(i,[3]),'.','MarkerSize',20,'MarkerEdgeColor', c3, 'MarkerFaceColor', 'w')
hold on
plot([-0.2 0.7],[0 0],'k--')

labelpoints(lead_vp1(i,[3]), lead_vp2(i,[3]), i, 'Color', c3,'FontSize', 8)
   
end

xlabel('LI_{21}-LI_{11} [w]','fontsize',14,'FontName', 'Times New Roman')
ylabel('LI_{12}-LI_{22} [w]','fontsize',14,'FontName', 'Times New Roman')
% title('mean of intervel 300ms before TC')
 xlim([-0.2 0.7])
ylim([-.2 .7])
axis square
box off
col1 =[.9 .9 .1];
col2 =[.1 .8 .1];

% another representation (Nash vs. LQG)
% LQG
 hold on
set(gcf,'pos',[100 100 300 300])
 plot(mean(lqg_lead_vp1,1),mean(lqg_lead_vp2,1),'.','MarkerSize',30,'MarkerEdgeColor', col1, 'MarkerFaceColor', 'w')
plot(mean(nash_lead_vp1,1),mean(nash_lead_vp2,1),'.','MarkerSize',30, 'MarkerEdgeColor', col2, 'MarkerFaceColor', 'w')
 
plot([0 0],[-0.2 0.7],'k--'), hold on
plot([-0.2 0.7],[0 0],'k--')
text(0.2,0.30,'No partner' ,'Color',col1,'FontSize',10,'FontName', 'Times New Roman')
text(-0.13,-0.06,'Nash' ,'Color',col2,'FontSize',10,'FontName', 'Times New Roman')
% legend([a1,a2],'2×LQG','Nash') 
% legend boxoff
% ylim([0 0.1])
% xlim([0 4])
%  xlabel('Difference in power at VP_{1} [w]','fontsize',14,'FontName', 'Times New Roman')
%  ylabel('Difference in power at VP_{2} [N]','fontsize',14,'FontName', 'Times New Roman')
xlim([-0.2 0.7])
ylim([-0.2 0.7])
% title('LQG: interaction force vs minimum distance ')
legend([a1,a2,a3],'H','VH','PV', 'Location', 'northeast')
legend boxoff
axis square
box off

name = sprintf('leadindex_expvssim');
fullFileName = fullfile(folder, name);
print(gcf,fullFileName,'-depsc', '-r300')


% stat
lead_1 = [lead_vp1(1:5,:); lead_vp1(7:11,:)];
lead_2 = [lead_vp2(1:5,:); lead_vp2(7:11,:)];
% l1(6,:)=[];
% l2(6,:)=[];

%% ttest at VP1
[h1, p1, ci1, stat1]=ttest2(lead_vp1(1:5,3), lead_vp1(7:11,3));
[h2, p2, ci2, stat2]=ttest2(lead_vp2(1:5,3), lead_vp2(7:11,3));

lead_vp = [lead_1_2_H; lead_1_2_VH];
save('lead_vp.txt','lead_vp', '-ascii')


%% H and VH

groupH = 1:5;
groupVH = 7:11;  % 6 removed

figure
set(gcf,'pos',[100 100 250 300])
title('Via point 1, H group','FontName', 'Times New Roman','fontsize',14)
hold on
errorbar(1+[-0.15  0.15], mean(lead1(groupH,[4 1])), ...
                         std(lead1(groupH,[4 1]))./sqrt(length(groupH)),'k.')
hold on
errorbar(2+[-0.15  0.15], mean(lead1(groupH,[5 2])), ...
                         std(lead1(groupH,[5 2]))./sqrt(length(groupH)),'k.')
hold on
errorbar(3+[-0.15  0.15], mean(lead1(groupH,[6 3])), ...
                         std(lead1(groupH,[6 3]))./sqrt(length(groupH)),'k.')

h=bar([mean(lead1(groupH,4)) mean(lead1(groupH,1)); ...
       mean(lead1(groupH,5)) mean(lead1(groupH,2)); ...
       mean(lead1(groupH,6)) mean(lead1(groupH,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'},'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
%  ylim([-.2  .15])
%ylim(yrange)
legend(h,{'S_1','S_2'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])

figure
set(gcf,'pos',[100 100 250 300])
title('Via point 1, VH group','FontName', 'Times New Roman','fontsize',14)
hold on
errorbar(1+[-0.15  0.15], mean(lead1(groupVH,[4 1])), ...
                         std(lead1(groupVH,[4 1]))./sqrt(length(groupVH)),'k.')
hold on
errorbar(2+[-0.15  0.15], mean(lead1(groupVH,[5 2])), ...
                         std(lead1(groupVH,[5 2]))./sqrt(length(groupVH)),'k.')
hold on
errorbar(3+[-0.15  0.15], mean(lead1(groupVH,[6 3])), ...
                         std(lead1(groupVH,[6 3]))./sqrt(length(groupVH)),'k.')

h=bar([mean(lead1(groupVH,4)) mean(lead1(groupVH,1)); ...
       mean(lead1(groupVH,5)) mean(lead1(groupVH,2)); ...
       mean(lead1(groupVH,6)) mean(lead1(groupVH,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'},'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14,'FontName', 'Times New Roman','fontsize',14)
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
%  ylim([-.2  .15])
%ylim(yrange)
legend(h,{'S_1','S_2'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])


%
figure
set(gcf,'pos',[100 100 250 300])
title('Via point 2, H group','FontName', 'Times New Roman','fontsize',14)
hold on
errorbar(1+[-0.15  0.15], mean(lead2(groupH,[4 1])), ...
                         std(lead2(groupH,[4 1]))./sqrt(length(groupH)),'k.')
hold on
errorbar(2+[-0.15  0.15], mean(lead2(groupH,[5 2])), ...
                         std(lead2(groupH,[5 2]))./sqrt(length(groupH)),'k.')
hold on
errorbar(3+[-0.15  0.15], mean(lead2(groupH,[6 3])), ...
                         std(lead2(groupH,[6 3]))./sqrt(length(groupH)),'k.')

h=bar([mean(lead2(groupH,4)) mean(lead2(groupH,1)); ...
       mean(lead2(groupH,5)) mean(lead2(groupH,2)); ...
       mean(lead2(groupH,6)) mean(lead2(groupH,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'},'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
%  ylim([-.2  .15])
%ylim(yrange)
legend(h,{'S_1','S_2'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])

figure
set(gcf,'pos',[100 100 250 300])
title('Via point 2, VH group','FontName', 'Times New Roman','fontsize',14)
hold on
errorbar(1+[-0.15  0.15], mean(lead2(groupVH,[4 1])), ...
                         std(lead2(groupVH,[4 1]))./sqrt(length(groupVH)),'k.')
hold on
errorbar(2+[-0.15  0.15], mean(lead2(groupVH,[5 2])), ...
                         std(lead2(groupVH,[5 2]))./sqrt(length(groupVH)),'k.')
hold on
errorbar(3+[-0.15  0.15], mean(lead2(groupVH,[6 3])), ...
                         std(lead2(groupVH,[6 3]))./sqrt(length(groupVH)),'k.')

h=bar([mean(lead2(groupVH,4)) mean(lead2(groupVH,1)); ...
       mean(lead2(groupVH,5)) mean(lead2(groupVH,2)); ...
       mean(lead2(groupVH,6)) mean(lead2(groupVH,3))])
box off
 set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Mid', 'Late'},'FontName', 'Times New Roman')
%  xlabel('group','fontsize',14)
 ylabel('Leadership Index','fontsize',14,'FontName', 'Times New Roman')
%  sigstar({[1-0.225,1+0.225],[2-0.225,2+0.225],[1,2]},[0.005, 0.005,0.05])
%  ylim([-.2  .15])
%ylim(yrange)
legend(h,{'S_1','S_2'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.9 0.9 0.9])



%% vp1 and vp2

figure
set(gcf,'pos',[100 100 250 350])
title('VP_{1}')
hold on
errorbar(2+[-0.15 0.15], mean(lead_vp1(subVH,[1 3])), ...
                         std(lead_vp1(subVH,[1 3]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(lead_vp1(subH,[1 3])), ...
                         std(lead_vp1(subH,[1 3]))./sqrt(length(subH)),'k.')

h=bar([mean(lead_vp1(subH,1)) mean(lead_vp1(subH,3)); ...
     mean(lead_vp1(subVH,1)) mean(lead_vp1(subVH,3))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('Group')
ylabel('Power [N]')
 ylim([-0.4 .4])
%ylim(yrange)
legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])


figure
set(gcf,'pos',[100 100 250 350])
title('VP_{2}')
hold on
errorbar(2+[-0.15 0.15], mean(lead_vp2(subVH,[1 3])), ...
                         std(lead_vp2(subVH,[1 3]))./sqrt(length(subVH)),'k.')
hold on
errorbar(1+[-0.15 0.15], mean(lead_vp2(subH,[1 3])), ...
                         std(lead_vp2(subH,[1 3]))./sqrt(length(subH)),'k.')

h=bar([mean(lead_vp2(subH,1)) mean(lead_vp2(subH,3)); ...
     mean(lead_vp2(subVH,1)) mean(lead_vp2(subVH,3))])
box off
set(gca,'xtick',[1 2],'xticklabel',{'H', 'VH'})
xlabel('Group')
ylabel('Power [N]')
 ylim([-0.4 .4])
%ylim(yrange)
%legend(h,{'early','late'})
legend boxoff
set(h(1),'facecol','w')
set(h(2),'facecol',[0.8 0.8 0.8])

%% Stats analysis on lead index(epoch wise)
group = zeros(subjs,1);
StatMatrix1 = [group, lead1(:,[1 2 3]),lead1(:,[4 5 6])];
StatMatrix1(subH,1)=1;
StatMatrix1(subVH,1)=2;
StatMatrix1(subPV,1)=3;
save('lead_par1_VP1_VP2.txt','StatMatrix1', '-ascii')
group = zeros(subjs,1);
StatMatrix2 = [group, lead2(:,[1 2 3]),lead2(:,[4 5 6])];
StatMatrix2(subVH,1)=2;
StatMatrix2(subPV,1)=3;
save('lead_par2_VP1_VP2.txt','StatMatrix2', '-ascii')

%% for anova R
subjects=[(1:15)';(1:15)';(1:15)'];
exposure = [ones(1, 15)'; 2*ones(1,15)'; 3*ones(1,15)'];

DV11 = [lead1(subH,1);lead1(subVH, 1);lead1(subPV, 1);lead1(subH,2); lead1(subVH, 2);lead1(subPV, 2);lead1(subH,3); lead1(subVH, 3);lead1(subPV, 3)];
DV22 = [lead2(subH,4); lead2(subVH, 4);lead2(subPV, 4);lead2(subH,5); lead2(subVH, 5);lead2(subPV, 5);lead2(subH,6); lead2(subVH, 6);lead2(subPV, 6)];
DV12 = [lead1(subH,4);lead1(subVH, 4);lead1(subPV, 4);lead1(subH,5); lead1(subVH, 5);lead1(subPV, 5);lead1(subH,6); lead1(subVH, 6);lead1(subPV, 6)];
DV21 = [lead2(subH,1); lead2(subVH, 1);lead2(subPV, 1);lead2(subH,2); lead2(subVH, 2);lead2(subPV, 2);lead2(subH,3); lead2(subVH, 3);lead2(subPV, 3)];
HVH = [ones(1, 5)'; 2*ones(1,5)';3*ones(1,5)'];
group1 = [HVH;HVH;HVH];

for_anova11= [DV11 group1 exposure subjects];
for_anova22= [DV22 group1 exposure subjects];
for_anova21= [DV21 group1 exposure subjects];
for_anova12= [DV12 group1 exposure subjects];

save('lead11_anova.txt','for_anova11', '-ascii')
save('lead22_anova.txt','for_anova22', '-ascii')
save('lead12_anova.txt','for_anova12', '-ascii')
save('lead21_anova.txt','for_anova21', '-ascii')

%% figure 4

c1= [0.2 0.5 0.5];
c2= [0.5 0.2 0.5];
c3= [0.5 0.5 0.2];
% name = sprintf('score1_stat');
% fullFileName = fullfile(folder, name);
figure
title('VP_{1}')
% title('MD_{12}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
% title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(1+[-0.225 0 0.225], [mean(lead_vp1(subH,1)) mean(lead_vp1(subVH,1)) mean(lead_vp1(subPV,1))], ...
                         [std(lead_vp1(subH,1))./sqrt(length(subH)) std(lead_vp1(subVH,1))./sqrt(length(subVH)) std(lead_vp1(subPV,1))./sqrt(length(subPV))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [mean(lead_vp1(subH,2)) mean(lead_vp1(subVH,2)) mean(lead_vp1(subPV,2))], ...
                         [std(lead_vp1(subH,2))./sqrt(length(subH)) std(lead_vp1(subVH,2))./sqrt(length(subVH)) std(lead_vp1(subPV,2))./sqrt(length(subPV))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [mean(lead_vp1(subH,3)) mean(lead_vp1(subVH,3)) mean(lead_vp1(subPV,3))], ...
                         [std(lead_vp1(subH,3))./sqrt(length(subH)) std(lead_vp1(subVH,3))./sqrt(length(subVH)) std(lead_vp1(subPV,3))./sqrt(length(subPV))],'k.')
 hold on
h=bar([mean(lead_vp1(subH,1)) mean(lead_vp1(subVH,1)) mean(lead_vp1(subPV,1)) ; ...
      mean(lead_vp1(subH,2)) mean(lead_vp1(subVH,2)) mean(lead_vp1(subPV,2));
      mean(lead_vp1(subH,3)) mean(lead_vp1(subVH,3)) mean(lead_vp1(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
%   sigstar({[2-0.225,2+0.225], [3,3-0.225], [3-0.225,3+0.225]},[0.05, 0.05, 0.05])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'},'FontName', 'Times New Roman')
% xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Difference in power [w]','fontsize',14,'FontName', 'Times New Roman')
ylim([-0.4 0.4])

legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)

figure
title('VP_{2}')
% title('MD_{12}','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
% title('Subject 1','fontsize', 14,'FontName', 'Times New Roman', 'FontWeight','Normal')
hold on
set(gcf,'pos',[100 100 250 300])
errorbar(1+[-0.225 0 0.225], [mean(lead_vp2(subH,1)) mean(lead_vp2(subVH,1)) mean(lead_vp2(subPV,1))], ...
                         [std(lead_vp2(subH,1))./sqrt(length(subH)) std(lead_vp2(subVH,1))./sqrt(length(subVH)) std(lead_vp2(subPV,1))./sqrt(length(subPV))],'k.')
hold on
errorbar(2+[-0.225 0 0.225], [mean(lead_vp2(subH,2)) mean(lead_vp2(subVH,2)) mean(lead_vp2(subPV,2))], ...
                         [std(lead_vp2(subH,2))./sqrt(length(subH)) std(lead_vp2(subVH,2))./sqrt(length(subVH)) std(lead_vp2(subPV,2))./sqrt(length(subPV))],'k.')
hold on
errorbar(3+[-0.225 0 0.225], [mean(lead_vp2(subH,3)) mean(lead_vp2(subVH,3)) mean(lead_vp2(subPV,3))], ...
                         [std(lead_vp2(subH,3))./sqrt(length(subH)) std(lead_vp2(subVH,3))./sqrt(length(subVH)) std(lead_vp2(subPV,3))./sqrt(length(subPV))],'k.')
 hold on
h=bar([mean(lead_vp2(subH,1)) mean(lead_vp2(subVH,1)) mean(lead_vp2(subPV,1)) ; ...
      mean(lead_vp2(subH,2)) mean(lead_vp2(subVH,2)) mean(lead_vp2(subPV,2));
      mean(lead_vp2(subH,3)) mean(lead_vp2(subVH,3)) mean(lead_vp2(subPV,3))])
 
 
%  sigstar({[1-0.15,1+0.15],[2-0.15,2+0.15], [1+0.15,2+0.15]},[0.0005, 0.0005, 0.005])
%   sigstar({[2-0.225,2+0.225], [3,3-0.225], [3-0.225,3+0.225]},[0.05, 0.05, 0.05])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'Early', 'Middle', 'Late'},'FontName', 'Times New Roman')
% xlabel('Group', 'fontsize', 14,'FontName', 'Times New Roman')
ylabel('Difference in power [w]','fontsize',14,'FontName', 'Times New Roman')

  ylim([-0.4 0.4])
legend(h,{'H','VH','PV'}, 'Location', 'southeast')
legend boxoff
% set(h(2),'facecol','w')
% set(h(3),'facecol',[0.9 0.9 0.9])
% set(h(4),'facecol',[0.4 0.4 0.4])
% set(h(1),'facecol',[0.4 0.9 0.4])
% set(h(5),'facecol',[.7 .7 .1])
set(h(1),'facecol',c1, 'edgecolor', c1)
set(h(2),'facecol',c2, 'edgecolor', c2)
set(h(3),'facecol',c3, 'edgecolor', c3)


