%% ANOVA
% by Giulia Sedda

clear all
close all

%% General settings
datadir = '../../Data/';
resdir = [datadir, 'results/'];

n_H = 10;
n_VH = 10;
subjects_TE = n_H + n_VH;
H = 1:n_H;
VH = n_H+1:n_H+n_VH;
% dyads_h = 1:(subjects_TE-n_VH*2)/2;
% dyads_h = 1:(n_VH*2-subjects_TE)/2;
% dyads_a = size(dyads_h')+1:size(dyads_h')+n_VH;

dyads_H = [1 2 3 4 5];
dyads_VH = [6 7 8 9 10];
dyads = (n_H+n_VH)/2;

phases_d = {'baseline','training','aftereffect'}; 
order_d = {1,2:11,12:13};  % target sets in experiment 
tsets_d = 13;
trials = 12;

%% Minimum VP distances settings

load([resdir, 'md11_statmat.txt']);
load([resdir, 'md11_statmat_nofor.txt']);
load([resdir, 'md12_statmat.txt']);
load([resdir, 'md12_statmat_nofor.txt']);
load([resdir, 'md21_statmat.txt']);
load([resdir, 'md21_statmat_nofor.txt']);
load([resdir, 'md22_statmat.txt']);
load([resdir, 'md22_statmat_nofor.txt']);

% md11_st_nofor_d = md11_statmat_nofor;
% md12_st_nofor_d = md12_statmat_nofor;
% md22_st_nofor_d = md22_statmat_nofor;
% md21_st_nofor_d = md21_statmat_nofor;
% 
% % load per-trial data
% load ([resdir,'md11.mat']);
% md11_d = (ind_ts);
% load ([resdir,'md12.mat']);
% md12_d = (ind_ts);
% load ([resdir,'md21.mat']);
% md21_d = (ind_ts);
% load ([resdir,'md22.mat']);
% md22_d = (ind_ts);
% load([resdir, 'forces.mat']);

%% Effort settings

% load per-epoch effort data
load([resdir, 'eff1_statmat.txt']);
load([resdir, 'eff1_statmat_catch.txt']);
load([resdir, 'eff2_statmat.txt']);
load([resdir, 'eff2_statmat_catch.txt']);
load([resdir, 'R_statmat.txt']);
load([resdir, 'R_statmat_catch.txt']);

% load per-trial data
load ([resdir,'eff1.mat']);
eff1_ts_d = sqrt(ind_ts);
load ([resdir,'eff2.mat']);
eff2_ts_d = sqrt(ind_ts);
load([resdir, 'forces.mat']);
load ([resdir,'R.mat']);
R_ts_d = sqrt(ind_ts);

eff1_ts_d(find(forces==0))= nan;
eff2_ts_d(find(forces==0))= nan;

%% Delta crossing time settings

load([resdir, 'tc11_statmat.txt']); % S1 crossing VP1
load([resdir, 'tc12_statmat.txt']); % S1 crossing VP2
load([resdir, 'tc21_statmat.txt']); % S2 crossing VP1
load([resdir, 'tc22_statmat.txt']); % S2 crossing VP2

load([resdir, 'tc11_statmat_nofor.txt']); % crossing time
load([resdir, 'tc12_statmat_nofor.txt']);
load([resdir, 'tc21_statmat_nofor.txt']);
load([resdir, 'tc22_statmat_nofor.txt']);

tc11_statmat_nofor_d = tc11_statmat_nofor;
tc12_statmat_nofor_d = tc12_statmat_nofor;
tc21_statmat_nofor_d = tc21_statmat_nofor;
tc22_statmat_nofor_d = tc22_statmat_nofor;

load([resdir, 'ts1_statmat.txt']); % start time: when each subject start to move the robot
load([resdir, 'ts2_statmat.txt']);
load([resdir, 'te1_statmat.txt']); % end time: when each subject end the movement
load([resdir, 'te2_statmat.txt']);

% ts1d = ts1_statmat(:,2:end);
% ts2d = ts2_statmat(:,2:end);
% te1d = te1_statmat(:,2:end);
% te2d = te2_statmat(:,2:end);

% % load per-trial data
% load ([resdir,'tc11.mat']);
% tc11_ts_d = (ind_ts);
% load ([resdir,'tc12.mat']);
% tc12_ts_d = (ind_ts);
% load ([resdir,'tc21.mat']);
% tc21_ts_d = (ind_ts);
% load ([resdir,'tc22.mat']);
% tc22_ts_d = (ind_ts);
% 
% % Load per-trial data
% load ([resdir,'ts1.mat']);
% ts1_ts_d = ind_ts;
% load ([resdir,'ts2.mat']);
% ts2_ts_d = ind_ts;
% load ([resdir,'te1.mat']);
% te1_ts_d = ind_ts;
% load ([resdir,'te2.mat']);
% te2_ts_d = ind_ts;

%% Power settings

% load per-epoch data
load([resdir, 'pow1_statmat.txt']);
load([resdir, 'pow1_statmat_catch.txt']);
load([resdir, 'pow2_statmat.txt']);
load([resdir, 'pow2_statmat_catch.txt']);

% load per-trial data
% load ([resdir,'pow1.mat']);
% pow1_ts_d = sqrt(ind_ts);
% load ([resdir,'pow2.mat']);
% pow2_ts_d = sqrt(ind_ts);
% load([resdir, 'forces.mat']);
% 
% pow1_ts_d(find(forces==0)) = nan;
% pow2_ts_d(find(forces==0)) = nan;

pow1_st_d = pow1_statmat;
pow2_st_d = pow2_statmat;



% Indicators for subjects
%% Minimum VP distances per-subject, per-epoch
md11_statmat(7,:) = [];
md12_statmat(7,:) = [];
md22_statmat(7,:) = [];
md21_statmat(7,:) = [];

md11_d_t = md11_statmat(:,2:end);
md12_d_t = md12_statmat(:,2:end);
md21_d_t = md21_statmat(:,2:end);
md22_d_t = md22_statmat(:,2:end);



md12c_d_t = md12_statmat_nofor(:,2:end);
md21c_d_t = md21_statmat_nofor(:,2:end);

md_avg = 0.5*(md12_d_t+md21_d_t);

label = {'C','A'};

%% Own VIA-point
MD_ii = [md22_d_t; md11_d_t]; 
group_MD_ii = [md22_statmat(:,1); md11_statmat(:,1)];
% for i=dyads_a
%     group_MD_ii(i,1) = 1;
% end
subjects_MD_ii = (1:length(MD_ii))';
t_MD_ii = table(label(group_MD_ii)',MD_ii(:,2),MD_ii(:,11),'VariableNames',{'group','pre','post'});
time_MD_ii = table([1 2]','VariableNames',{'Time'});

rm_MD_ii = fitrm(t_MD_ii,'post-pre~group','WithinDesign',time_MD_ii);
Apre_MD_ii = [ MD_ii(:,2) group_MD_ii ones(length(group_MD_ii),1) subjects_MD_ii];
Amid_MD_ii = [ MD_ii(:,6) group_MD_ii 2*ones(length(group_MD_ii),1) subjects_MD_ii];
Apos_MD_ii = [ MD_ii(:,11) group_MD_ii 3*ones(length(group_MD_ii),1) subjects_MD_ii];
A_MD_ii = [Apre_MD_ii; Amid_MD_ii; Apos_MD_ii];
save('ownvia_fourcolumn.txt','A_MD_ii', '-ascii')
BWAOV2(A_MD_ii,0.05);

%% partner VP
MD_ij = [md21_d_t; md12_d_t]; 
group_MD_ij = [md22_statmat(:,1); md11_statmat(:,1)];
%gg = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2];
subjects_MD_ij = (1:length(MD_ij))';
t_MD_ij = table(label(group_MD_ij)',MD_ij(:,2),MD_ij(:,11),'VariableNames',{'group','pre','post'});
time_MD_ij = table([1 2]','VariableNames',{'Time'});

rm_MD_ij = fitrm(t_MD_ij,'post-pre~group','WithinDesign',time_MD_ij);
Apre_MD_ij = [ MD_ij(:,2) group_MD_ij ones(length(group_MD_ij),1) subjects_MD_ij];
Amid_MD_ij = [ MD_ij(:,6) group_MD_ij 2*ones(length(group_MD_ij),1) subjects_MD_ij];
Apos_MD_ij = [ MD_ij(:,11) group_MD_ij 3*ones(length(group_MD_ij),1) subjects_MD_ij];
A_MD_ij = [Apre_MD_ij; Amid_MD_ij; Apos_MD_ij];
save('parvia_fourcolumn.txt','A_MD_ij', '-ascii')
BWAOV2(A_MD_ij);


%% partner VP subject 1 betwwen H and VH
MD_ij = [md12_d_t]; 
group_MD = [md11_statmat(:,1)];
subjects_MD_ij = (1:length(MD_ij(:,1)))';
% t_MD_ij = table(label(group_MD_ij)',MD_ij(:,2),MD_ij(:,11),'VariableNames',{'group','pre','post'});
time_MD_ij = table([1 2]','VariableNames',{'Time'});

rm_MD_ij = fitrm(t_MD_ij,'post-pre~group','WithinDesign',time_MD_ij);
Apre_MD_ij = [ MD_ij(:,2) group_MD ones(length(group_MD),1) subjects_MD_ij];
Amid_MD_ij = [ MD_ij(:,6) group_MD 2*ones(length(group_MD),1) subjects_MD_ij];
Apos_MD_ij = [ MD_ij(:,11) group_MD 3*ones(length(group_MD),1) subjects_MD_ij];
A_MD_ij = [Apre_MD_ij; Amid_MD_ij; Apos_MD_ij];
save('subject_1_fourcolumn.txt','A_MD_ij', '-ascii')
BWAOV2(A_MD_ij);

%% partner VP subject 2 betwwen H and VH
MD_ij = [md21_d_t]; 
group_MD = [md22_statmat(:,1)];
subjects_MD_ij = (1:length(MD_ij(:,1)))';
%t_MD_ij = table(label(group_MD_ij)',MD_ij(:,2),MD_ij(:,11),'VariableNames',{'group','pre','post'});
time_MD_ij = table([1 2]','VariableNames',{'Time'});

rm_MD_ij = fitrm(t_MD_ij,'post-pre~group','WithinDesign',time_MD_ij);
Apre_MD_ij = [ MD_ij(:,2) group_MD ones(length(group_MD),1) subjects_MD_ij];
Amid_MD_ij = [ MD_ij(:,6) group_MD 2*ones(length(group_MD),1) subjects_MD_ij];
Apos_MD_ij = [ MD_ij(:,11) group_MD 3*ones(length(group_MD),1) subjects_MD_ij];
A_MD_ij = [Apre_MD_ij; Amid_MD_ij; Apos_MD_ij];
save('subject_2_fourcolumn.txt','A_MD_ij', '-ascii')
BWAOV2(A_MD_ij);




%% partner VP mean of md12 and md21 betwwen H and VH
MD_ij = [md_avg]; 
group_MD = [md22_statmat(:,1)];
subjects_MD_ij = (1:length(MD_ij(:,1)))';
%t_MD_ij = table(label(group_MD_ij)',MD_ij(:,2),MD_ij(:,11),'VariableNames',{'group','pre','post'});
time_MD_ij = table([1 2]','VariableNames',{'Time'});

rm_MD_ij = fitrm(t_MD_ij,'post-pre~group','WithinDesign',time_MD_ij);
Apre_MD_ij = [ MD_ij(:,2) group_MD ones(length(group_MD),1) subjects_MD_ij];
Amid_MD_ij = [ MD_ij(:,6) group_MD 2*ones(length(group_MD),1) subjects_MD_ij];
Apos_MD_ij = [ MD_ij(:,11) group_MD 3*ones(length(group_MD),1) subjects_MD_ij];
A_MD_ij = [Apre_MD_ij; Amid_MD_ij; Apos_MD_ij];
save('subject_avg_fourcolumn.txt','A_MD_ij', '-ascii')
BWAOV2(A_MD_ij);

%% partner VP catch trial 
MDc_ij = [md21c_d_t; md12c_d_t]; 
group_MDc_ij = [md22_statmat(:,1); md11_statmat(:,1)];
%gg = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2];
subjects_MDc_ij = (1:length(MDc_ij))';
t_MDc_ij = table(label(group_MDc_ij)',MDc_ij(:,2),MDc_ij(:,11),'VariableNames',{'group','pre','post'});
time_MDc_ij = table([1 2]','VariableNames',{'Time'});

rm_MDc_ij = fitrm(t_MDc_ij,'post-pre~group','WithinDesign',time_MDc_ij);
Apre_MDc_ij = [ MDc_ij(:,3) group_MDc_ij ones(length(group_MDc_ij),1) subjects_MDc_ij];
Amid_MDc_ij = [ MDc_ij(:,6) group_MDc_ij 2*ones(length(group_MDc_ij),1) subjects_MDc_ij];
Apos_MDc_ij = [ MDc_ij(:,10) group_MDc_ij 3*ones(length(group_MDc_ij),1) subjects_MDc_ij];
A_MDc_ij = [Apre_MDc_ij; Amid_MDc_ij; Apos_MDc_ij];
save('catch_fourcolumn.txt','A_MDc_ij', '-ascii')
BWAOV2(A_MDc_ij);

%% Interaction Force per-subject, per-epoch

eff1_d = eff1_statmat(:,2:end);
eff2_d = eff2_statmat(:,2:end);
IF = 0.5*(eff1_d+eff2_d);
len_IF = size(IF);

eff_group = eff1_statmat(:,1);
label = {'C','A'};
t_eff = table(label(eff_group)',IF(:,2),IF(:,6),IF(:,11),'VariableNames',{'group','pre','mid','post'});
time_eff = table([1 2 3]','VariableNames',{'Time'});
subjects_eff = (1:len_IF(1))';

% rm_eff = fitrm(t_eff,'post-mid-pre~group','WithinDesign',time_eff);
Apre_eff = [ IF(:,2) eff_group ones(len_IF(1),1) subjects_eff];
Amid_eff = [ IF(:,6) eff_group 2*ones(len_IF(1),1) subjects_eff];
Apos_eff = [ IF(:,11) eff_group 3*ones(len_IF(1),1) subjects_eff];
A_eff = [Apre_eff; Amid_eff; Apos_eff];
save('intforce_epoch3.txt','A_eff', '-ascii');

BWAOV2(A_eff,0.05);
    
%% Delta crossing time per-subject, per-epoch

tc11d = tc11_statmat(:,2:end);
tc12d = tc12_statmat(:,2:end);
tc21d = tc21_statmat(:,2:end);
tc22d = tc22_statmat(:,2:end);

DTC_1 = tc21d - tc11d; %
DTC_2 = tc12d - tc22d; %
len_DTC = size(DTC_1);

%% DTC_1 
group_DTC_1 = tc11_statmat(:,1);
% label = {'C','A'};
subjects_DTC_1 = (1:len_DTC(1))';
% t_DTC_1 = table(label(group_DTC_1)',DTC_1(:,2),DTC_1(:,11),'VariableNames',{'group','pre','post'});
% time_DTC_1 = table([1 2]','VariableNames',{'Time'});

% rm_DTC_1 = fitrm(t_DTC_1,'post-pre~group','WithinDesign',time_DTC_1);
Apre_DTC_1 = [ DTC_1(:,2) group_DTC_1 ones(len_DTC(1),1) subjects_DTC_1];
Amid_DTC_1 = [ DTC_1(:,6) group_DTC_1 2*ones(len_DTC(1),1) subjects_DTC_1];
Apos_DTC_1 = [ DTC_1(:,11) group_DTC_1 3*ones(len_DTC(1),1) subjects_DTC_1];
A_DTC_1 = [Apre_DTC_1;Amid_DTC_1; Apos_DTC_1];
save('dct1_fourcolumn.txt','A_DTC_1', '-ascii')
BWAOV2(A_DTC_1,0.05);

%% DCT_2
group_DTC_2 = tc22_statmat(:,1);
%label = {'C','A'};
subjects_DTC_2 = (1:len_DTC(1))';
% t_DTC_2 = table(label(group_DTC_2)',DTC_2(:,2),DTC_2(:,11),'VariableNames',{'group','pre','post'});
% time_DTC_2 = table([1 2]','VariableNames',{'Time'});
% 
% rm_DTC_2 = fitrm(t_DTC_2,'post-pre~group','WithinDesign',time_DTC_2);
Apre_DTC_2 = [ DTC_2(:,2) group_DTC_2 ones(len_DTC(1),1) subjects_DTC_2];
Amid_DTC_2 = [ DTC_2(:,6) group_DTC_2 2*ones(len_DTC(1),1) subjects_DTC_2];
Apos_DTC_2 = [ DTC_2(:,11) group_DTC_2 3*ones(len_DTC(1),1) subjects_DTC_2];
A_DTC_2 = [Apre_DTC_2; Amid_DTC_2; Apos_DTC_2];
save('dct2_fourcolumn.txt','A_DTC_2', '-ascii')
BWAOV2(A_DTC_2,0.05);

%%


