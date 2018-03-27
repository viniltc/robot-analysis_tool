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
load([resdir, 'speed1_statmat.txt']);
load([resdir, 'speed1_statmat_catch.txt']);
load([resdir, 'speed2_statmat.txt']);
load([resdir, 'speed2_statmat_catch.txt']);
load([resdir, 'R_statmat.txt']);
load([resdir, 'R_statmat_catch.txt']);

load([resdir, 'eff1_statmat.txt']);
load([resdir, 'eff1_statmat_catch.txt']);
load([resdir, 'eff2_statmat.txt']);
load([resdir, 'eff2_statmat_catch.txt']);

group = speed1_statmat(:,1);

% load per-trial data
load ([resdir,'speed1.mat']);
speed1_ts = sqrt(ind_ts);
load ([resdir,'speed2.mat']);
speed2_ts = sqrt(ind_ts);
load ([resdir,'R.mat']);
R_ts = sqrt(ind_ts);
load([resdir, 'forces.mat']);

speed1_ts(find(forces==0))= nan;
speed2_ts(find(forces==0))= nan;

plotlabel = 'speed1,speed2 [m/s]';
yrange = [0 1]; % range of variation of effort


%% Plot per-subject, per-epoch
speed1 = sqrt(speed1_statmat(:,2:end));
speed2 = sqrt(speed2_statmat(:,2:end));
R = sqrt(R_statmat(:,2:end));
subjs=size(speed1,1);
tsets=size(speed1,2);
allvals1 = [];
allvals2 = [];
for subj=1:subjs
    figure
    set(gcf,'pos',[100 100 300 400])
    patch([1.5 11.5 11.5 1.5],yrange([1 1 2 2]),[0.8 0.8 0.8],'edgecol',[0.8 0.8 0.8]);
    xlabel('epochs')
    ylabel(plotlabel)
    title(sprintf('dyad%d',subj))
       % line(1:tsets,speed1(subj,:),'col','r','marker','s');
       % line(1:tsets,speed2(subj,:),'col','b','marker','s');
        line(1:tsets,R(subj,:),'col','k','marker','s');
       %line(1:tsets,R(subj,:),'col','g','marker','s');
       %line(1:tsets,xcorr(speed1(subj,:),speed2(subj,:)),'col','g','marker','s');
    xlim([1 tsets])
    ylim(yrange)
    
    allvals1 = [allvals1; speed1(subj,:)];
    allvals2 = [allvals2; speed2(subj,:)];

end

 

