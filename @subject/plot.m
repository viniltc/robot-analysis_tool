
% function[lead11 lead12 lead21 lead22]= plot(subj,displist,figdir,epsdir,varargin)
function[lead11 lead12 lead21 lead22 cc1 cc2]= plot(subj,displist,figdir,epsdir,varargin)
% plots SUBJECT data and saves the figures
% displ_list: list of graph types
% figdir: directory where figs are saved
% epsdir: directory where eps figs are saved
% (C) V. Sanguineti (2008)

pow_matrix1 = ones(2000, 156)*nan; % 20 sec (100Hz sampling)
pow_matrix2 = ones(2000, 156)*nan;
tc_matrix1 = ones(2, 156)*nan; 
tc_matrix2 = ones(2, 156)*nan;
tc_mean = zeros(4,13);
dur_trial = zeros(2,156);
perc_trial = zeros(2,156);
perc_vp = zeros(4,156);
mean_vp = zeros(4,156);
allperc_trial = [];
allperc_vp = [];
allmean_vp = [];

timeshift = 30; % 100ms
n_trial = 156;
samp_freq = 100.0;


Red = [.9 0 .4];
Blue = [0 .4 .9];

% filename1 = sprintf('Lead_index1.txt');
% filename2 = sprintf('Lead_index2.txt');
% file1 = fopen(filename1,'wt');
% file2 = fopen(filename2,'wt');

% keyboard

hsvcolormap=zeros(64,3); % 64 x 3 standard colormap dimension
hsvcolormap(:,3)=1; % brillantezza max
hsvcolormap(1:31,1)=0.8333; % magenta
hsvcolormap(34:end,1)=0.5; % ciano
for i=1:(32-1)
    hsvcolormap(i+1,2) = 1-(i*0.032); % 1/31
    hsvcolormap(i+32+1,2) = (i*0.032);
end
hsvcolormap(1,2) = 1;
hsvcolormap(end,2) = 1;
hsvcolormap(32:33,1:2)=0; % bianco
newcolormap = hsv2rgb(hsvcolormap);
% newcolormap(1,1:3) = 0;

% keyboard
tsno = 1;
trial = 12;
for p = 1:length(subj.phases)
  for r = 1:length(subj.order{p})
      tset = subj.tset{tsno}; 
      for i = 1:length(displist)
             switch(displist{i})
               case 'traj'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_trajs(tset,tsno,varargin{1});
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
                   
               case 'trajforce'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_onetraj(tset,tsno,varargin{1});
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
                   
               case 'power'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_powers(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'min_viadist'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_vpdistances(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'int_force'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_intforce(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'sumpower'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_sumpowers(tset,tsno,lim);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
           
               case 'peakacc'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_peakaccs(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'speed'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_speeds(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'jerk'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_jerks(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case '1traj'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_1trajs(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'vpdist'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_vpdistances(tset,tsno); % this is for dyads...
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'variability'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_variability(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'speed_corr'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_correlations(tset,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])
               case 'scores'
                   figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
                   spcs = find(figname==' ');
                   figname(spcs)=[];
                   plot_scores(tset1,tsno);
                   saveas(gcf,[figdir,figname],'fig');
                   preprint;
                   eval(['print -depsc ', epsdir,figname])                   
               case 'leadership'
                    [pow_matrix1,pow_matrix2,tc_matrix1,tc_matrix2, perc_trial, perc_vp, mean_vp]=plot_leaderships(tset,tsno,pow_matrix1,pow_matrix2,tc_matrix1,tc_matrix2, perc_trial, perc_vp, mean_vp);
%                  keyboard
                    dur_trial(1,:) = sum(~(isnan(pow_matrix1())));
                    dur_trial(2,:) = sum(~(isnan(pow_matrix2())));
                    pow_matrix1(find(isnan(pow_matrix1)))=0;
                    pow_matrix2(find(isnan(pow_matrix2)))=0;
                    
%                     tc1_shift(tc1_shift==0)= mean(tc1_shift(1,:),2);
%                     pow_matrix1(pow_matrix1==0)=mean;
%                     pow_matrix2(find(isnan(pow_matrix2)))=0;
% for i=1:2000                    
% pow_matrix1(pow_matrix1(i,:)==0)=mean(pow_matrix1(i,:),2)    
% end



% for i=1:2000
% 
%  if pow_matrix1(i,:)==0 
%    pow_matrix1(i,:)=[];
%  end 
% end 

                    pow_matrix1(pow_matrix1==0)=(mean(mean(pow_matrix1,1),2));
                    pow_matrix2(pow_matrix2==0)=(mean(mean(pow_matrix2,1),2));
                    allperc_trial = [allperc_trial; perc_trial]; %add 2 rows at a time
                    allperc_vp = [allperc_vp; perc_vp]; %add 4 rows at a time
                    allmean_vp = [allmean_vp; mean_vp]; %add 4 rows at a time
                    for t=1:13
                        
                        tc_mean(1,t) = nanmean(tc_matrix1(1,1+(12*(t-1)):12+(12*(t-1))));
                        tc_mean(2,t) = nanmean(tc_matrix1(2,1+(12*(t-1)):12+(12*(t-1))));
                        tc_mean(3,t) = nanmean(tc_matrix2(1,1+(12*(t-1)):12+(12*(t-1))));
                        tc_mean(4,t) = nanmean(tc_matrix2(2,1+(12*(t-1)):12+(12*(t-1))));
                        dur_mean(1,t) = nanmean(dur_trial(1,1+(12*(t-1)):12+(12*(t-1))));
                        dur_mean(2,t) = nanmean(dur_trial(2,1+(12*(t-1)):12+(12*(t-1))));
                    end
                    t_trial = [12:12:121];
                    t_trial = t_trial-6;
                    if (tsno == 12)
                        figname = sprintf('%s%s%d%s',subj.name,displist{i},1)     
                        spcs = find(figname==' ');
                        figname(spcs)=[];
                        figure
 
                        set(gcf,'pos',[100 100 300 150])
                        colormap(newcolormap);        
%                         image((1:12:120),[0 4],pow_matrix1(1:400,12:132), 'CDataMapping', 'scaled') % guardiamo solo i primi 6 sec.
                        pow_matrix1(find(isnan(pow_matrix1)))=0;
                        pow_matrix1(pow_matrix1==0) = [];
                        pow_matrix1(pow_matrix1==0)=(mean(mean(pow_matrix1,1),2));
                        
%                         keyboard
%                         image([1 120],[0 4],pow_matrix1(1:400,12:132), 'CDataMapping', 'scaled')
%                         image([1 120],[0 4],pow_matrix1(1:400,12:132), 'CDataMapping', 'scaled')
                        image([0 2.5],[1 12],pow_matrix1(1:250,120:132)', 'CDataMapping', 'scaled')
                        image([0 2.5],[1 120],pow_matrix1(1:250,12:132)', 'CDataMapping', 'scaled')

                        
                        %  1 120......12 132
                        axis xy
%                       colorbar
                        ax = gca;
                        ax.CLim=[-0.2 0.2];
                        skipLabel = 2;
                        ax.YTick = ax.YTick(1:skipLabel:end)
                        hold on
%                         plot(t_trial,tc_mean(1,2:11),'b-', 'LineWidth',1);
%                         plot(t_trial,tc_mean(3,2:11),'r--', 'LineWidth',1);
    
%  keyboard
                        plot(tc_mean(1,2:11),t_trial,'--','Color',Blue,'LineWidth', 1);
                        plot(tc_mean(4,2:11),t_trial,'--','Color',Red,'LineWidth', 1);
%                           plot(tc_matrix1(1,120:132),(1:13),'--','Color',Blue, 'LineWidth',2);
%                           plot(tc_matrix2(2,120:132),(1:13),'--','Color',Red,'LineWidth',2);
%                           
% %                           e1=plot(tc_matrix1(1,12:132),(1:121),'--','Color',Blue, 'LineWidth',1);
% %                           e2=plot(tc_matrix2(2,12:132),(1:121),'--','Color',Red,'LineWidth',1);
%                           plot(double(int32(tc_matrix1(1,120:132)*samp_freq)-timeshift)/samp_freq,(1:13),'--','Color',Blue, 'LineWidth',1);
%                           plot(double(int32(tc_matrix2(2,120:132)*samp_freq)-timeshift)/samp_freq,(1:13),'--','Color',Red, 'LineWidth',1);
%                           set(gca, 'YTick', []);
%                         hold on
%                         e3=plot(t_trial,dur_mean(1,2:11)./100,'k--');
                        %       l1=legend([e2,e1,e3],'own vp','other vp','avg duration');
                        %        set(l1,'Location','northeast','FontName', 'Times New Roman')
                        xlabel('Time [s]','FontName', 'Times New Roman','fontsize',14)
                        ylabel('Training','FontName', 'Times New Roman','fontsize',14)
                        title('S_{1}','FontName', 'Times New Roman','Color','b','fontsize',14, 'FontWeight','Normal')
                        saveas(gcf,[figdir,figname],'fig');
                        preprint;
                        eval(['print -depsc ', epsdir,figname])

                        figname = sprintf('%s%s%d%s',subj.name,displist{i},2)
                        spcs = find(figname==' ');
%                         spdstrval=['Leadership - ', num2str(phaseno)];
                        
                        figname(spcs)=[];
                        figure
                        
%                         set(gcf,'pos',[100 100 370 150])
                        set(gcf,'pos',[100 100 370 150])
                        colormap(newcolormap);
                        pow_matrix2(find(isnan(pow_matrix2)))=0;
                        pow_matrix2(pow_matrix2==0) = [];
                        pow_matrix2(pow_matrix2==0)=(mean(mean(pow_matrix2,1),2));
                      %  pow_matrix2(pow_matrix2==0) = pow_matrix2(mean(pow_matrix2));
%                         image((1:12:120),[0 4],pow_matrix2(1:400,12:132), 'CDataMapping', 'scaled')
%                         image([1 120],[0 4],pow_matrix2(1:400,12:132), 'CDataMapping', 'scaled')
%                         image([1 120],[0 4],pow_matrix2(1:400,12:132), 'CDataMapping', 'scaled')
                        image([0 2.5],[1 12],pow_matrix2(1:250,120:132)', 'CDataMapping', 'scaled')
                        image([0 2.5],[1 120],pow_matrix2(1:250,12:132)', 'CDataMapping', 'scaled')
                        axis xy
                      colorbar
                        ax = gca;
                        ax.CLim=[-0.2 0.2];
                        skipLabel = 2;

                        ax.YTick = ax.YTick(1:skipLabel:end)
                        hold on
%                         e1=plot(t_trial,tc_mean(2,2:11),'b--', 'LineWidth',1);
%                         e2=plot(t_trial,tc_mean(4,2:11),'r-', 'LineWidth',1);
                         plot(tc_mean(1,2:11),t_trial,'--','Color',Blue,'LineWidth', 1);
                         plot(tc_mean(4,2:11),t_trial,'--','Color',Red,'LineWidth', 1);
%                           plot(tc_matrix1(1,120:132),(1:13),'--','Color',Blue, 'LineWidth',2);
%                           plot(tc_matrix2(2,120:132),(1:13),'--','Color',Red,'LineWidth',2);
%                           
% %                           e1=plot(tc_matrix1(1,12:132),(1:121),'--','Color',Blue, 'LineWidth',1);
% %                           e2=plot(tc_matrix2(2,12:132),(1:121),'--','Color',Red,'LineWidth',1);
%                           plot(double(int32(tc_matrix1(1,120:132)*samp_freq)-timeshift)/samp_freq,(1:13),'--','Color',Blue, 'LineWidth',1);
%                           plot(double(int32(tc_matrix2(2,120:132)*samp_freq)-timeshift)/samp_freq,(1:13),'--','Color',Red, 'LineWidth',1);

                        set(gca, 'YTick', []);
%                         hold on
%                         e3=plot(t_trial,dur_mean(2,2:11)./100,'k--');
                        %       l1=legend([e2,e1,e3],'own vp','other vp','avg duration');
                        %       set(l1,'Location','northeast','FontName', 'Times New Roman')
                        xlabel('Time [s]','FontName', 'Times New Roman','fontsize',14)
%                       ylabel('Late training','FontName', 'Times New Roman','fontsize',14)
                        title('S_{2}','FontName', 'Times New Roman','Color','r','fontsize',14, 'FontWeight','Normal')
                        saveas(gcf,[figdir,figname],'fig');
                        preprint;
                        eval(['print -depsc ', epsdir,figname])
                        
%                         figname = sprintf('%s%s%d%s',subj.name,displist{i},1,'LI1')     
%                         spcs = find(figname==' ');
%                         figname(spcs)=[];
%                         figure
%                         set(gcf,'pos',[100 100 200 200])                             
%                         plot(mean_vp(1,12:132), 'b-')
%                         hold on
%                         plot(mean_vp(2,12:132), 'r--')
%                         axis xy
%                         hold on
%                         legend('S1','S2')
%                         xlabel('trials','FontName', 'Times New Roman')
%                         ylabel('power [W]','FontName', 'Times New Roman')
%                         title('LI1 - interval on VP1 tc','FontName', 'Times New Roman','Color','b')
% %                         saveas(gcf,[figdir,figname],'fig');
% %                         preprint;
% %                         eval(['print -depsc ', epsdir,figname])
% 
%                         figname = sprintf('%s%s%d%s',subj.name,displist{i},1,'LI2')     
%                         spcs = find(figname==' ');
%                         figname(spcs)=[];
%                         figure
%                         set(gcf,'pos',[100 100 200 200])                             
%                         plot(mean_vp(4,12:132), 'r-')
%                         hold on
%                         plot(mean_vp(3,12:132), 'b--')
%                         axis xy
%                         hold on
%                         legend('S2','S1')
%                         xlabel('trials','FontName', 'Times New Roman')
%                         ylabel('power [W]','FontName', 'Times New Roman')
%                         title('LI2 - interval on VP2 tc','FontName', 'Times New Roman','Color','r')
% %                         saveas(gcf,[figdir,figname],'fig');
% %                         preprint;
% %                         eval(['print -depsc ', epsdir,figname])
% %                   keyboard
               

                   
                   
       
%                         image((1:12:120),[0 4],pow_matrix1(1:400,12:132), 'CDataMapping', 'scaled') % guardiamo solo i primi 6 sec.
                       
                       
                        tc1 = int32(samp_freq*tc_matrix1);
                        tc2 = int32(samp_freq*tc_matrix2);
%                         tc1(tc1==0)= mean(tc1(1,:),2);
%                         tc2(tc2==0)= mean(tc2(1,:),2);
                        tc1(tc1==0)= nan;
                        tc2(tc2==0)= nan;
                        
                        tc1_shift = tc1-timeshift; %% 100ms time before tc1
                        tc2_shift = tc2-timeshift; %% 100ms time before tc2
                        
                        tc1_shift(tc1_shift==0)= nanmean(tc1_shift(1,:),2); 
                        tc2_shift(tc2_shift==0)= nanmean(tc2_shift(1,:),2);
                        tc1_shift(tc1_shift<0)= nanmean(tc1_shift(1,:),2);
                        tc2_shift(tc2_shift<0)= nanmean(tc2_shift(1,:),2);
                        
%                         tc1_shift(tc1_shift==0)= nan; 
%                         tc2_shift(tc2_shift==0)= nan;
%                         tc1_shift(tc1_shift<0)= nan;
%                         tc2_shift(tc2_shift<0)= nan;
                        
                        pow1_tc1 = zeros(1,n_trial);
                        pow2_tc2 = zeros(1,n_trial);
                        pow1_tc2 = zeros(1,n_trial);
                        pow2_tc1 = zeros(1,n_trial);
                        
                        for i= 1:n_trial
%                         pow1_tc1(1,i)= nanmean(pow_matrix1(tc1_shift(1, i):tc1(1,i),i),1);
% %                       pow2_tc2(1,i)= pow_matrix2(tc2(1, i),i);
%                         pow2_tc2(1,i)= nanmean(pow_matrix2(tc2_shift(1, i):tc2(1,i),i),1);
%                         pow1_tc2(1,i)= nanmean(pow_matrix1(tc2_shift(1, i):tc2(1,i),i),1);
%                         pow2_tc1(1,i)= nanmean(pow_matrix2(tc1_shift(1, i):tc1(1,i),i),1);
                        
                        pow1_tc1(1,i)= nanmean(pow_matrix1(tc1_shift(1, i):tc1(1,i),i),1);
%                       pow2_tc2(1,i)= pow_matrix2(tc2(1, i),i);
                        pow2_tc2(1,i)= nanmean(pow_matrix2(tc2_shift(2, i):tc2(2,i),i),1);
                        pow1_tc2(1,i)= nanmean(pow_matrix1(tc2_shift(1, i):tc2(1,i),i),1);
                        pow2_tc1(1,i)= nanmean(pow_matrix2(tc1_shift(2, i):tc1(2,i),i),1);
                        
%                         pow1_tc1(1,i)= (pow_matrix1(tc1_shift(1, i),i));
% %                       pow2_tc2(1,i)= pow_matrix2(tc2(1, i),i);
%                         pow2_tc2(1,i)= (pow_matrix2(tc2_shift(1, i),i));
%                         pow1_tc2(1,i)= (pow_matrix1(tc2_shift(1, i),i));
%                         pow2_tc1(1,i)= (pow_matrix2(tc1_shift(1, i),i));
                        end
                        
                       early_epoch = 12:24;
                       mid_epoch   = 60:72;
                       late_epoch  = 120:132;
                        
                       pre11 = nanmean(pow1_tc1(:,early_epoch)); 
                       mid11 = nanmean(pow1_tc1(:, mid_epoch)); 
                       post11 = nanmean(pow1_tc1(:,late_epoch)); 
                       
                       pre12 = nanmean(pow1_tc2(:,early_epoch)); 
                       mid12 = nanmean(pow1_tc2(:,mid_epoch)); 
                       post12 = nanmean(pow1_tc2(:,late_epoch));
                      
                       pre21 = nanmean(pow2_tc1(:,early_epoch)); 
                       mid21 = nanmean(pow2_tc1(:,mid_epoch)); 
                       post21 = nanmean(pow2_tc1(:,late_epoch)); 
                      
                       pre22 = nanmean(pow2_tc2(:,early_epoch)); 
                       mid22 = nanmean(pow2_tc2(:,mid_epoch)); 
                       post22 = nanmean(pow2_tc2(:,late_epoch));
                       
                       cc1 = pow_matrix1;
                       cc2 = pow_matrix2;
                         
                    lead11 = [pre11 mid11 post11];
                    lead12 = [pre12 mid12 post12];
                    lead21 = [pre21 mid21 post21];
                    lead22 = [pre22 mid22 post22];
%          keyboard          

%             
%                     
%                     pre11  = mean(mean_vp(1, 12:24));
%                     mid11  = mean(mean_vp(1, 60:72));
%                     post11 = mean(mean_vp(1, 120:132));
%                     pre12  = mean(mean_vp(2, 12:24));
%                     mid12  = mean(mean_vp(2, 60:72));
%                     post12 = mean(mean_vp(2, 120:132));
%                     pre21  = mean(mean_vp(3, 12:24));
%                     mid21  = mean(mean_vp(3, 60:72));
%                     post21 = mean(mean_vp(3, 120:132));
%                     pre22  = mean(mean_vp(4, 12:24));
%                     mid22  = mean(mean_vp(4, 60:72));
%                     post22 = mean(mean_vp(4, 120:132));

                  
 
                    end
                    
             end          
      end   % of i
    tsno = tsno+1;
   end % of r
 end % of p

end


%% test

% function plot(subj,displist,figdir,epsdir,varargin)
% % plots SUBJECT data and saves the figures
% % displ_list: list of graph types
% % figdir: directory where figs are saved
% % epsdir: directory where eps figs are saved
% % (C) V. Sanguineti (2008)
% 
% pow_matrix1 = ones(2000, 156)*nan; % matrice per tutti i trial fino a 20 secondi
% pow_matrix2 = ones(2000, 156)*nan;
% tc_matrix1 = ones(2, 156)*nan; % matrice per tutti i trial fino a 20 secondi
% tc_matrix2 = ones(2, 156)*nan;
% 
% hsvcolormap=zeros(64,3); % 64 x 3 standard colormap dimension
% hsvcolormap(:,3)=1; % brillantezza max
% hsvcolormap(1:31,1)=0.8333; % magenta
% 
% hsvcolormap(34:end,1)=0.5; % ciano
% for i=1:(32-1)
%     hsvcolormap(i+1,2) = 1-(i*0.032); % 1/31
%     hsvcolormap(i+32+1,2) = (i*0.032);
% end
% hsvcolormap(1,2) = 1;
% hsvcolormap(end,2) = 1;
% hsvcolormap(32:33,1:2)=0; % bianco
% newcolormap = hsv2rgb(hsvcolormap);
% 
% 
% % newcolormap(1,1:3) = 0;
% % keyboard
% tsno = 1;
% trial = 12;
% for p = 1:length(subj.phases)
%   for r = 1:length(subj.order{p})
%    tset = subj.tset{tsno}; 
%    for i = 1:length(displist)
%         if (strcmp(displist{i},'leadership'))
% %           figname = sprintf('%s%s%d%s',subj.name,subj.phases{p},subj.order{p}(r),displist{i})
% % 
% %             spcs = find(figname==' ');
% %             figname(spcs)=[];
% %             switch(displist{i})
% %                case 'traj',  plot_trajs(tset,tsno,varargin{1});
% %                case 'power',plot_powers(tset,tsno,lim);
% %                case 'minviadist',plot_vpdistances(tset,tsno,lim);
% %                case 'intforce',plot_intforce(tset,tsno,lim);           
% %                case 'sumpower',plot_sumpowers(tset,tsno,lim);
% %                case 'trackerror',plot_tracking_errors(tset,tsno,lim,varargin{1});
% % 
% %                case 'peakacc', plot_peakaccs(tset,tsno);
% %                case 'speed',plot_speeds(tset,tsno);
% %                case 'jerk',plot_jerks(tset,tsno);    
% %                case '1traj', plot_1trajs(tset,tsno);
% %                case 'vpdist',plot_vpdistances(tset,tsno); % this is for dyads...
% %                case 'variability',plot_variability(tset,tsno);
% %                case 'speed_corr',plot_correlations(tset,tsno);
% %                case 'scores',plot_scores(tset1,tsno);
% %             end
% %            title(subj.name);  
% %            saveas(gcf,[figdir,figname],'fig');
% %            preprint;
% %            eval(['print -depsc ', epsdir,figname])
% %        else       
%         switch(displist{i})
%             case 'leadership', [pow_matrix1, pow_matrix2,tc_matrix1,tc_matrix2]=plot_leaderships(tset,tsno,pow_matrix1,pow_matrix2,tc_matrix1,tc_matrix2);
%         end
%        end % of if
%     end   % of i
%     tsno = tsno+1;
%    end % of r
%  end % of p
%  
%  
% 
%  tc_mean = zeros(4,13);
%  dur_trial = zeros(2,156);
% for i = 1:length(displist)
%     if (strcmp(displist{i},'leadership'))
% %         perc1 = sum(sum(pow_matrix1>0));
% %         perctot1 = sum(sum(~(isnan(pow_matrix1))));
% %         perc2 = sum(sum(pow_matrix2>0));
% %         perctot2 = sum(sum(~(isnan(pow_matrix2))));
%       dur_trial(1,:) = sum(~(isnan(pow_matrix1())));
%       dur_trial(2,:) = sum(~(isnan(pow_matrix2())));
%       pow_matrix1(find(isnan(pow_matrix1)))=0;
%       pow_matrix2(find(isnan(pow_matrix2)))=0;
%       for t=1:13
%           tc_mean(1,t) = nanmean(tc_matrix1(1,1+(12*(t-1)):12+(12*(t-1))));
%           tc_mean(2,t) = nanmean(tc_matrix1(2,1+(12*(t-1)):12+(12*(t-1))));
%           tc_mean(3,t) = nanmean(tc_matrix2(1,1+(12*(t-1)):12+(12*(t-1))));
%           tc_mean(4,t) = nanmean(tc_matrix2(2,1+(12*(t-1)):12+(12*(t-1))));
%           dur_mean(1,t) = nanmean(dur_trial(1,1+(12*(t-1)):12+(12*(t-1))));
%           dur_mean(2,t) = nanmean(dur_trial(2,1+(12*(t-1)):12+(12*(t-1))));
%       end
%        t_trial = [12:12:121];
%        t_trial = t_trial-6;
%        figname = sprintf('%s%s%d%s',subj.name,displist{i},1)     
%         spcs = find(figname==' ');
%         figname(spcs)=[];
%         figure
%         set(gcf,'pos',[100 100 200 200])
%         colormap(newcolormap);        
% 
%         image((1:12:120),[0 4],pow_matrix1(1:400,12:132), 'CDataMapping', 'scaled') % guardiamo solo i primi 3 sec.
%         axis xy
%         colorbar
%         ax = gca;
%         ax.CLim=[-0.6 0.6];
%         hold on
%         e1=plot(t_trial,tc_mean(1,2:11),'b-');
%         hold on
%         e2=plot(t_trial,tc_mean(2,2:11),'r-');
%         hold on
%         e3=plot(t_trial,dur_mean(1,2:11)./100,'k--');
%  %       l1=legend([e2,e1,e3],'own vp','other vp','avg duration');
% %        set(l1,'Location','northeast','FontName', 'Times New Roman')
%         xlabel('epochs','FontName', 'Times New Roman')
%         ylabel('time [s]','FontName', 'Times New Roman')
%         title('S1','FontName', 'Times New Roman','Color','b')
%         saveas(gcf,[figdir,figname],'fig');
%         preprint;
%         eval(['print -depsc ', epsdir,figname])
%        
%         figname = sprintf('%s%s%d%s',subj.name,displist{i},2)
%         spcs = find(figname==' ');
%         figname(spcs)=[];
%         figure
%         set(gcf,'pos',[100 100 200 200])
%         colormap(newcolormap);
%         pow_matrix2(find(isnan(pow_matrix2)))=0;
%         image((1:12:120),[0 4],pow_matrix2(1:400,12:132), 'CDataMapping', 'scaled')
%         axis xy
%         colorbar
%         ax = gca;
%         ax.CLim=[-0.6 0.6];
%         hold on
%         e1=plot(t_trial,tc_mean(3,2:11),'b-');
%         hold on
%         e2=plot(t_trial,tc_mean(4,2:11),'r-');
%         hold on
%         e3=plot(t_trial,dur_mean(2,2:11)./100,'k--');
%  %       l1=legend([e2,e1,e3],'own vp','other vp','avg duration');
%  %       set(l1,'Location','northeast','FontName', 'Times New Roman')
%         xlabel('epochs','FontName', 'Times New Roman')
%         ylabel('time [s]','FontName', 'Times New Roman')
%         title('S2','FontName', 'Times New Roman','Color','r')
%         saveas(gcf,[figdir,figname],'fig');
%         preprint;
%         eval(['print -depsc ', epsdir,figname])
%     end
% end
%            
% 
% 
% end





