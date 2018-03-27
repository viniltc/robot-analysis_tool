function plot_polar(exp,indlist,indlabel)
% Draws polar plots of average behaviors of individual indicators
% indlist: indicators to be displayed
% indlabel: labels of indicators
% (C) V. Sanguineti (2008)

%thstep = 2*pi./exp.Ntargets;
%th=(0:exp.Ntargets-1)*thstep;
th = atan2(exp.task.targets(:,2), exp.task.targets(:,1))';
Ntargets = length(th);

for ind = 1:length(indlist)
    indname = indlist{ind};
    indlab = indlabel{ind};
    
    for s = 1:length(exp.subj)
            [indmat,indmat_catch,inds,indmat_catch_nofor,indmat_catch_norot]=get_indmatrix(exp.subj(s),indname);
            indmat_bydir(s,:) = indmat;  % this is Nreps x Ntargets
            indmat_bydir_catch(s,:) = indmat_catch;
            indmat_ts(s,:) = inds;
            ind_bydir_nofor(s,:) = indmat_catch_nofor; %ADD
            ind_bydir_norot(s,:) = indmat_catch_norot; %ADD
    end

    [Nsu,Nts]=size(indmat_bydir);
    
    % take average over reps
    indavg_null = zeros(Nsu,Ntargets,2);
    indavg_fieldrot = zeros(Nsu,Ntargets,2); %add
    indavg_after = zeros(Nsu,Ntargets,2);
    indavg_rot = zeros(Nsu,Ntargets,2); %%add
    for s = 1:Nsu
        cnt_null = 0;
        cnt_fieldrot = 0;
        cnt_after = 0;
        cnt_rot = 0;
        
        ind_null = [];
        ind_fieldrot = [];
        ind_after = [];
        ind_rot= [];
        
        s_ind_null = [];
        s_ind_fieldrot = [];
        s_ind_after = [];
        s_ind_rot= [];
        
        cv_ind_null = [];
        cv_ind_fieldrot = [];
        cv_ind_after = [];
        cv_ind_rot = [];
        
        
        cont_null = [];
        cont_rot = [];
        cont_rot+forc = [];
        
            
        
        for ts = 1:Nts
          indsub = (indmat_bydir{s,ts});
          inda(ts,:) = nanmean(indsub);
          s_inda(ts,:) = nanstd(indsub);
          cv_inda(ts,:)=s_inda(ts,:)./abs(inda(ts,:));
          
%           switch ts % totally ad-hoc....
%               case {1,2,3,4}
%                   suff = exp.protocol.phases{1};
%               case {5,6,7,8,9,10}
%                   suff = exp.protocol.phases{2};
%               case {11,12}
%                   suff = exp.protocol.phases{3};
%           end        
          % suff=exp.protocol.phasesphaseno{ts};
          a=size(exp.protocol.phases);
          for i=1:a(1,2)
          ts
          
          if(strmatch('1_null',exp.protocol.phases{i}))
             cnt_null = cnt_null+1;
             ind_null(cnt_null,:)=inda(ts,:);
             s_ind_null(cnt_null,:)=s_inda(ts,:);
             cv_ind_null(cnt_null,:)=cv_inda(ts,:);
                        disp('entrato in null');
                                cont_null = cont_null+1;


         elseif (strmatch('3_rot+force',exp.protocol.phases{i}))
             cnt_fieldrot = cnt_fieldrot+1;
             ind_fieldrot(cnt_fieldrot,:)=inda(ts,:);
             s_ind_fieldrot(cnt_fieldrot,:)=s_inda(ts,:);
             cv_ind_fieldrot(cnt_fieldrot,:)=cv_inda(ts,:);
                        disp('entrato in rot+f');
                                cont_rot+forc = cont_rot+forc+1;


    
          elseif  (strmatch('2_rot',exp.protocol.phases{i}))
              cnt_rot=cnt_rot+1;
              ind_rot(cnt_rot,:)=inda(ts,:);
           s_ind_rot(cnt_rot,:)=s_inda(ts,:);
           cv_ind_rot(cnt_rot,:)=cv_inda(ts,:);
           disp('entrato in rot');
           cont_rot = cont_rot+1;
           
           
           elseif  (strmatch('4_rot',exp.protocol.phases{i}))
              cnt_rot=cnt_rot+1;
              ind_rot(cnt_rot,:)=inda(ts,:);
           s_ind_rot(cnt_rot,:)=s_inda(ts,:);
           cv_ind_rot(cnt_rot,:)=cv_inda(ts,:);
           disp('entrato in rot');
           cont_rot = cont_rot+1;
          else
           
             cnt_after = cnt_after+1;
             ind_after(cnt_after,:)=inda(ts,:);
             s_ind_after(cnt_after,:)=s_inda(ts,:);
             cv_ind_after(cnt_after,:)=cv_inda(ts,:);
             
    
          end  
          end
        end
        
        keyboard;
      
        indavg_null(s,:,1) = ind_null(1,:);
        indavg_fieldrot(s,:,1) =  ind_fieldrot(1,:);
        indavg_after(s,:,1) = ind_after(1,:);
        indavg_rot(s,:,1)=ind_rot(1,:);
          
        indavg_null(s,:,2) = ind_null(end,:);
        indavg_fieldrot(s,:,2) =  ind_fieldrot(end,:);
        indavg_after(s,:,2) = ind_after(end,:);
        indavg_rot(s,:,2)= ind_rot(end,:);
        
        
        indstd_null(s,:,1) = s_ind_null(1,:);
        indstd_fieldrot(s,:,1) =  s_ind_fieldrot(1,:);
        indstd_after(s,:,1) = s_ind_after(1,:);
        indstd_rot(s,:,1)=s_ind_rot(1,:);
        
        indstd_null(s,:,2) = s_ind_null(end,:);
        indstd_fieldrot(s,:,2) =  s_ind_fieldrot(end,:);
        indstd_after(s,:,2) = s_ind_after(end,:);
         indstd_rot(s,:,2)=s_ind_rot(end,:);

        
        indcv_null(s,:,1) = cv_ind_null(1,:);
        indcv_fieldrot(s,:,1) =  cv_ind_fieldrot(1,:);
        indcv_after(s,:,1) = cv_ind_after(1,:);
        indcv_rot(s,:,1)=cv_ind_rot(1,:);

        
        indcv_null(s,:,2) = cv_ind_null(end,:);
        indcv_fieldrot(s,:,2) =  cv_ind_fieldrot(end,:);
        indcv_after(s,:,2) = cv_ind_after(end,:);
        indcv_rot(s,:,2)=cv_ind_rot(end,:);

    
    end
    
    rl=get_range(indname);

    figure
    inc_i = mean(indavg_null(find(exp.groups==1),:,1),1);
    inc_f = mean(indavg_null(find(exp.groups==1),:,2),1);
    Ncon = sum(exp.groups==1);
    line((-rl(1)+inc_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',3);
    line((-rl(1)+inc_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',1);
    thr = (0:50)*2*pi/50;
    
    if any(exp.groups==2)
     hold on
     Npat = sum(exp.groups==2);
     inp_i = mean(indavg_null(find(exp.groups==2),:,1),1);
     inp_f = mean(indavg_null(find(exp.groups==2),:,2),1);
     line((-rl(1)+inp_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',3);
     line((-rl(1)+inp_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',1);
     
    end
   
    for i=1:Ntargets
        line((rl(2)-rl(1)).*[0 cos(th(i))], ...
             (rl(2)-rl(1)).*[0 sin(th(i))],'col','k','lines','-');
    end
    line(-rl(1).*cos(thr),-rl(1).*sin(thr),'col','k','linew',2);
    set(gca,'xlim',[rl(1)-rl(2) rl(2)-rl(1)],'ylim',[rl(1)-rl(2) rl(2)-rl(1)],'vis','off');
    set(gcf, 'pos', [100 100 250 250])
    set(gcf,'name',[indlabel{ind}, ': NULL'])
    
    
    figure
    inc_i = mean(indavg_fieldrot(find(exp.groups==1),:,1),1);
    inc_f = mean(indavg_fieldrot(find(exp.groups==1),:,2),1);
    Ncon = sum(exp.groups==1);
    line((-rl(1)+inc_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',3);
    line((-rl(1)+inc_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',1);
    
    if any(exp.groups==2)
     hold on
     Npat = sum(exp.groups==2);
     inp_i = mean(indavg_fieldrot(find(exp.groups==2),:,1),1);
     inp_f = mean(indavg_fieldrot(find(exp.groups==2),:,2),1);
     line((-rl(1)+inp_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',3);
     line((-rl(1)+inp_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',1);
     end
    for i=1:Ntargets
        line((rl(2)-rl(1)).*[0 cos(th(i))], ...
             (rl(2)-rl(1)).*[0 sin(th(i))],'col','k','lines','-');
    end
    line(-rl(1).*cos(thr),-rl(1).*sin(thr),'col','k','linew',2);
    set(gca,'xlim',[rl(1)-rl(2) rl(2)-rl(1)],'ylim',[rl(1)-rl(2) rl(2)-rl(1)],'vis','off');
    set(gcf, 'pos', [400 100 250 250])
    set(gcf,'name',[indlabel{ind}, ': FIELD'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inc_i = mean(indavg_rot(find(exp.groups==1),:,1),1);
    inc_f = mean(indavg_rot(find(exp.groups==1),:,2),1);
    Ncon = sum(exp.groups==1);
    line((-rl(1)+inc_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','g','linew',3);
    line((-rl(1)+inc_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','g','linew',1);
    
    if any(exp.groups==2)
     hold on
     Npat = sum(exp.groups==2);
     inp_i = mean(indavg_rot(find(exp.groups==2),:,1),1);
     inp_f = mean(indavg_rot(find(exp.groups==2),:,2),1);
     line((-rl(1)+inp_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',3);
     line((-rl(1)+inp_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',1);
     end
    for i=1:Ntargets
        line((rl(2)-rl(1)).*[0 cos(th(i))], ...
             (rl(2)-rl(1)).*[0 sin(th(i))],'col','k','lines','-');
    end
    line(-rl(1).*cos(thr),-rl(1).*sin(thr),'col','k','linew',2);
    set(gca,'xlim',[rl(1)-rl(2) rl(2)-rl(1)],'ylim',[rl(1)-rl(2) rl(2)-rl(1)],'vis','off');
    set(gcf, 'pos', [400 100 250 250])
    set(gcf,'name',[indlabel{ind}, ': FIELD'])
    
    
    
    
    
    
    figure
    inc_i = mean(indavg_after(find(exp.groups==1),:,1),1);
    inc_f = mean(indavg_after(find(exp.groups==1),:,2),1);
    Ncon = sum(exp.groups==1);
    
    line((-rl(1)+inc_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',3);
    line((-rl(1)+inc_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (-rl(1)+inc_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',1);
    
    if any(exp.groups==2)
     hold on
     Npat = sum(exp.groups==2);
     inp_i = mean(indavg_after(find(exp.groups==2),:,1),1);
     inp_f = mean(indavg_after(find(exp.groups==2),:,2),1);
     line((-rl(1)+inp_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',3);
     line((-rl(1)+inp_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (-rl(1)+inp_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',1);
        end
    for i=1:Ntargets
        line((rl(2)-rl(1)).*[0 cos(th(i))], ...
             (rl(2)-rl(1)).*[0 sin(th(i))],'col','k','lines','-');
    end
    line(-rl(1).*cos(thr),-rl(1).*sin(thr),'col','k','linew',2);
    set(gca,'xlim',[rl(1)-rl(2) rl(2)-rl(1)],'ylim',[rl(1)-rl(2) rl(2)-rl(1)],'vis','off');
    set(gcf, 'pos', [700 100 250 250])
    set(gcf,'name',[indlabel{ind}, ': ROT'])
    
    % Variability (STD)
    
    figure
    
    inc_i = mean(indstd_null(find(exp.groups==1),:,1),1);
    inc_f = mean(indstd_null(find(exp.groups==1),:,2),1);
    Ncon = sum(exp.groups==1);
    line((inc_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (inc_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',3);
    line((inc_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (inc_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',1);
    
    if any(exp.groups==2)
     hold on
     Npat = sum(exp.groups==2);
     inp_i = mean(indstd_null(find(exp.groups==2),:,1),1);
     inp_f = mean(indstd_null(find(exp.groups==2),:,2),1);
     line((inp_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (inp_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',3);
     line((inp_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (inp_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',1);
     end
    for i=1:Ntargets
        line((rl(2)).*[0 cos(th(i))], ...
             (rl(2)).*[0 sin(th(i))],'col','k','lines','-');
    end
    set(gca,'xlim',2*[-rl(2) rl(2)],'ylim',2*[-rl(2) rl(2)],'vis','off');
    set(gcf, 'pos', [100 100 250 250])
    set(gcf,'name',[indlabel{ind}, ' VARIABILITY: NULL'])
    
    
    
    
    figure
    inc_i = mean(indstd_fieldrot(find(exp.groups==1),:,1),1);
    inc_f = mean(indstd_fieldrot(find(exp.groups==1),:,2),1);
    Ncon = sum(exp.groups==1);
    line((inc_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (inc_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',3);
    line((inc_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (inc_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',1);
    
    if any(exp.groups==2)
     hold on
     Npat = sum(exp.groups==2);
     inp_i = mean(indstd_fieldrot(find(exp.groups==2),:,1),1);
     inp_f = mean(indstd_fieldrot(find(exp.groups==2),:,2),1);
     line((inp_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (inp_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',3);
     line((inp_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (inp_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',1);
     end
    for i=1:Ntargets
        line((rl(2)).*[0 cos(th(i))], ...
             (rl(2)).*[0 sin(th(i))],'col','k','lines','-');
    end
    set(gca,'xlim',2*[-rl(2) rl(2)],'ylim',2*[-rl(2) rl(2)],'vis','off');
    set(gcf, 'pos', [400 100 250 250])
    set(gcf,'name',[indlabel{ind}, ' VARIABILITY: FIELD'])
    
    
    
    figure
    inc_i = mean(indstd_after(find(exp.groups==1),:,1),1);
    inc_f = mean(indstd_after(find(exp.groups==1),:,2),1);
    
    Ncon = sum(exp.groups==1);
    line((inc_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (inc_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',3);
    line((inc_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
         (inc_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','b','linew',1);
    
    if any(exp.groups==2)
     hold on
     Npat = sum(exp.groups==2);
     inp_i = mean(indstd_after(find(exp.groups==2),:,1),1);
     inp_f = mean(indstd_after(find(exp.groups==2),:,2),1);
     line((inp_f([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (inp_f([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',3);
     line((inp_i([1:Ntargets 1])).*cos(th([1:Ntargets 1])),...
          (inp_i([1:Ntargets 1])).*sin(th([1:Ntargets 1])),'color','r','linew',1);
        end
 
    for i=1:Ntargets
        line((rl(2)).*[0 cos(th(i))], ...
             (rl(2)).*[0 sin(th(i))],'col','k','lines','-');
    end
    set(gca,'xlim',[-2*rl(2) 2*rl(2)],'ylim',[-2*rl(2) 2*rl(2)],'vis','off');
    set(gcf, 'pos', [700 100 250 250])
    set(gcf,'name',[indlabel{ind}, ' VARIABILITY: AFTER'])
    
    pause
end

 
%crated new camp field+rot and rot ###under construction ###