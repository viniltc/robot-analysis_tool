function plot_indicator(exp,indlist,indlabel,indcatch,epsdir,figdir)
% plots average behaviors of individual indicators and saves the figures
% indlist: indicators to be displayed
% indlabel: labels of indicators
% displ_list: list of graph types
% figdir: directory where figs are saved
% epsdir: directory where eps figs are saved
% (C) V. Sanguineti (2008)


for ind = 1:length(indlist)
    indname = indlist{ind};
    indlab = indlabel{ind};
    
    for s = 1:length(exp.subj)
            [indmat,indmat_catch,inds,indmat_catch_nofor, indmat_catch_norot]=get_indmatrix(exp.subj(s),indname);
            indmat_bydir(s,:) = indmat;  % this is Nreps x Ntargets
            indmat_bydir_catch(s,:) = indmat_catch;
            indmat_ts(s,:) = inds;
            ind_bydir_nofor(s,:) = indmat_catch_nofor; %add
            ind_bydir_norot(s,:) = indmat_catch_norot; %add
    end
    
    [Nsu,Nts]=size(indmat_bydir);
    
    lastnull = length(exp.protocol.order);
    
    % take average over dirs 
    for s = 1:Nsu % loop over subjects...
        Ncnt = 0;
        Ngrp = 1;  % 1 points per each tset
        
        for ts = 1:Nts % loop over target sets...
          
          ts
          
          % take out outliers
          o_med = nanmedian(indmat_bydir{s,ts});
          o_std = nanstd(indmat_bydir{s,ts});
          [is_outl1,is_outl2] = find(abs(indmat_bydir{s,ts} - o_med(ones(size(indmat_bydir{s,ts},1),1),:))> ...
                                     3*(o_std(ones(size(indmat_bydir{s,ts},1),1),:)));
          indmat_bydir{s,ts}(is_outl1,is_outl2)=NaN*indmat_bydir{s,ts}(is_outl1,is_outl2);
         
          
          indsub = abs(indmat_bydir{s,ts}); % this is used for adaptation plots: size is nrep X ndirs
         
          size(indsub)
          Nreps = size(indsub,1);
          Navg = round(Nreps/Ngrp);
          
          % this is a binned version of indsub: size is now 1 x Ngrp
          indavg(s,Ncnt+(1:Ngrp)) = nanmean(reshape((nanmean(indsub')),Navg,Ngrp)); 
          %keyboard;
          if indname == indcatch  
              keyboard;
             % this is error at last phases (10:12) of null field: size is 1 x ndirs  
             meanlast = nanmedian((indmat_bydir{s,lastnull}(end+(-1:0),:)));  
            
             % this is like indsub, but with sign...: size is nrep x ndirs
             indsub_noa = (indmat_bydir{s,ts})-meanlast(ones(Nreps,1),:);
             indsub_noa_catch = indmat_bydir_catch{s,ts}(:,1)-meanlast(ones(size(indmat_bydir_catch{s,ts},1),1),1);
            
             % Now take indsub_noa and bin it (Ngrp is bin size)
             % Size is now 1 x Ngrp...
             indavg_noa(s,Ncnt+(1:Ngrp)) = nanmean(reshape(nanmean(indsub_noa'),Navg,size(indsub_noa,1)/Navg));
             indavg_noa_catch(s,Ncnt+(1:Ngrp)) = nanmean(indsub_noa_catch');
%             if isnan(nanmean(indsub_noa_catch'))
%                 indavg_noa_catch(s,Ncnt+(1:Ngrp)) = indavg_noa(s,Ncnt+(1:Ngrp));
%             end
%             
            
            li_noa(s,Ncnt+(1:Ngrp))=abs(indavg_noa_catch(s,Ncnt+(1:Ngrp)))./...
                                    abs(indavg_noa(s,Ncnt+(1:Ngrp))-indavg_noa_catch(s,Ncnt+(1:Ngrp)));

          end
          Ncnt = Ncnt+Ngrp;
          
                   
        end
    end
    
    inc = nanmean(indavg(find(exp.groups==1),:));
    s_inc = nanstd(indavg(find(exp.groups==1),:));
    
    if indname == indcatch
     inc_noa = nanmean(indavg_noa(find(exp.groups==1),:));
     s_inc_noa = nanstd(indavg_noa(find(exp.groups==1),:));
     
     inc_noa_catch = nanmean(indavg_noa_catch(find(exp.groups==1),:));
     s_inc_noa_catch = nanstd(indavg_noa_catch(find(exp.groups==1),:));
     
     inc_li =nanmean(li_noa(find(exp.groups==1),:));
     s_inc_li =nanstd(li_noa(find(exp.groups==1),:));    
    end
    Ncon = sum(exp.groups==1);
    
    % Figure with absolute value of indicator
    figure
    h(1)=line(1:length(inc),inc,'col','b','LineWidth',2);%,'lines','-' );
    line(1:length(inc),inc-s_inc./sqrt(Ncon),'col','b','LineWidth',1);%,'lines','-')
    line(1:length(inc),inc+s_inc./sqrt(Ncon),'col','b','LineWidth',1);%,'lines','-')
    hold on
    if any(exp.groups==2)
     Npat = sum(exp.groups==2);
     inp = nanmean(indavg(find(exp.groups==2),:));
     s_inp = nanstd(indavg(find(exp.groups==2),:));
    
     if indname == indcatch
    
         inp_noa = nanmean(indavg_noa(find(exp.groups==2),:));
         s_inp_noa = nanstd(indavg_noa(find(exp.groups==2),:));
        
         inp_noa_catch = nanmean(indavg_noa_catch(find(exp.groups==2),:));
         s_inp_noa_catch = nanstd(indavg_noa_catch(find(exp.groups==2),:));
         
         inp_li =nanmean(li_noa(find(exp.groups==2),:));
         s_inp_li =nanstd(li_noa(find(exp.groups==2),:));    
    
     end
     
     h(2)=line(1:length(inp),inp,'col','r','LineWidth',2);%,'lines','--' );
     line(1:length(inp),inp-s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
     line(1:length(inp),inp+s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
     legend(h,upper(exp.protocol.groupdesc{1}),upper(exp.protocol.groupdesc{2}));
     legend boxoff
    end
    xlabel('MOVEMENT SET')
    ylabel(indlabel{ind})
    set(gca, 'box', 'off')
    %axis square
    set(gcf, 'pos', [100 100 450 300])
    set(gcf,'name',indlabel{ind})
    
    
    saveas(gcf,[figdir,upper(indname)],'fig');

    %preprint;
    eval(['print -depsc ',epsdir,upper(indname)])
    
    
    % Figure with signed value of indicator (and catch trial)
    if indname == indcatch
      figure
      plot(inc_noa,'b','LineWidth',2 )
      hold on, plot(inc_noa_catch,'b','LineWidth',2, 'lines',':' )
      
      if any(exp.groups==2)
         hold on, plot(inp_noa,'r','LineWidth',2 )
         hold on, plot(inp_noa_catch,'r','LineWidth',2, 'lines',':' )
      end
      xlabel('MOVEMENT SET')
      ylabel(indlabel{ind})
      set(gca, 'box', 'off')
      %axis square
      set(gcf, 'pos', [100 100 350 300])
      set(gcf,'name',indlabel{ind})
    
      saveas(gcf,[figdir,upper(indname),'_withcatch'],'fig');
      preprint;
      eval(['print -depsc ',epsdir,upper(indname),'_withcatch'])
    
      
      % Figure with learning index
      figure % learning index plot
      set(gcf, 'pos', [100 100 350 300])
      isokc = find(~isnan(inc_li));

      nli = length(isokc)+1;
      
      plot([0 inc_li(isokc)],'b','LineWidth',2 )
      hold on, plot([0 inc_li(isokc)-s_inc_li(isokc)./sqrt(Ncon)],'b','LineWidth',1 )
      hold on, plot([0 inc_li(isokc)+s_inc_li(isokc)./sqrt(Ncon)],'b','LineWidth',1 )
     
      if any(exp.groups==2)

         isokp = find(~isnan(inp_li));
         hold on, plot([0 inp_li(isokp)],'r','LineWidth',2 )
         hold on, plot([0 inp_li(isokp)-s_inp_li(isokp)./sqrt(Npat)],'r','LineWidth',1 )
         hold on, plot([0 inp_li(isokp)+s_inp_li(isokp)./sqrt(Npat)],'r','LineWidth',1 )  
      end
      xlabel('MOVEMENT SET')
      ylabel('LEARNING INDEX')
      set(gca, 'box', 'off','ylim',[0 1])
      
      
      %axis square
      set(gcf, 'pos', [100 100 350 300])
      set(gca,'xlim',[1 nli],'xtick',2:nli,'xticklabel','F')
      set(gcf,'name',indlabel{ind})
    
      saveas(gcf,[figdir,upper(indname), '_learningindex'],'fig');
      preprint;
      eval(['print -depsc ',epsdir,upper(indname),'_learningindex'])
    
      
      
      li_bysub_c = nanmean(li_noa(find(exp.groups==1),isokc)');
      s_li_bysub_c = nanstd(li_noa(find(exp.groups==1),isokc)')./sqrt(length(isokc));
  
      if any(exp.groups==2)
          li_bysub_p = nanmean(li_noa(find(exp.groups==2),isokp)');
          s_li_bysub_p = nanstd(li_noa(find(exp.groups==2),isokp)')./sqrt(length(isokp));
      
          
          figure % learning index histogram
          set(gcf, 'pos', [100 100 350 300]) 
          
          hl=errorbar([1:Ncon ],[li_bysub_c],[s_li_bysub_c],'k.');
          set(hl, 'markers',20);
          hold on
          
          hl=errorbar([(Ncon+2+(1:Npat))],[li_bysub_p],[s_li_bysub_p],'k.');
          set(hl, 'markers',20);
          keyboard
          hlb = bar([li_bysub_c; 0; 0; li_bysub_p]);
          hp = findobj(hlb, 'type','patch');
          set(hp,'facecol',[0.8 0.8 0.8]);
          set(gca,'xtick',[Ncon/2 Ncon+2+Npat/2],'xticklabel',{upper(exp.disease_desc{1}),upper(exp.disease_desc{2})},...
                  'xlim',[0 Ncon+Npat+2+1],'box','off');
          ylabel('LEARNING INDEX');
          
          saveas(gcf,[figdir,upper(indname), '_bar_learningindex'],'fig');
          preprint;
          eval(['print -depsc ',epsdir,upper(indname),'_bar_learningindex'])
      end      
    end
    %pause
    %close all 
    
    
end
