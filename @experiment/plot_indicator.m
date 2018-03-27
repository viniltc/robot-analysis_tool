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
    
    % takes average over dirs (only if ndirs>1)
    for s = 1:Nsu % loop over subjects...
        Ncnt = 0;
        Ngrp = 1;  % 1 points per each tset
        
        for ts = 1:Nts % loop over target sets...
          
 
          
          % take out outliers
          o_med = nanmedian(indmat_bydir{s,ts});
          o_std = nanstd(indmat_bydir{s,ts});
          [is_outl1,is_outl2] = find(abs(indmat_bydir{s,ts} - o_med(ones(size(indmat_bydir{s,ts},1),1),:))> ...
                                     3*(o_std(ones(size(indmat_bydir{s,ts},1),1),:)));
          indmat_bydir{s,ts}(is_outl1,is_outl2)=NaN*indmat_bydir{s,ts}(is_outl1,is_outl2);
         
          
          indsub = abs(indmat_bydir{s,ts}); % this is used for adaptation plots: size is nrep X ndirs
          
          % Creates a binned version of indsub: size is now 1 x Ngrp
          Nreps = size(indsub,1);
          Navg = round(Nreps/Ngrp);
          if size(indsub,2)>1
             indavg(s,Ncnt+(1:Ngrp)) = nanmean(reshape((nanmean(indsub')),Navg,Ngrp)); 
          else
             indavg(s,Ncnt+(1:Ngrp)) = nanmean(reshape(indsub,Navg,Ngrp)); 
          end
          
          % This part is to calculate learning index
          if ~isempty(indcatch) & indname == indcatch  
             
             % this is error at last phases (10:12) of null field: size is 1 x ndirs  
             meanlast = nanmedian((indmat_bydir{s,lastnull}(end+(-1:0),:)));  
            
             % this is like indsub, but with sign...: size is nrep x ndirs
             indsub_noa = (indmat_bydir{s,ts})-meanlast(ones(Nreps,1),:);
             indsub_noa_catch = indmat_bydir_catch{s,ts}(:,1)-meanlast(ones(size(indmat_bydir_catch{s,ts},1),1),1);
             indsub_noa_catchnofor=ind_bydir_nofor{s,ts}(:,1)-meanlast(ones(size(ind_bydir_nofor{s,ts},1),1),1);
             indsub_noa_catchnorot=ind_bydir_norot{s,ts}(:,1)-meanlast(ones(size(ind_bydir_norot{s,ts},1),1),1);
            
             % Now take indsub_noa and bin it (Ngrp is bin size) 
             % Size is now 1 x Ngrp...
             indavg_noa(s,Ncnt+(1:Ngrp)) = nanmean(reshape(nanmean(indsub_noa'),Navg,size(indsub_noa,1)/Navg));
             indavg_noa_catch(s,Ncnt+(1:Ngrp)) = nanmean(indsub_noa_catch');
             indavg_noa_catchnofor(s,Ncnt+(1:Ngrp)) = nanmean(indsub_noa_catchnofor');
             indavg_noa_catchnorot(s,Ncnt+(1:Ngrp)) = nanmean(indsub_noa_catchnorot');
             
%              
%              if isnan(nanmean(indsub_noa_catch'))
%                  indavg_noa_catch(s,Ncnt+(1:Ngrp)) = 0;%indavg_noa(s,Ncnt+(1:Ngrp));
%              end
%              
%                if isnan(nanmean(indsub_noa_catchnofor'))
%                  indavg_noa_catchnofor(s,Ncnt+(1:Ngrp)) =  indavg_noa(s,Ncnt+(1:Ngrp));
%                end
%              
%                  if isnan(nanmean(indsub_noa_catchnorot'))
%                  indavg_noa_catchnorot(s,Ncnt+(1:Ngrp)) = indavg_noa(s,Ncnt+(1:Ngrp));
%              end
% %             
            
            li_noa(s,Ncnt+(1:Ngrp))=abs(indavg_noa_catch(s,Ncnt+(1:Ngrp)))./...
                                    abs(indavg_noa(s,Ncnt+(1:Ngrp))-indavg_noa_catch(s,Ncnt+(1:Ngrp)));
            li_noanf(s,Ncnt+(1:Ngrp))=abs(indavg_noa_catchnofor(s,Ncnt+(1:Ngrp)))./...
                                    abs(indavg_noa(s,Ncnt+(1:Ngrp))-indavg_noa_catchnofor(s,Ncnt+(1:Ngrp)));
            li_noanr(s,Ncnt+(1:Ngrp))=abs(indavg_noa_catchnorot(s,Ncnt+(1:Ngrp)))./...
                                    abs(indavg_noa(s,Ncnt+(1:Ngrp))-indavg_noa_catchnorot(s,Ncnt+(1:Ngrp)));

          end
          Ncnt = Ncnt+Ngrp;
          
              
        end  % of tset
    end % of subject
    
    inc = nanmean(indavg(find(exp.groups==1),:));
    s_inc = nanstd(indavg(find(exp.groups==1),:));
    
    if ~isempty(indcatch) &  indname == indcatch
     inc_noa = nanmean(indavg_noa(find(exp.groups==1),:));
     s_inc_noa = nanstd(indavg_noa(find(exp.groups==1),:));
     
     inc_noa_catch = nanmean(indavg_noa_catch(find(exp.groups==1),:));
     s_inc_noa_catch = nanstd(indavg_noa_catch(find(exp.groups==1),:));
     
     inc_noa_catchnofor = nanmean(indavg_noa_catchnofor(find(exp.groups==1),:));
     s_inc_noa_catchfor = nanstd(indavg_noa_catchnofor(find(exp.groups==1),:));
     
     
     inc_noa_catchnorot = nanmean(indavg_noa_catchnorot(find(exp.groups==1),:));
     s_inc_noa_catchnorot = nanstd(indavg_noa_catchnorot(find(exp.groups==1),:));
     
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
    
     if ~isempty(indcatch) & indname == indcatch
    
         inp_noa = nanmean(indavg_noa(find(exp.groups==2),:));
         s_inp_noa = nanstd(indavg_noa(find(exp.groups==2),:));
        
         inp_noa_catch = nanmean(indavg_noa_catch(find(exp.groups==2),:));
         s_inp_noa_catch = nanstd(indavg_noa_catch(find(exp.groups==2),:));
         
         inp_noa_catchnorot = nanmean(indavg_noa_catchnorot(find(exp.groups==2),:));
         s_inp_noa_catchnorot = nanstd(indavg_noa_catchnorot(find(exp.groups==2),:));
       
         inp_noa_catchnofor = nanmean(indavg_noa_catchnofor(find(exp.groups==2),:));
         s_inp_noa_catchnofor = nanstd(indavg_noa_catchnofor(find(exp.groups==2),:));
         
         inp_li =nanmean(li_noa(find(exp.groups==2),:));
         s_inp_li =nanstd(li_noa(find(exp.groups==2),:)); 
         
                  inp_linr =nanmean(li_noanr(find(exp.groups==2),:));
         s_inp_linr =nanstd(li_noanr(find(exp.groups==2),:)); 
         
                  inp_linf =nanmean(li_noanf(find(exp.groups==2),:));
         s_inp_linf =nanstd(li_noanf(find(exp.groups==2),:)); 
    
     end
     
     h(2)=line(1:length(inp),inp,'col','r','LineWidth',2);%,'lines','--' );
     line(1:length(inp),inp-s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
     line(1:length(inp),inp+s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
     
     
%       line(1:length(inp),inp-s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
%      line(1:length(inp),inp+s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
%      
%           line(1:length(inp),inp-s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
%      line(1:length(inp),inp+s_inp./sqrt(Npat),'col','r','LineWidth',1);%,'lines','--')
     
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
    if ~isempty(indcatch) & indname == indcatch
      figure
      
%      keyboard;
      plot(inc_noa,'b','LineWidth',2 )
      hold on;
              isfin = find(isfinite(inc_noa_catch));
              a=size(isfin);
              for x=1:a(1,2)
                  keyboard;
              i=isfin(1,x)
              
               plot(inc_noa_catch(1,i),'b','LineWidth',2, 'lines',':' )
               end
%                keyboard;
%               end
%                    isfin = find(isfinite(inc_noa_catchnofor));
%                if (isfin>0)
%       plot(inc_noa,'b','LineWidth',2 )
%       hold on, plot(inc_noa_catchnofor,'c','LineWidth',2, 'lines',':' )
%                end


% 
%                 isfin = find(isfinite(inc_noa_catchnorot));
%                if (isfin>0)
%                plot(inc_noa,'b','LineWidth',2 )
%       hold on, plot(inc_noa_catchnorot,'g','LineWidth',2, 'lines',':' )
%                end
%       
       isfin = find(isfinite(inc_noa_catchnofor));
              a=size(isfin);
              for x=1:a(1,2)

              i=isfin(1,x)
              
               plot(inc_noa_catchnofor(1,i),'c','LineWidth',2, 'lines',':' )
               end
               
                isfin = find(isfinite(inc_noa_catchnorot));
              a=size(isfin);
              for x=1:a(1,2)

              i=isfin(1,x)
              
               plot(inc_noa_catchnorot(1,i),'g','LineWidth',2, 'lines',':' )
               end

      
      if any(exp.groups==2)
         hold on, plot(inp_noa,'r','LineWidth',2 )
         
         hold on, plot(inp_noa_catch,'r','LineWidth',2, 'lines',':' )
         
          hold on, plot(inp_noa_catchnofor,'r','LineWidth',2, 'lines',':' )
          
          hold on, plot(inp_noa_catchnorot,'r','LineWidth',2, 'lines',':' )
      end
      xlabel('MOVEMENT SET')
      ylabel(indlabel{ind})
      set(gca, 'box', 'off')
      %axis square
      set(gcf, 'pos', [100 100 350 300])
      set(gcf,'name',indlabel{ind})
    
      saveas(gcf,[figdir,upper(indname),'_withcatch'],'fig');
%      preprint;
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
%      preprint;
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
