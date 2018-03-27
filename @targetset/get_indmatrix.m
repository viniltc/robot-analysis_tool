function [indmat,indmat_catch, inds, indmat_catch_nofor, indmat_catch_norot]=get_indmatrix(tset,indname)
% gets matrix of values for indicator INDNAME, for that target set
% (C) V. Sanguineti, 2008

inds = get_indicator(tset,indname)';
indmat = NaN*ones(1,tset.Ntargets);
indmat_catch = NaN*ones(1,tset.Ntargets);
indmat_catch_nofor= NaN*ones(1,tset.Ntargets);
indmat_catch_norot= NaN*ones(1,tset.Ntargets); 

forces = get_forces(tset);

force_trial = any(forces)%
rotations = get_rotations(tset);
rotation_trial = not(sum(rotations+forces)==36);%any(rotations)%



% from now on work in progress to manage rotation ct's....

Nrepmax = tset.Ntrials./tset.Ntargets;
tgt_cnt=zeros(tset.Ntargets,1);
tgt_cnt_catch=zeros(tset.Ntargets,1);
tgt_cnt_catch_norot=zeros(tset.Ntargets,1);
tgt_cnt_catch_nofor=zeros(tset.Ntargets,1);
% indmat = zeros(Nrepmax,tset.Ntargets);

% % % % for trial = 1:tset.Ntrials
% % % %     tgno = tset.targetno(trial);
% % % %     if force_trial & rotation_trial
% % % %         %keyboard
% % % %         if (forces(trial) & rotations(trial))
% % % %             tgt_cnt(tgno)=tgt_cnt(tgno)+1;
% % % %             indmat(tgt_cnt(tgno),tgno)=inds(trial);
% % % %         elseif forces(trial) & ~rotations(trial)  %%ct sulla rot!
% % % %             tgt_cnt_catch_norot(tgno)=tgt_cnt_catch_norot(tgno)+1;
% % % %             indmat_catch_norot(tgt_cnt_catch_norot(tgno),tgno)=inds(trial);
% % % %           %  keyboard;
% % % %         elseif ~forces(trial) & rotations(trial) %%ct sulla for!
% % % %           %  keyboard;
% % % %             tgt_cnt_catch_nofor(tgno)=tgt_cnt_catch_nofor(tgno)+1;
% % % %             indmat_catch_nofor(tgt_cnt_catch_nofor(tgno),tgno)=inds(trial);
% % % %         else
% % % %             tgt_cnt_catch(tgno)=tgt_cnt_catch(tgno)+1;
% % % %             indmat_catch(tgt_cnt_catch(tgno),tgno)=inds(trial);
% % % %         end
% % % %     elseif force_trial & ~rotation_trial %nn c entra mai?
% % % %         if forces(trial)
% % % %             tgt_cnt(tgno)=tgt_cnt(tgno)+1;
% % % %             indmat(tgt_cnt(tgno),tgno)=inds(trial);
% % % %         else
% % % %             tgt_cnt_catch(tgno)=tgt_cnt_catch(tgno)+1;
% % % %             indmat_catch(tgt_cnt_catch(tgno),tgno)=inds(trial);
% % % %         end   
% % % %     elseif ~force_trial & rotation_trial  %tutte tranne dalla 6 alla 13
% % % %         trial
% % % %         tgno 
% % % %         tgt_cnt
% % % %         keyboard
% % % %         if rotations(trial)
% % % %             tgt_cnt(tgno)=tgt_cnt(tgno)+1;
% % % %             indmat(tgt_cnt(tgno),tgno)=inds(trial);
% % % %         else
% % % %             tgt_cnt_catch(tgno)=tgt_cnt_catch(tgno)+1;
% % % %             indmat_catch(tgt_cnt_catch(tgno),tgno)=inds(trial);
% % % %         end
% % % %     else    
% % % %         tgt_cnt(tgno)=tgt_cnt(tgno)+1;
% % % %         indmat(tgt_cnt(tgno),tgno)=inds(trial);      
% % % %     end
% % % %      
% % % %   
% % % % end
% keyboard;
for trial = 1:tset.Ntrials
    tgno = tset.targetno(trial);
    if force_trial & rotation_trial
        if (forces(trial) & rotations(trial))
            tgt_cnt(tgno)=tgt_cnt(tgno)+1;
            indmat(tgt_cnt(tgno),tgno)=inds(trial);
        elseif forces(trial) & ~rotations(trial)  %%ct sulla rot!
            tgt_cnt_catch_norot(tgno)=tgt_cnt_catch_norot(tgno)+1;
            indmat_catch_norot(tgt_cnt_catch_norot(tgno),tgno)=inds(trial);
          %  keyboard;
        elseif ~forces(trial) & rotations(trial) %%ct sulla for!
          %  keyboard;
            tgt_cnt_catch_nofor(tgno)=tgt_cnt_catch_nofor(tgno)+1;
            indmat_catch_nofor(tgt_cnt_catch_nofor(tgno),tgno)=inds(trial);
        else
            tgt_cnt_catch(tgno)=tgt_cnt_catch(tgno)+1;
            indmat_catch(tgt_cnt_catch(tgno),tgno)=inds(trial);
        end
    end
    if ~force_trial & rotation_trial %TUTTE TRANNE BLOCCO 3 DALLA 6 ALLA 13
        if rotations(trial)
            tgt_cnt(tgno)=tgt_cnt(tgno)+1;
            indmat(tgt_cnt(tgno),tgno)=inds(trial);
        else
            tgt_cnt_catch(tgno)=tgt_cnt_catch(tgno)+1;
            indmat_catch(tgt_cnt_catch(tgno),tgno)=inds(trial);
        end  
    end
   
    if ~force_trial & ~rotation_trial  %
            tgt_cnt(tgno)=tgt_cnt(tgno)+1;
            indmat(tgt_cnt(tgno),tgno)=inds(trial);
    end

     
  
end    
% keyboard;
