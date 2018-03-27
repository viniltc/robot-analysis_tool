function [indmat,indmat_catch,inds]=get_indmatrix(tset,indname)
% gets matrix of values for indicator INDNAME, for that target set
% (C) V. Sanguineti, 2008

inds = get_indicator(tset,indname)';
indmat = NaN*ones(1,tset.Ntargets);
indmat_catch = NaN*ones(1,tset.Ntargets);

forces = get_forces(tset);
force_trial = any(forces);%
rotations = get_rotations(tset);
rotation_trial = any(rotations);%

% from now on work in progress to manage rotation ct's....

Nrepmax = tset.Ntrials./tset.Ntargets;
tgt_cnt=zeros(tset.Ntargets,1);
tgt_cnt_catch=zeros(tset.Ntargets,1);
%indmat = zeros(Nrepmax,tset.Ntargets);
for trial = 1:tset.Ntrials
    tgno = tset.targetno(trial);
    if force_trial
        if forces(trial)
            tgt_cnt(tgno)=tgt_cnt(tgno)+1;
            indmat(tgt_cnt(tgno),tgno)=inds(trial);
        else
            tgt_cnt_catch(tgno)=tgt_cnt_catch(tgno)+1;
            indmat_catch(tgt_cnt_catch(tgno),tgno)=inds(trial);
        end
    else
            tgt_cnt(tgno)=tgt_cnt(tgno)+1;
            indmat(tgt_cnt(tgno),tgno)=inds(trial);      
    end
  
end
    

