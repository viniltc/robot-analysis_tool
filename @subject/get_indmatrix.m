function [ind_bydir,ind_bydir_catch,all_inds,ind_bydir_nofor,ind_bydir_norot] = get_indmatrix(subj,indicator)
% gets matrix of values for indicator INDNAME, for all target sets of that subject
% (C) V. Sanguineti, 2008

all_inds = [];
ntsets = length([subj.order{:}]);

for tsno = 1:ntsets
        [indmat,indmat_catch, inds, indmat_catch_nofor, indmat_catch_norot]=get_indmatrix(subj.tset{tsno},indicator);
        ind_bydir{1,tsno} = indmat;  % this is Nreps x Ntargets
        ind_bydir_catch{1,tsno} = indmat_catch;
        all_inds = [all_inds inds];
            ind_bydir_nofor{1,tsno} = indmat_catch_nofor;
            ind_bydir_norot{1,tsno} = indmat_catch_norot;
end