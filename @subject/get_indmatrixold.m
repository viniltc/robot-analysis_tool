function [ind_bydir,ind_bydir_catch,all_inds] = get_indmatrix(subj,indicator)
% gets matrix of values for indicator INDNAME, for all target sets of that subject
% (C) V. Sanguineti, 2008

all_inds = [];
ntsets = length([subj.order{:}]);

for tsno = 1:ntsets
        [indmat,indmat_catch,inds]=get_indmatrix(subj.tset{tsno},indicator);
        ind_bydir{1,tsno} = indmat;  % this is Nreps x Ntargets
        ind_bydir_catch{1,tsno} = indmat_catch;
        all_inds = [all_inds inds];
end