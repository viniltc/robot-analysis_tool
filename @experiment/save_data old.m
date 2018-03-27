function save_data(exp,resdir,indlist)
% SAVE_DATA  Saves on file data from one EXPERIMENT
% An EXPERIMENT object stores information about one specific experiment,
% involving one experimental protocol and multiple SUBJECTS
% See also: SUBJECT, EXPERIMENT
% (C) V. Sanguineti 2008

forces = exp.forces;
dirs = exp.dirs;
save([resdir,'forces.mat'],'forces');
save([resdir,'dirs.mat'],'dirs');

ntsets = length([exp.protocol.order{:}]);
nsubj  = length(exp.subj);
Ntargets = size(exp.task.targets,1);

for ind = 1:length(indlist)
        indname = indlist{ind};
        for s = 1:nsubj
            [indmat,indmat_catch,inds]=get_indmatrix(exp.subj(s),indname);
            ind_bydir(s,:) = indmat;  % this is Nreps x Ntargets
            ind_bydir_catch(s,:) = indmat_catch;
            ind_ts(s,:) = inds;
        end

        eval(['save(''',resdir,indname,'.mat'',''ind_ts'');']);
        eval(['save(''',resdir,indname,'_bydir.mat'',''ind_bydir'');']);
        eval(['save(''',resdir,indname,'_bydir_catch.mat'',''ind_bydir_catch'');']);  
        
        for tset = 1:ntsets
         for subj=1:nsubj
            mind(subj,tset,:) = mean(ind_bydir{subj,tset});
            mind_catch(subj,tset,:) = (ind_bydir_catch{subj,tset});

         end
        end

        
        mind_statmat = exp.groups;
        mind_statmat_catch = exp.groups;
        
        for tset = 1:ntsets
           mind_statmat = [mind_statmat reshape(mind(:,tset,:),nsubj,Ntargets)];
           mind_statmat_catch = [mind_statmat_catch reshape(mind_catch(:,tset,:),nsubj,Ntargets)];
       
        end

        fname = [resdir, '\', indname,'_statmat.txt'];
        dlmwrite(fname,mind_statmat,'\t');
        
        fname_c = [resdir, '\', indname,'_statmat_catch.txt'];
        dlmwrite(fname_c,mind_statmat_catch,'\t');
        
 end
 
     