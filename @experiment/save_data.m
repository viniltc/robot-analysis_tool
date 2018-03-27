% function save_data(exp,resdir,indlist)
% % SAVE_DATA  Saves on file data from one EXPERIMENT
% % An EXPERIMENT object stores information about one specific experiment,
% % involving one experimental protocol and multiple SUBJECTS
% % See also: SUBJECT, EXPERIMENT
% % (C) V. Sanguineti 2008
% 
% forces = exp.forces;
% dirs = exp.dirs;
% rotations = exp.rotations;
% 
% save([resdir,'forces.mat'],'forces');
% save([resdir,'dirs.mat'],'dirs');
% save([resdir,'rotations.mat'],'rotations');
% 
% ntsets = length([exp.protocol.order{:}]);
% nsubj  = length(exp.subj);
% Ntargets = size(exp.task.targets,1);
% 
% for ind = 1:length(indlist)
%         indname = indlist{ind};
%         for s = 1:nsubj
%             [indmat,indmat_catch,inds,indmat_catch_nofor, indmat_catch_norot]=get_indmatrix(exp.subj(s),indname);
%             ind_bydir(s,:) = indmat;  % this is Nreps x Ntargets
%             ind_bydir_catch(s,:) = indmat_catch;
%             ind_ts(s,:) = inds;
%             ind_bydir_nofor(s,:) = indmat_catch_nofor;
%             ind_bydir_norot(s,:) = indmat_catch_norot;
%             
%         end
% 
%         eval(['save(''',resdir,indname,'.mat'',''ind_ts'');']);
%         eval(['save(''',resdir,indname,'_bydir.mat'',''ind_bydir'');']);
%         eval(['save(''',resdir,indname,'_bydir_catch.mat'',''ind_bydir_catch'');']);  
%         eval(['save(''',resdir,indname,'_bydir_nofor.mat'',''ind_bydir_nofor'');']);
%         eval(['save(''',resdir,indname,'_bydir_norot.mat'',''ind_bydir_norot'');']);
% 
% 
%         
%         for tset = 1:ntsets
%          for subj=1:nsubj
%             mind(subj,tset,:) = mean(ind_bydir{subj,tset});
%            % keyboard;
%             %phase=1;
%             %while((tset-1)>exp.protocol.phases{phase}(end))
%             %    phase=phase+1;
%             %end
%             %if not(strfind(exp.protocol.phases{phase},'null'))
%                 mind_catch(subj,tset,:) = nanmean(ind_bydir_catch{subj,tset}); 
%                 mind_nofor(subj,tset,:) = nanmean(ind_bydir_nofor{subj,tset});
%                 mind_norot(subj,tset,:) = nanmean(ind_bydir_norot{subj,tset});
%             %end
%          end
%         end
% 
%         
%         mind_statmat = exp.groups;
%         mind_statmat_catch = exp.groups;
%         mind_statmat_norot = exp.groups;
%         mind_statmat_nofor = exp.groups;
%         
%         for tset = 1:ntsets
%            mind_statmat = [mind_statmat reshape(mind(:,tset,:),nsubj,Ntargets)];
%            %phase=1;
%            %while((tset-1)>exp.protocol.phases{phase}(end))
%            %     phase=phase+1;
%            %end
%            %if not(strfind(exp.protocol.phases{phase},'null'))
%                mind_statmat_catch = [mind_statmat_catch reshape(mind_catch(:,tset,:),nsubj,Ntargets)];
%                mind_statmat_nofor = [mind_statmat_nofor reshape(mind_nofor(:,tset,:),nsubj,Ntargets)];
%                mind_statmat_norot = [mind_statmat_norot reshape(mind_norot(:,tset,:),nsubj,Ntargets)];
%            %end
% 
%        
%         end
%         
% %          for i=1:nsubj
% %   % plot_learning_index(exp.subj(i),indname,'learning index')
% % 
% %  % plot_catch_force(exp.subj(i),indname,indname)
% %   % displ_list: list of graph types
% % % figdir: directory where figs are saved
% % % epsdir: directory where eps figs are saved
% % %figdir='\data\figs'
% %    %plot(exp.subj(i),indlabel,figdir,epsdir)
% %  
% % %  plot_adaptation(exp.subj(i),indname,indname)
% %          end
%         fname = [resdir, indname,'_statmat.txt'];
%         dlmwrite(fname,mind_statmat,'\t');
%         
%         fname_c = [resdir, indname,'_statmat_catch.txt'];
%         dlmwrite(fname_c,mind_statmat_catch,'\t');
%         
%         fname_c = [resdir, indname,'_statmat_nofor.txt'];
%         dlmwrite(fname_c,mind_statmat_nofor,'\t');
%         
%         fname_c = [resdir, indname,'_statmat_norot.txt'];
%         dlmwrite(fname_c,mind_statmat_norot,'\t');
%         
%  end
 
     

function save_data(scores,exp,resdir,indlist)
% SAVE_DATA  Saves on file data from one EXPERIMENT
% An EXPERIMENT object stores information about one specific experiment,
% involving one experimental protocol and multiple SUBJECTS
% See also: SUBJECT, EXPERIMENT
% (C) V. Sanguineti 2008

forces = exp.forces;
dirs = exp.dirs;
rotations = exp.rotations;

save([resdir,'forces.mat'],'forces');
save([resdir,'dirs.mat'],'dirs');
save([resdir,'rotations.mat'],'rotations');
save([resdir,'scores.mat'],'scores');

ntsets = length([exp.protocol.order{:}]);
nsubj  = length(exp.subj);
Ntargets = size(exp.task.targets,1);

for ind = 1:length(indlist)
        indname = indlist{ind};
        for s = 1:nsubj
            [indmat,indmat_catch,inds,indmat_catch_nofor, indmat_catch_norot]=get_indmatrix(exp.subj(s),indname);
            ind_bydir(s,:) = indmat;  % this is Nreps x Ntargets
            ind_bydir_catch(s,:) = indmat_catch;
            ind_ts(s,:) = inds;
            ind_bydir_nofor(s,:) = indmat_catch_nofor;
            ind_bydir_norot(s,:) = indmat_catch_norot;
            
        end

        eval(['save(''',resdir,indname,'.mat'',''ind_ts'');']);
        eval(['save(''',resdir,indname,'_bydir.mat'',''ind_bydir'');']);
        eval(['save(''',resdir,indname,'_bydir_catch.mat'',''ind_bydir_catch'');']);  
        eval(['save(''',resdir,indname,'_bydir_nofor.mat'',''ind_bydir_nofor'');']);
        eval(['save(''',resdir,indname,'_bydir_norot.mat'',''ind_bydir_norot'');']);


        
        for tset = 1:ntsets
         for subj=1:nsubj
            mind(subj,tset,:) = nanmean(ind_bydir{subj,tset});
           % keyboard;
            %phase=1;
            %while((tset-1)>exp.protocol.phases{phase}(end))
            %    phase=phase+1;
            %end
            %if not(strfind(exp.protocol.phases{phase},'null'))
                mind_catch(subj,tset,:) = nanmean(ind_bydir_catch{subj,tset}); 
                mind_nofor(subj,tset,:) = nanmean(ind_bydir_nofor{subj,tset});
                mind_norot(subj,tset,:) = nanmean(ind_bydir_norot{subj,tset});
            %end
         end
        end

        
        mind_statmat = exp.groups;
        mind_statmat_catch = exp.groups;
        mind_statmat_norot = exp.groups;
        mind_statmat_nofor = exp.groups;
        
        for tset = 1:ntsets
           mind_statmat = [mind_statmat reshape(mind(:,tset,:),nsubj,Ntargets)];
           %phase=1;
           %while((tset-1)>exp.protocol.phases{phase}(end))
           %     phase=phase+1;
           %end
           %if not(strfind(exp.protocol.phases{phase},'null'))
               mind_statmat_catch = [mind_statmat_catch reshape(mind_catch(:,tset,:),nsubj,Ntargets)];
               mind_statmat_nofor = [mind_statmat_nofor reshape(mind_nofor(:,tset,:),nsubj,Ntargets)];
               mind_statmat_norot = [mind_statmat_norot reshape(mind_norot(:,tset,:),nsubj,Ntargets)];
           %end

       
        end
        
%          for i=1:nsubj
%   % plot_learning_index(exp.subj(i),indname,'learning index')
% 
%  % plot_catch_force(exp.subj(i),indname,indname)
%   % displ_list: list of graph types
% % figdir: directory where figs are saved
% % epsdir: directory where eps figs are saved
% %figdir='\data\figs'
%    %plot(exp.subj(i),indlabel,figdir,epsdir)
%  
% %  plot_adaptation(exp.subj(i),indname,indname)
%          end
        fname = [resdir, indname,'_statmat.txt'];
        dlmwrite(fname,mind_statmat,'\t');
        
        fname_c = [resdir, indname,'_statmat_catch.txt'];
        dlmwrite(fname_c,mind_statmat_catch,'\t');
        
        fname_c = [resdir, indname,'_statmat_nofor.txt'];
        dlmwrite(fname_c,mind_statmat_nofor,'\t');
        
          fname_c = [resdir, indname,'_statmat_norot.txt'];
        dlmwrite(fname_c,mind_statmat_norot,'\t');
        
 end
 
     