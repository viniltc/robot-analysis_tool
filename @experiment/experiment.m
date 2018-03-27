% function exp=experiment(subjects,protocol, task, varargin)
% 
% % EXPERIMENT() Constructor of EXPERIMENT object
% % An EXPERIMENT object stores information about one specific experiment,
% % involving one experimental protocol and multiple SUBJECTS
% % See also: SUBJECT
% % (C) V. Sanguineti 2008-2015
% 
% exp.subjnames=subjects.names;
% exp.datadir=subjects.datadir;
% exp.groups=subjects.groups;
% 
% % keyboard;
% exp.protocol=protocol;
% exp.task = task; 
% 
% forces = [];
% dirs = [];
% rotations = [];
% 
% for s = 1:length(subjects.names)
%     subjname=subjects.names{s}
%     switch(nargin)
%         case 3, 
%         exp.subj(s)=subject(subjname,exp.datadir,protocol,task);
%         case 4,
%         exp.subj(s)=subject(subjname,exp.datadir,protocol,task,varargin{1});
%     end
%     
%     forces = [forces;get_forces(exp.subj(s))];
%     dirs = [dirs;get_dirs(exp.subj(s))];
%     rotations = [rotations;get_rotations(exp.subj(s))];
% 
% end
% exp.dirs=dirs;
% exp.forces=forces;
% exp.rotations=rotations;
% 
% exp = class(exp,'experiment');


function [scores,exp]=experiment(subjects,protocol, task, varargin)

% EXPERIMENT() Constructor of EXPERIMENT object
% An EXPERIMENT object stores information about one specific experiment,
% involving one experimental protocol and multiple SUBJECTS
% See also: SUBJECT
% (C) V. Sanguineti 2008-2015

exp.subjnames=subjects.names;
exp.datadir=subjects.datadir;
exp.groups=subjects.groups;

%keyboard;
exp.protocol=protocol;
exp.task = task; 

% keyboard

forces = [];
dirs = [];
rotations = [];
scores=[];
for s = 1:length(subjects.names)
    subjname=subjects.names{s}
    switch(nargin)
        case 3, 
        [tset1,exp.subj(s)]=subject(subjname,exp.datadir,protocol,task);
        case 4,
        [tset1,exp.subj(s)]=subject(subjname,exp.datadir,protocol,task,varargin{1});
    end
    scores=[scores; tset1];   
    forces = [forces;get_forces(exp.subj(s))];
    dirs = [dirs;get_dirs(exp.subj(s))];
    rotations = [rotations;get_rotations(exp.subj(s))];

end
exp.dirs=dirs;
exp.forces=forces;
exp.rotations=rotations;

exp = class(exp,'experiment');
 
%  