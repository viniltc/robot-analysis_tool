% function subj=subject(subjname,pathname,protocol,task,varargin)
% % Constructor of a SUBJECT object (set of target sets for a gine subject) 
% % One subject is characterized by: 
% %   subject name
% %   pathname
% %   protocol pars
% %   task pars
% %   varargin is for legacy data
% % (C) V. Sanguineti (2008-2015)
% 
% if isfield(protocol,'data_type')
%     DATATYPE = protocol.data_type;
% else
%     DATATYPE = 'H3D'; % this is the default
% end
% 
% if ismac | isunix
%     sepstr = '/';
% elseif ispc
%     sepstr='\';
% end
% 
% 
% subj.name = subjname;
% subj.pathname = pathname;
% 
% 
% subj.phases = protocol.phases;
% for p=1:length(protocol.phases)
%     subj.order{p} = 1:length(protocol.order{p});
% end
% subj.data_format = protocol.data_format;
% 
% 
% tsets = [protocol.order{:}];
% 
% forces = [];
% dirs = [];
% rotations = [];
% for tsno  = 1:length(tsets)  
%     tsno
%     % Change lines below if different data format
%     switch DATATYPE
%         case 'H3D', 
%             nomefile = sprintf('%s_%1d.dat',protocol.dataprefix,tsets(tsno));
%         case 'RTLAB',
%             nomefile = sprintf('%s_%1d.mat',protocol.dataprefix,tsets(tsno));        
%         case 'CHAI3D',
%             nomefile  = sprintf('%s_%1d.dat',protocol.dataprefix,tsets(tsno));
%             nomefile1 = sprintf('%s_%1d.pro',protocol.dataprefix,tsets(tsno));%%%%%%
%     
%     end
%     str = sprintf('%s%s',pathname,subjname,sepstr,nomefile)
%     str1 = sprintf('%s%s',pathname,subjname,sepstr,nomefile1)%%%%%%
%     
%     
%     
%     switch DATATYPE
%         case 'CHAI3D',  % This format is used for dyad data
%             fid = fopen(str,'r');
%             data = fread(fid,[14 inf],'double')';
%             fclose(fid);
%             data = data';
%             
%             %%%%
%             fid1= fopen(str1,'r');
%             data1 = fread(fid1,[3 inf],'double')';
%             fclose(fid1);
%             data1=data1';
%        
%         case 'H3D',   
%             data = load(str);
%             % merge timestamp columns (assumes they come at first 2
%             % columns)
%            % keyboard
%            % In new versions of Humour save node, time is stored in the
%            % first column only
%            if(data(2,1)==fix(data(2,1)))
%             data(:,1) = data(:,1)+data(:,2);
%            else
%             data=[data(:,1) data(:,1) data(:,2:end)];   
%            end
%             data = data';%(:,[1 3:end])'; % take out second time field and transpose
% 
%             case 'RTLAB',
%             load(str)
%             nomedata=who('-file', str);
%             eval(['data=', nomedata{1},';']);
%     end
% 
%     
%     switch(task.type)
%         case 'bimanual',
%             if nargin==4
%                 tset = targetset(data,protocol.data_format,task,{'Ltraj', 'Rtraj'});       
%             else
%                 tset = targetset(data,protocol.data_format,task,{'Ltraj', 'Rtraj'},varargin{1});       
%             end
%         case 'unimanual', 
%             if nargin==4
%                 tset = targetset(data,protocol.data_format,task,'traj');       
%             else
%                 tset = targetset(data,protocol.data_format,task,'traj',varargin{1});
%             end
%         case 'dyad', 
%             if nargin==4
%                 tset = targetset(data,protocol.data_format,task,{'traj1','traj2'});
%                 
%                 %%%
%                 %tset1 = targetset(data1,protocol.data_format1,task,{'traj1','traj2'});
%             else
%                 tset = targetset(data,protocol.data_format,task,{'traj1','traj2'},varargin{1});
%                 
%                %%%
%                % tset1 = targetset(data1,protocol.data_format1,task,{'traj1','traj2'},varargin{1});
%             end
%             
%             
%             
%     end    
%     forces = [forces get_forces(tset)];
%     dirs = [dirs get_dirs(tset)];
%     rotations = [rotations get_rotations(tset)];
%     subj.tset{tsno}=tset;
%     
% end    
% subj.forces = forces; 
% subj.dirs = dirs; 
% subj.rotations = rotations; 
% 
% 
% subj=class(subj,'subject');

function [tset1,subj]=subject(subjname,pathname,protocol,task,varargin)
% Constructor of a SUBJECT object (set of target sets for a gine subject) 
% One subject is characterized by: 
%   subject name
%   pathname
%   protocol pars
%   task pars
%   varargin is for legacy data
% (C) V. Sanguineti (2008-2015)

if isfield(protocol,'data_type')
    DATATYPE = protocol.data_type;
else
    DATATYPE = 'H3D'; % this is the default
end

if ismac | isunix
    sepstr = '/';
elseif ispc
    sepstr='\';
end


subj.name = subjname;
subj.pathname = pathname;


subj.phases = protocol.phases;
for p=1:length(protocol.phases)
    subj.order{p} = 1:length(protocol.order{p});
end
subj.data_format = protocol.data_format;

tsets = [protocol.order{:}];

forces = [];
dirs = [];
rotations = [];
tset1 = [];
for tsno  = 1:length(tsets)  
    tsno
    % Change lines below if different data format
    switch DATATYPE
        case 'H3D', 
            nomefile = sprintf('%s_%1d.dat',protocol.dataprefix,tsets(tsno));
        case 'RTLAB',
            nomefile = sprintf('%s_%1d.mat',protocol.dataprefix,tsets(tsno));        
        case 'CHAI3D',
            nomefile  = sprintf('%s_%1d.dat',protocol.dataprefix,tsets(tsno));
            nomefile1 = sprintf('%s_%1d.pro',protocol.dataprefix,tsets(tsno));%%%%%%
    
    end
    str  =  sprintf('%s%s',pathname,subjname,sepstr,nomefile)
    str1 =  sprintf('%s%s',pathname,subjname,sepstr,nomefile1)
    
%     keyboard 
    
    switch DATATYPE
        case 'CHAI3D',  % This format is used for dyad data
            fid = fopen(str,'r');
            data = fread(fid,[14 inf],'double')';
            fclose(fid);
            data = data';

            fid1=fopen(str1,'r');
            data1 = fscanf(fid1,'%d', [5,12])'; %
            fclose(fid1);
            data1 = data1(:,1:3); %sliced data , score1, score2 over trial
            data1 = data1';
            
         
            
        case 'H3D',   
            data = load(str);
            % merge timestamp columns (assumes they come at first 2
            % columns)
           % keyboard
           % In new versions of Humour save node, time is stored in the
           % first column only
           if(data(2,1)==fix(data(2,1)))
            data(:,1) = data(:,1)+data(:,2);
           else
            data=[data(:,1) data(:,1) data(:,2:end)];   
           end
            data = data';%(:,[1 3:end])'; % take out second time field and transpose

            case 'RTLAB',
            load(str)
            nomedata=who('-file', str);
            eval(['data=', nomedata{1},';']);
    end

    
    switch(task.type)
        case 'bimanual',
            if nargin==4
                tset = targetset(data,protocol.data_format,task,{'Ltraj', 'Rtraj'});       
            else
                tset = targetset(data,protocol.data_format,task,{'Ltraj', 'Rtraj'},varargin{1});       
            end
        case 'unimanual', 
            if nargin==4
                tset = targetset(data,protocol.data_format,task,'traj');       
            else
                tset = targetset(data,protocol.data_format,task,'traj',varargin{1});
            end
        case 'dyad', 
            if nargin==4
                tset = targetset(data,protocol.data_format,task,{'traj1','traj2'});              
                %tset1 = targetset(data1,protocol.data_format1,task,{'traj1','traj2'}); % dati .pro
            else
                tset = targetset(data,protocol.data_format,task,{'traj1','traj2'},varargin{1});
                %tset1 = targetset(data1,protocol.data_format1,task,{'traj1','traj2'},varargin{1}); % dati .pro
            end
            
            
            
    end    
    forces = [forces get_forces(tset)];
    dirs = [dirs get_dirs(tset)];
    rotations = [rotations get_rotations(tset)];
    subj.tset{tsno}=tset;
    tset1=[tset1 data1];
    
end    
subj.forces = forces; 
subj.dirs = dirs; 
subj.rotations = rotations; 


subj=class(subj,'subject');


