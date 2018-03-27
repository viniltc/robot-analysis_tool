function tset=targetset(data,data_format,task,trajfield,varargin)
% Constructor of a TARGETSET (block of movement trials) object
% One target set is characterized by: 
%   data: data matrix
%   data_format: data format
%   task: structure with task parameters 
%   trajfield: name of data format entry containing trajectory samples
%   varargin: function handle for function containing definition of the
%   'runs' array (used for legacy data)
% (C) V. Sanguineti (2008-2015)

tset.fc=task.fc;
tset.Ntargets = size(task.targets,1);
tset.targets = task.targets;
tset.tgsize = task.tgsize;
tset.expmode = task.expmode;
tset.type = task.type; % unimanual or bimanual or dyad
tset.viapoints = task.viapoints; % if any

% filter params
tset.forder=6;
tset.flength=17;   % 170 ms window, this is about 11 Hz

tset.forder = 4;
tset.flength=37;   % 270 ms window, this is about 7.5 Hz


target_angles = atan2(task.targets(:,2), task.targets(:,1))';

% Assuming first row in data is time...
rawtime = data(1,:)-data(1,1);

maxtime = max(rawtime);

% Resample at uniform sampling rate (in case it is not...)
newtime = 0:(1/tset.fc):maxtime;
xy_target_f = zeros(size(data,1),length(newtime));
xy_target_f(1,:) = newtime; % assuming that first row is always time...

xy_target_f(2:(end-1),:) = interp1(rawtime,data(2:(end-1),:)',newtime)';

% this is state:
xy_target_f(end,:) = interp1(rawtime,data(end,:)',newtime,'nearest')';

% Generates one tset 'field' per element in the data_format cell array
[ncomps,nsamples] = size(xy_target_f);


for comp=1:ncomps
    compname = data_format{comp};
    eval(['d.',compname ,' =xy_target_f(',num2str(comp),',:)'';']);
end

% This is to make sure that there is a 'traj' field (useful to extend to
% bimanual targetset)
% eval(['d.traj_x = d.',trajfield,'_x;']);
% eval(['d.traj_y = d.',trajfield,'_y;']);

% just in case 'force' has not been set
if ~isfield(d,'force') && isfield(d,'force_x') && isfield(d,'force_y')
    d.force = (sqrt(d.force_x.^2+d.force_y.^2))>0;
end


time = d.time-d.time(1); % time starts from first sample in the target set

switch task.type
    case 'bimanual', 
        traj = sav_golay([d.Ltraj_x d.Ltraj_y d.Rtraj_x d.Rtraj_y],tset.flength,tset.forder,0,tset.fc); % this is trajectory (x,y)
        vel = sav_golay([d.Ltraj_x d.Ltraj_y d.Rtraj_x d.Rtraj_y],tset.flength,tset.forder,1,tset.fc); % this is velocity (x,y)
        acc = sav_golay([d.Ltraj_x d.Ltraj_y d.Rtraj_x d.Rtraj_y],tset.flength,tset.forder,2,tset.fc); % this is acceleration (x,y)
        jerk = sav_golay([d.Ltraj_x d.Ltraj_y d.Rtraj_x d.Rtraj_y],tset.flength,tset.forder,3,tset.fc); % this is jerk (x,y)
        nhands = 2;
    case 'dyad',
        traj = sav_golay([d.traj1_x d.traj1_y d.traj1_z d.traj2_x d.traj2_y d.traj2_z],tset.flength,tset.forder,0,tset.fc); % this is trajectory (x,y)
        vel = sav_golay([d.traj1_x d.traj1_y d.traj1_z  d.traj2_x d.traj2_y d.traj2_z],tset.flength,tset.forder,1,tset.fc); % this is velocity (x,y)
        acc = sav_golay([d.traj1_x d.traj1_y d.traj1_z  d.traj2_x d.traj2_y d.traj2_z],tset.flength,tset.forder,2,tset.fc); % this is acceleration (x,y)
        jerk = sav_golay([d.traj1_x d.traj1_y d.traj1_z  d.traj2_x d.traj2_y d.traj2_z],tset.flength,tset.forder,3,tset.fc); % this is jerk (x,y)
        force = sav_golay([d.force1_x d.force1_y d.force1_z  d.force2_x d.force2_y d.force2_z],tset.flength,tset.forder,0,tset.fc); % this is jerk (x,y)   
        nhands = 2;   
    otherwise, 
        traj = sav_golay([d.traj_x d.traj_y],tset.flength,tset.forder,0,tset.fc); % this is trajectory (x,y)
        vel = sav_golay([d.traj_x d.traj_y],tset.flength,tset.forder,1,tset.fc); % this is velocity (x,y)
        acc = sav_golay([d.traj_x d.traj_y],tset.flength,tset.forder,2,tset.fc); % this is acceleration (x,y)
        jerk = sav_golay([d.traj_x d.traj_y],tset.flength,tset.forder,3,tset.fc); % this is jerk (x,y)
        nhands = 1;
end

% Need to identify single 'runs' (typically, from target onset to movement end)
% Assume (default situation) that there is a 'run' variable which is <>0 during
% the part of the recording which is relevant to the analysis, '0' otherwise
% In particular, runs is 1 in forward movements, and -1 in backward
% movements

% In legacy data, this information has to be inferred from other fields
% (typically, the 'state') as specified in the function 'get_runs')

if ~isfield(d,'runs')
    get_runs = varargin{1};
    d.runs = feval(get_runs, d);
else
    d.runs = d.runs';
    %disp('passato correttamente ->leggo runs');
    d.runs = (d.runs==1) - (d.runs==-1);
end

% switch exptype
% case 1, % MSROBOT
%     runs = notused_15'>0;  % 1 when center-out, -1 when out-center, 0 otherwise
%     runs = runs - (1-runs);
% case {0,3,4}, % EEGROBOT
%     runs = state'>1;  % 1 when center-out, 0 otherwise
% case 2, % Nabeel
%     runs = state'>4 & state' <=8; % this is used in Nabeel data
%     runs = state'>3 & state' <=6; % this is used in Nabeel data
% end

fwd_runs = d.runs>0;
bwd_runs = d.runs<0;

fwd_start = find([0 diff(fwd_runs)>0]);
fwd_end = find([0 diff(fwd_runs)<0]);

bwd_start = find([0 diff(bwd_runs)<0]);
bwd_end = find([0 diff(bwd_runs)>0]);

% In some data, runs is 1 at the very beginning. In these cases, need
% to add 1
if fwd_runs(1)
    fwd_start = [1 fwd_start];
end

if bwd_runs(1)
    bwd_start = [1 bwd_start];
end


% In some data, runs stays 1 at the very end. In these cases, need
% to add 1 at  the end
if fwd_runs(end)
    fwd_end = [fwd_end length(fwd_runs)];
end


if bwd_runs(end)
    bwd_end = [bwd_end length(bwd_runs)];
end

switch tset.expmode  
    case 'fwdbwd' % both forward and backward movements
        beg_traj = sort([fwd_start bwd_start]);
	    end_traj = sort([fwd_end bwd_end]);
    case 'fwd' % only forward (center-out) movements
	    beg_traj = fwd_start;
	    end_traj = fwd_end;    
   case 'bwd' % only backward (out-center) movements
	    beg_traj = bwd_start;
	    end_traj = bwd_end;
end


tset.Ntrials = length(beg_traj);

 if((size(end_traj+beg_traj))~=0)
if (end_traj(end)-beg_traj(end))==1
      tset.Ntrials = length(beg_traj)-1;
end
 end
% if there is a 'target' field in our data....
% keyboard;
if isfield(d,'targetpos_x') && isfield(d,'targetpos_y')  
    %disp('pass corr leggo targetpos_xy');
    xyT = [d.targetpos_x d.targetpos_y];
    xyT_exists = 1;
else
    xyT_exists=0;
end

for trial = 1:tset.Ntrials
 %   trial;
    onetime = time(beg_traj(trial):end_traj(trial))-time(beg_traj(trial));    
    onetraj = traj(beg_traj(trial):end_traj(trial),:);
    for i=1:size(traj,2) % 2-dim for unimanual, 4-dim for bimanual, 4-dim for dyad
      onetraj(:,i) =onetraj(:,i)-onetraj(1,i);
    end

    onevel = vel(beg_traj(trial):end_traj(trial),:);
    oneacc = acc(beg_traj(trial):end_traj(trial),:);
    onejerk = jerk(beg_traj(trial):end_traj(trial),:);
    
    switch task.type
        case 'dyad',  % this is a 6D vector
            oneforce = force(beg_traj(trial):end_traj(trial),:);
        otherwise,  % this is a scalar
            oneforce = d.force(beg_traj(trial):end_traj(trial));
    end
    
    onestate = d.state(beg_traj(trial):end_traj(trial)); %il onestate (1,1) nn è intero in ts=1
    oneruns = d.runs(beg_traj(trial):end_traj(trial)); %il oneruns (1,1) nn è intero in ts1

    if isfield(d,'rotation')
        onerot = d.rotation(beg_traj(trial):end_traj(trial));
        %disp('pass corr leggo rotationt');
    end
        %xyT_exists=0;
    if xyT_exists
        disp('xyT_Exists');
        onexyT = xyT(beg_traj(trial):end_traj(trial),:); %prende solo il target negativo -0.12 0 ? in ts=1
     
        if oneruns(1)
          angle = mean(atan2(onexyT(:,2),onexyT(:,1))); % perchè angolo viene pigreco?!
         % dist = (xyT(end_traj(trial),1)-tset.targets(:,1)).^2+(xyT(end_traj(trial),2)-tset.targets(:,2)).^2;

          dist = (xyT(beg_traj(trial),1)-tset.targets(:,1)).^2+(xyT(beg_traj(trial),2)-tset.targets(:,2)).^2; 
        else
          prevxyT = xyT(beg_traj(trial)-1,:);
           angle = mean(atan2(onexyT(:,2)-prevxyT(2),onexyT(:,1)-prevxyT(1)));    
          dist = (xyT(end_traj(trial),1)-prevxyT(1)-tset.targets(:,1)).^2+(xyT(end_traj(trial),2)-prevxyT(2)-tset.targets(:,2)).^2;
        end
                      

        if(angle<0) 
          angle=angle+2*pi;
        end
        theangle(trial) = angle;
        
        %disp('incredibilmente ho visto la posizione del target!');

        [m,tset.targetno(trial)]=min(dist);  %m è 2,38 .. sicuramente errato..
%keyboard;
    else % I don't have target position so lets try to guess from trajectory...
        
        %disp('non ho visto la posizione del target!');
        dist = (traj(end_traj(trial),1)-traj(beg_traj(trial),1)-tset.targets(:,1)).^2+...
               (traj(end_traj(trial),2)-traj(beg_traj(trial),2)-tset.targets(:,2)).^2;
        [m,tset.targetno(trial)]=min(dist);

        angle = target_angles(tset.targetno(trial));         
        onexyT = ones(length(onetime),1)*tset.targets(tset.targetno(trial),:);
    end
  
   % NB: in visuomotor rotation trials there is a 'rotation' field
   if isfield(d,'rotation')
       
       tset.traj{trial} = trajectory(nhands,onetime,onetraj,onevel,oneacc,onejerk,onexyT,onestate,oneforce,angle,onerot,tset.viapoints);   
   else
       tset.traj{trial} = trajectory(nhands,onetime,onetraj,onevel,oneacc,onejerk,onexyT,onestate,oneforce,angle,0.*onestate,tset.viapoints);
   end
   
   
  
end
tset.Ntrials
%plot_tset(time,d.state,traj,xyT,d.runs,beg_traj, end_traj)
% keyboard
tset;
% clear classes;
tset=class(tset,'targetset');
