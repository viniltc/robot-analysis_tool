function tset=targetset(data,data_format,targets,tgsize,fc)
% Constructor of a TARGETSET (block of movement trials) object
% One target set is characterized by: 
%   data: data matrix
%   data_format: data format
%   targets  (Ntg x 2)
%   tgsize: target size
%   fc: nominal sampling rate
% (C) V. Sanguineti (2008)

exptype = 1; % MS Asy data
exptype = 0; % EEG Babiloni data
exptype = 2; % Nabeel data
exptype = 3; % Bimanual data
exptype = 4; % Bimanual control data


tset.fc=fc;
tset.Ntargets = size(targets,1);
tset.targets = targets;
tset.tgsize = tgsize;

% filter params
tset.flength=27;   % 270 ms window, this is about 7.5 Hz
tset.forder=6;
%tset.flength=17;   % 170 ms window, this is about 11 Hz

target_angles = atan2(targets(:,2), targets(:,1))';

% Assuming first row in data is time...
rawtime = data(1,:)-data(1,1);
maxtime = max(rawtime);

% Resample at uniform sampling rate (in case it is not...)
newtime = 0:(1/fc):maxtime;
xy_target_f = zeros(size(data,1),length(newtime));
xy_target_f(1,:) = newtime; % assuming that first row is always time...
xy_target_f(2:end,:) = interp1(rawtime,data(2:end,:)',newtime)';

% Generates one tset 'field' per element in the data_format cell array
[ncomps,nsamples] = size(xy_target_f);
for comp=1:ncomps
    compname = data_format{comp};
    eval([compname ,' =xy_target_f(',num2str(comp),',:)'';']);
end
% just in case 'force' has not been set
if ~exist('force','var') & exist('force_x','var') & exist('force_y','var')
    force = (sqrt(force_x.^2+force_y.^2))>0;
end

time = time-time(1); % time starts from first sample in the target set
xyT = [targetpos_x targetpos_y];
traj = sav_golay([traj_x traj_y],tset.flength,tset.forder,0,tset.fc); % this is trajectory (x,y)
vel = sav_golay([traj_x traj_y],tset.flength,tset.forder,1,tset.fc); % this is velocity (x,y)
acc = sav_golay([traj_x traj_y],tset.flength,tset.forder,2,tset.fc); % this is acceleration (x,y)
jerk = sav_golay([traj_x traj_y],tset.flength,tset.forder,3,tset.fc); % this is jerk (x,y)


% Need to identify single 'runs' (typically, from target onset to movement end)
% Assume (default situation) that there is a 'run' variable which is '1' during
% part of the recording which is relevant to the analysis, '0' otherwise

% In legacy data, this information has to be inferred from other fields
% (typically, the 'state')

switch exptype
case 1, % MSROBOT
    center_out = notused_15'>0;  % 1 when center-out, 0 otherwise
    target_transition = [find([0 diff(center_out)]) length(center_out)];  % these are times of target transition
case {0,3,4}, % EEGROBOT
    center_out = state'>1;  % 1 when center-out, 0 otherwise
    target_transition = [find([0 diff(center_out)]) length(center_out)];  % these are times of target transition
case 2, % Nabeel
    center_out = state'>4 & state' <=8; % this is used in Nabeel data
    center_out = state'>3 & state' <=6; % this is used in Nabeel data
    target_transition = [find([0 diff(center_out)])];
end
    
% In some data, center_out is 1 at the very beginning. In these cases, need
% to add 1
if center_out(1)
    target_transition = [1 target_transition];
end

switch exptype  % forward and backward movements
    case 1 % forward and backward movements
        beg_traj = [target_transition(1:end)];
	    end_traj = [target_transition(2:end)-1 length(xy_target_f)];
	    tset.Ntrials = length(target_transition)-1;
    case 0 % only forward (center-out) movements
	    beg_traj = [target_transition(1:2:(end-1))];
	    end_traj = [target_transition(2:2:end)];
	    tset.Ntrials = length(beg_traj);
    case 2 % Nabeel
        beg_traj = [target_transition(1:2:end-1)];
	    end_traj = [target_transition(2:2:end)];
        
        if (end_traj(end)-beg_traj(end))==1
            tset.Ntrials = length(beg_traj)-1;
        else
	        tset.Ntrials = length(beg_traj);   
        end
    case 3 % only forward (center-out) movements
	    beg_traj = [target_transition(1:4:(end-1))];
	    end_traj = [target_transition(2:4:end)];
	    tset.Ntrials = length(beg_traj);

    case 4 %  forward and backward movements
	    beg_traj = [target_transition(1:2:(end-1))];
	    end_traj = [target_transition(2:2:end)];
	    tset.Ntrials = length(beg_traj);

end

for trial = 1:tset.Ntrials
    onetime = time(beg_traj(trial):end_traj(trial))-time(beg_traj(trial));
    onetraj = traj(beg_traj(trial):end_traj(trial),:);
    
    if exptype == 4 %bimanual control data consider forward and backward movts...
      onetraj(:,1) =onetraj(:,1)-onetraj(1,1);
      onetraj(:,2) =onetraj(:,2)-onetraj(1,2);
    end
    
    onevel = vel(beg_traj(trial):end_traj(trial),:);
    oneacc = acc(beg_traj(trial):end_traj(trial),:);
    onejerk = jerk(beg_traj(trial):end_traj(trial),:);
    onexyT = xyT(beg_traj(trial):end_traj(trial),:);
    
    onestate = state(beg_traj(trial):end_traj(trial));
    oneforce = force(beg_traj(trial):end_traj(trial));
    onecenter_out = center_out(beg_traj(trial):end_traj(trial));

    if onecenter_out(1)
      angle = mean(atan2(onexyT(:,2),onexyT(:,1)));
      dist = (xyT(end_traj(trial),1)-targets(:,1)).^2+(xyT(end_traj(trial),2)-targets(:,2)).^2;
    else
      prevxyT = xyT(beg_traj(trial)-1,:);
      angle = mean(atan2(onexyT(:,2)-prevxyT(2),onexyT(:,1)-prevxyT(1)));    
      dist = (xyT(end_traj(trial),1)-prevxyT(1)-targets(:,1)).^2+(xyT(end_traj(trial),2)-prevxyT(2)-targets(:,2)).^2;
    end
    if(angle<0) 
      angle=angle+2*pi;
    end
    theangle(trial) = angle;
    
    
    [m,tset.targetno(trial)]=min(dist); 
    
    
    tset.traj{trial} = trajectory(onetime,onetraj,onevel,oneacc,onejerk,onexyT,onestate,oneforce,angle);
    
    
end
tset.Ntrials
%plot_tset(time,state,traj,xyT,center_out,beg_traj, end_traj)
%keyboard
tset=class(tset,'targetset');
