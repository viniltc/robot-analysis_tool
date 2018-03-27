function traj = trajectory(nhands,time,pos,vel,acc,jerk,xyT,state,force,angle,rotation,viapoints)
% Creates a TRAJECTORY object
% Takes the following vectors as inputs:
% time (Tx1): time instants
% pos, vel, acc, jerk (Tx(nhands*trsize) ): trajectory and its time derivatives
% xyT (Tx2): target position in time
% state (Tx(nhands*trsize)): 'state' of FSM definining the protocol
% force (Tx1): force (1=on, 0=off) as set by the robot or (in
% dyad experiments) interaction force vector (Tx(nhands*trsize))
% angle (1x1): target direction 
% rotation (Tx1) : rotation angle
% (C) V. Sanguineti (2008, 2015)

traj.nhands = nhands; % one if unimanual, 2 if bimanual or dyad

traj.trsize = size(pos,2)/nhands;  % 2 in planar movements, 3 in 3D movements

traj.time = time;

% The latter may be struct if bimanual experiment....
traj.pos  = pos;
traj.vel = vel;
traj.acc = acc;
traj.jerk = jerk;

traj.xyT = xyT;
traj.state = state; 
traj.intforce = force; % this is interaction force (a 2D, 4D or 6D variable)
traj.force = sqrt(sum(force.^2,2))>0.5; % this is presence of force (a boolean variable)
traj.angle = angle;

traj.rotation = rotation;
traj.viapoints = viapoints;


traj.fc = 1./mean(diff(time));
traj.thtype = 'abs';
traj.tol = 0.02;
%traj.thtype = 'rel';
%traj.tol = 0.1;

for hand = 1:traj.nhands
  if traj.trsize == 2 % planar movements
    speed = sqrt(traj.vel(:,2*(hand-1)+1).^2+traj.vel(:,2*(hand-1)+2).^2);
    distance = sqrt((traj.pos(:,2*(hand-1)+1)-traj.xyT(end,1)).^2+...
                   (traj.pos(:,2*(hand-1)+2)-traj.xyT(end,2)).^2);
  else
    speed = sqrt(traj.vel(:,3*(hand-1)+1).^2+traj.vel(:,3*(hand-1)+2).^2+traj.vel(:,3*(hand-1)+3).^2);
    distance = sqrt((traj.pos(:,3*(hand-1)+1)-traj.xyT(end,1)).^2+...
                   (traj.pos(:,3*(hand-1)+2)-traj.xyT(end,2)).^2+...
                   (traj.pos(:,3*(hand-1)+3)-traj.xyT(end,3)).^2);
  end
  
  [traj.pkspeed(hand),traj.pkind(hand)]=max(speed); %determina il valore massimo e l'indice relativo
  distance = distance./max(distance);
  i2 = find(distance < 0.3);
  i1 = find(distance > 0.7);
   
  if isempty(i2)
       traj.t2(hand) = length(distance);
  else
       traj.t2(hand) = min(i2);
  end
  
  if isempty(i1)
       traj.t1(hand) = 1;
  else
       traj.t1(hand) = max(i1);
  end
%   
  
  
  % Calculate the reaction time (onset of movement with respect to onset of
  % target
  
  %traj.pkind=min(findpeaks(speed,5,2*ttol)); %this find first peak of speed... Useful if subj moves back
  %traj.pkspeed = speed(traj.pkind);
  
  switch(traj.thtype)
  case 'abs',
   ttol = traj.tol;      %tol è la tolleranza: i suoi valori sono scelti in base a delle regole precise
  case 'rel',
   ttol = (traj.tol)*traj.pkspeed(hand);
  end  
  
  % Uses a second threshold (2*ttol) to restrict the search of movement start.
  % In this way, small bumps (below 2*ttol) at movement start are ignored (false starts), but greater bumps are not 
  aind = min(find(speed >2*ttol));        % aind is an overestimation of the true reaction time   
  rind = max(find(speed(1:aind) < ttol)); % then look below aind to determine true reaction time...
  if isempty(rind)
   rind = 1;
  end
  traj.rind(hand) = rind;
  
  traj.rtime(hand) = traj.time(rind)-traj.time(1);        % reaction  time
  traj.adur(hand)=traj.time(traj.pkind(hand))-traj.time(rind);  % acceleration time
  
  % We then calculate the time at which the movement ends...
  % Similar to reaction time, we want to ignore small bumps (final adjustments) 
  % but not large ones which may reflect sub-movements
  eind = max(find(speed>2*ttol));  % this underestimates the true end time...
  %sind=  traj.pkind(hand)-1+min(find(speed(traj.pkind(hand):end)<ttol));  % this fails if there are multiple peaks (eg bumps)
  sind=  eind-1+min(find(speed(eind:end)<ttol));  % this fails if there are multiple peaks (eg bumps)
  if isempty(sind)
    sind = length(speed);
  end
  traj.lind(hand) = sind;
  
  % intmax is lind plus 200 ms added at the end
  traj.intmax(hand) = min([traj.lind(hand)+round(0.2*traj.fc),size(traj.pos,1)]);
  

  traj.interval{hand} = (traj.rind(hand)):traj.intmax(hand); %traj.lind;

  % interval is based on lind...
  traj.interval{hand} = (traj.rind(hand)):traj.lind(hand); %traj.lind;

  traj.ddur(hand) = traj.time(sind)-traj.time(traj.pkind(hand)); % deceleration time
  traj.tdur(hand) =  traj.adur(hand) + traj.ddur(hand);          % total duration 
end
traj = class(traj,'trajectory');


