function forces = get_forces(tset)

for trial = 1:tset.Ntrials    
    %forces(1,trial)=any(isforce(tset.traj{trial}));
    forces(1,trial)=mean(isforce(tset.traj{trial}))>0.5;

end