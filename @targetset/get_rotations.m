function rotations = get_rotations(tset)

for trial = 1:tset.Ntrials    
    rotations(1,trial)=mean(isrotation(tset.traj{trial}))~=0;
end