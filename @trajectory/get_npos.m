function [ntime,npos] = get_npos( traj,ns )
npos = cell(traj.nhands,1);
for ha=1:traj.nhands
    time = traj.time(traj.interval{ha})-traj.time(traj.interval{ha}(1));
    pos = traj.pos(traj.interval{ha},3*(ha-1)+(1:3));
    % normalize time
    normtime=time./max(time);
    ntime = linspace(0,1,ns);
    npos{ha} = interp1(normtime,pos,ntime);
end

