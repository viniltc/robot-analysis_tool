function plot_tset(time,state,traj,xyT,center_out,beg_traj, end_traj)
figure
subplot(2,1,1)
plot(time,state)
hold on
plot(time,center_out,'r')
subplot(2,1,2)
plot(time,traj(:,1),'r',time,traj(:,2),'b')
hold on
plot(time,xyT(:,1),'r:',time,xyT(:,2),'b:')
stem(time(beg_traj),ones(1,length(beg_traj)),'c')
stem(time(end_traj),ones(1,length(end_traj)),'m')
