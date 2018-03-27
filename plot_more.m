function plot_more
% plots additional things on traj plots

thr = 0:pi/20:2*pi;
vpradius = 2.5/1000;
vp1 = [30 -30]/1000;
vp2 = [70  30]/1000; 
patch(vp1(1)+vpradius*cos(thr),vp1(2)+vpradius*sin(thr),'b','edgecol','k');
patch(vp2(1)+vpradius*cos(thr),vp2(2)+vpradius*sin(thr),'r','edgecol','k');

end

