function plot_peakaccs(tset,phaseno)
bwfigs = 0;
if bwfigs
    col_avg='k';
    col_fig='k';
    col_tgt='k';
else
    col_avg='b';
    col_fig='r';
    col_tgt='g';
end

tlim = [-6 6];
munit = 'm/s^2';
scalebar.center=[5 -5];
scalebar.xdim=1;
scalebar.ydim=1;

accstrval=sprintf('Peak acceleration - %d', phaseno);
figure
set(gcf,'Name',accstrval)
set(gcf,'pos',[0 0 150 150]);
axis('square');
set(gca,'xlim',tlim,'ylim',tlim);
    
xlabel(['a_x [', munit,']']);
ylabel(['a_y [', munit,']']);
    
%line(scalebar.center(1)+[0 scalebar.xdim],scalebar.center(2)+[0 0],'col','k');
%line(scalebar.center(1)+[0 0],scalebar.center(2)+[0 scalebar.ydim],'col','k');

line(tlim,[0,0],'col','k');
line([0,0],tlim,'col','k');
set(gca,'visible','off');

% th = (0:(tset.Ntargets-1))*(pi/tset.Ntargets);
% thr = 0:pi/20:2*pi;
% %for i = 1:tset.Ntargets
% %       line(tset.amplitude*cos(th(i))+tset.tgsize.*cos(thr),...
%            tset.amplitude*sin(th(i))+tset.tgsize.*sin(thr),'col',col_tgt);
%end

for trial = 1:tset.Ntrials
    [peak_accx,peak_accy]=get_peak_acceleration(tset.traj{trial});
    if isforce(tset.traj{trial})
     col = 'm';
    else
     col = 'k';
    end
    line(peak_accx,peak_accy,'col',col,'marker','.','lines','none');  
end     