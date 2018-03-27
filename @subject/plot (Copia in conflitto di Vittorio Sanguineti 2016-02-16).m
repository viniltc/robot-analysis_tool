function plot(subj,displist,figdir,epsdir,varargin)
% plots SUBJECT data and saves the figures
% displ_list: list of graph types
% figdir: directory where figs are saved
% epsdir: directory where eps figs are saved
% (C) V. Sanguineti (2008)

% ntsets = length([subj.order{:}]);
tsno = 1;
for p = 1:length(subj.phases)
  for r = 1:length(subj.order{p})
   tset = subj.tset{tsno}; 
   for i = 1:length(displist)
        figname = sprintf('%s_%s_%d_%s',...
                         subj.name,subj.phases{p},subj.order{p}(r),displist{i})
        spcs = find(figname==' ');
        figname(spcs)=[];
        switch(displist{i})
           case 'traj',  plot_trajs(tset,tsno,varargin{1});
           case 'peakacc', plot_peakaccs(tset,tsno);
           case 'speed',plot_speeds(tset,tsno);
           case 'int_force',plot_intforce(tset,tsno);
           case 'jerk',plot_jerks(tset,tsno);    
           case '1traj', plot_1trajs(tset,tsno);
           case 'power',plot_powers(tset,tsno);
           case 'vpdist',plot_vpdistances(tset,tsno); % this is for dyads...
           case 'variability',plot_variability(tset,tsno); 
           case 'min_viadist',plot_vpdistances(tset,tsno);
        end
        title(subj.name);
       saveas(gcf,[figdir,figname],'fig');
       preprint;
       eval(['print -depsc ', epsdir,figname])
    end   % of i
    tsno = tsno+1;
   end % of r
 end % of p
