function plot_adaptation(subj,indname,inddesc)
% plots time series (absolute value) for that indicator and that subject
% (C) V. Sanguineti (2008)

[indmat,indmat_catch,inds, indmat_catch_nofor, indmat_catch_norot] = get_indmatrix(subj,indname);
rl=get_range(indname);

figure
set(gcf, 'pos', [100 100 350 300])
%set(gca,'ylim',rl)
set(gcf,'name',inddesc)

tsno = 1;
bnd = [];
for p = 1:length(subj.phases)
 indm = [];
 for r = 1:length(subj.order{p})
      indm = [indm mean(median(abs(indmat{tsno+(r-1)}')))];
      
 end
 line(tsno+(1:length(subj.order{p})), indm,'col','b','LineWidth',2 )
 tsno = tsno+length(subj.order{p});
bnd(p) = tsno;
end  
yrange = get(gca,'ylim');
for p = 1:(length(subj.phases)-1)
   line(bnd(p)*[1 1]+0.5,yrange, 'col','k','lines','-');
end


xlabel('MOVEMENT SET')
ylabel(inddesc)
set(gca, 'box', 'off')
%axis square
