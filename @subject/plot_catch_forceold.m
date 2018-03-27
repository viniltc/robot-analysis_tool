function plot_catch_force(subj,indname,inddesc)
% plots time series (both force and catch trials) for that indicator and that subject
% (C) V. Sanguineti (2008)

[indmat,indmat_catch] = get_indmatrix(subj,indname);
rl=get_range(indname);


figure
set(gcf, 'pos', [100 100 350 300])
%set(gca,'ylim',rl)
set(gcf,'name',inddesc)

tsno = 1;
bnd = [];
for p = 1:length(subj.phases)
 indm = [];
 indm_c = [];
 for r = 1:length(subj.order{p})
    indm = [indm mean(median((indmat{tsno+(r-1)}')))];
    indm_c = [indm_c mean(median((indmat_catch{tsno+(r-1)}')))];
      
 end
 line(tsno+(1:length(subj.order{p})), indm,'col','b','LineWidth',2 )
 line(tsno+(1:length(subj.order{p})), indm_c,'col','r','LineWidth',2 )
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
