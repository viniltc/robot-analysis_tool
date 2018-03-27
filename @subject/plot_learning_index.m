function plot_learning_index(subj,indname,inddesc)
% plots learning index for that indicator and that subject
% (C) V. Sanguineti (2008)

[indmat,indmat_catch,inds,indmat_catch_nofor,indmat_catch_norot] = get_indmatrix(subj,indname);



figure
set(gcf, 'pos', [100 100 350 300])
%set(gca,'ylim',rl)
set(gcf,'name',inddesc)

tsno = 1;
bnd = [];
for p = 1:length(subj.phases)
 indm = [];
 indm_c = [];
 indm_nr = [];
 indm_nf = [];
 for r = 1:length(subj.order{p})
    indm = [indm mean(median((indmat{tsno+(r-1)}')))];
    indm_c = [indm_c mean(median((indmat_catch{tsno+(r-1)}')))];
    indm_nr = [indm_nr mean(median((indmat_catch_norot{tsno+(r-1)}')))];
    indm_nf = [indm_nf mean(median((indmat_catch_nofor{tsno+(r-1)}')))];

      
 end
 line(tsno+(1:length(subj.order{p})), -indm_c./abs(indm-indm_c),'col','b','LineWidth',2 )
 line(tsno+(1:length(subj.order{p})), -indm_nr./abs(indm-indm_nr),'col','g','LineWidth',2 )
 line(tsno+(1:length(subj.order{p})), -indm_c./abs(indm-indm_nf),'col','r','LineWidth',2 )


 tsno = tsno+length(subj.order{p});
 bnd(p) = tsno;
end  
yrange = get(gca,'ylim');
for p = 1:(length(subj.phases)-1)
   line(bnd(p)*[1 1]+0.5,yrange, 'col','k','lines','-');
end


xlabel('MOVEMENT SET')
ylabel('LEARNING INDEX')
ylabel(inddesc)
set(gca, 'box', 'off')
%axis square
set(gcf,'name',['LI '])
