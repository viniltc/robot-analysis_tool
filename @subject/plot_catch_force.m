function plot_catch_force(subj,indname,inddesc)
% plots time series (both force and catch trials) for that indicator and that subject
% (C) V. Sanguineti (2008)

[indmat,indmat_catch, inds, indmat_catch_nofor, indmat_catch_norot] = get_indmatrix(subj,indname);
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
 indm_nf = [];
 indm_nr = [];
 for r = 1:length(subj.order{p})
    indm = [indm mean(median((indmat{tsno+(r-1)}')))];
    indm_c = [indm_c mean(median((indmat_catch{tsno+(r-1)}')))];
    indm_nf = [indm_nf mean(median((indmat_catch_nofor{tsno+(r-1)}')))];
    indm_nr = [indm_nr mean(median((indmat_catch_norot{tsno+(r-1)}')))];


      
 end
if length(subj.order{p})==1
line(tsno+(1:length(subj.order{p})), indm,'col','k','LineWidth',2,'marker','*')
else
    line(tsno+(1:length(subj.order{p})), indm,'col','k','LineWidth',2 )
end
    
 line(tsno+(1:length(subj.order{p})), indm_c,'col','r','LineWidth',2 )
%  line(tsno+(1:length(subj.order{p})), indm_nf,'col','g','LineWidth',2 )
%   line(tsno+(1:length(subj.order{p})), indm_nr,'col','w','LineWidth',2 )
if p==3

line(tsno+(2:2:length(subj.order{p})),indm_nf(2:2:length(subj.order{p})),'col','r','linew',2,'marker','o')
line(tsno+(1:2:length(subj.order{p})),indm_nr(1:2:length(subj.order{p})),'col','m','linew',2,'marker','o')
end

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
