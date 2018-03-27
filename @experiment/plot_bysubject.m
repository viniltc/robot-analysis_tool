function plot_bysubject(exp, plotlist,figdir,epsdir,varargin)
% plots movements per subject and saves the figures
% plotlist: list of graph types
% figdir: directory where figs are saved
% epsdir: directory where eps figs are saved
% (C) V. Sanguineti (2008)

filename1 = sprintf('Lead_index_111.txt');
filename2 = sprintf('Lead_index_222.txt');
file1 = fopen(filename1,'wt');
file2 = fopen(filename2,'wt');
% filename1 = sprintf('powmatrix1.txt');
% filename2 = sprintf('powmatrix2.txt');
% file1 = fopen(filename1,'wt');
% file2 = fopen(filename2,'wt');


%  keyboard
for s  = 1:length(exp.subj)
    plot(exp.subj(s),plotlist,figdir,epsdir,varargin{1}); 
    %
    [lead11, lead12, lead21, lead22] =  plot(exp.subj(s),plotlist,figdir,epsdir,varargin{1}); 
    fprintf(file1,'%f\t%f\t%f\t%f\t%f\t%f\t\n',lead11, lead12);
    fprintf(file2,'%f\t%f\t%f\t%f\t%f\t%f\t\n',lead21, lead22);

%     [lead11, lead12, lead21, lead22, cc1, cc2] =  plot(exp.subj(s),plotlist,figdir,epsdir,varargin{1}); 
%     A = cc1  % Write this to file.
% 
% for ii = 1:size(A,1)
%     fprintf(file1,'%g\t',A(ii,:));
%     fprintf(file1,'\n');
% end

%      pause
end

fclose(file1);
fclose(file2);
