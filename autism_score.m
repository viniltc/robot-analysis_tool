

x=xlsread('C:\Users\Asus\Dropbox\pHHI sfn paper\experiment\asd_score.xlsx','Sheet1');

Chiara = x(:,1);
Alexis = x(:,2);
Francesco = x(:,3);
Sandeep = x(:,4);
Abdul = x(:,5);
Paolo = x(:,6);
Margerita = x(:,7);
Fabio = x(:,8);
Ninad = x(:,9);
Zeeshan = x(:,10);
Keerthi = x(:,11);
Pradeep = x(:,12);
Sand_kes = x(:,13);
Praveen = x(:,14);
Oussama = x(:,15);
Vishal = x(:,16);
Atal = x(:,17);
Jenny = x(:,18);
Gautham = x(:,19);
Rivo = x(:,20);
Luigi = x(:,21);
Jane = x(:,22);

asd_mat = [Francesco,Ninad,Zeeshan,Alexis,Paolo,Chiara,Margerita,Fabio,Pradeep,Abdul,Keerthi,Sandeep,Luigi,Jane,Praveen,Oussama,Vishal,Atal,Jenny,Sand_kes,Gautham,Rivo]
% 
 fid= fopen('asd_score_arr.txt','w');
 fprintf(fid,'%2d  %2d  %2d  %2d  %2d  %2d\n',asd_mat);

% figure
% %fig=[Abdul,Alexis,Atal,Chiara,Fabio,Francesco,Gautham,Jenny,Keerthi,Margerita,Ninad,Oussama,Paolo,Pradeep,Praveen,Sand_kes,Sandeep,Vishal,Zeeshan,Rivo,Luigi,Jane]
% fig=[Francesso,Ninad,Zeeshan,Alexis,Paolo,Chiara,Margerita,Fabio,Pradeep,Abdul,Keerthi,Sandeep,Luigi,Jane,Praveen,Oussama,Vishal,Atal,Jenny,Sand_kes,Gautham,Rivo]
% stem(fig')
% fig = gca;
% fig.XTick = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22];
% fig.XTickLabels = {'Abdul','Alexis','Atal','Chiara','Fabio','Francesco','Gautham','Jenny','Keerthi','Margerita','Ninad','Oussama','Paolo','Pradeep','Praveen','Sand kes','Sandeep','Vishal','Zeeshan', 'Rivo', 'Luigi','Jane'};
% fig.XTickLabelRotation = 90;
% grid on
% ylabel('Autism score')
% legend({'Social Skill','Attention to switching','Attention to details','Imagination','Communication','Overall'})
% legend ('boxoff')

load ('asd_score_arr.txt');
asd_score_arr
