function[score1,score2]=get_scores(traj,par)


 perf1 = 0.0025;
 perf2 = 0.02; 
 max_score = 100;
 
 K = 2*log(max_score*2-1)/(perf1-perf2); % slope of sigmoid 

% keyboard
 
% [md1,~]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,2)).^2+... %min dist da VP1 (x-y schermo)
%                    (traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,3)).^2));
% [md2,~]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,2)).^2+... %min dist da VP2
%                    (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,3)).^2));
%  md1= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,2))).^2);
%  md2= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,2))).^2);
%mean_x = min(sqrt(((traj.pos(traj.interval{1},3*(par-1)+2)) -...
%                (traj.pos(traj.interval{1},3*(par-1)+3))).^2));
%mean_y = min(sqrt(((traj.pos(traj.interval{2},3*(par-1)+2)) -...
%                (traj.pos(traj.interval{2},3*(par-1)+3))).^2));
par = 1;
[md1,~]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(1,2)).^2+...
                   (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(1,3)).^2));
[md2,~]= min(sqrt((traj.pos(traj.interval{par},3*(par-1)+2)-traj.viapoints(2,2)).^2+...
                   (traj.pos(traj.interval{par},3*(par-1)+3)-traj.viapoints(2,3)).^2));
               
 forcenorm = mean(sqrt(traj.intforce(:,1).^2+traj.intforce(:,2).^2));             

 score1_calculate = sqrt(md1)+0.5*sqrt(forcenorm); %
 score2_calculate = sqrt(md2)+0.5*sqrt(forcenorm); % 
 
 score1_calculate = sqrt(md1); %
 score2_calculate = sqrt(md2); % 
					
% keyboard
 score1 = max_score./(1+exp(K*(score1_calculate - (perf1+perf2)./2)));
 score2 = max_score./(1+exp(K*(score2_calculate - (perf1+perf2)./2)));

% score_1 = md2;
% score_2 = md1;

  
end
