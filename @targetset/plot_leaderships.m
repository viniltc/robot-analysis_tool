function [pow_matrix1, pow_matrix2,tc_matrix1,tc_matrix2, perc_trial, perc_vp, mean_vp]=plot_leaderships(tset,phaseno,pow_matrix1,pow_matrix2,tc_matrix1,tc_matrix2, perc_trial, perc_vp, mean_vp)
bwfigs = 0;
if bwfigs
    col_avg='k';
    col_fig='k';
    col_tgt='k';
    col_for = 'k';
    col_nul = 'k';
else
    col_avg='b';
    col_fig='r';
    col_tgt='g';
    col_for = 'm';
    col_nul = 'k';
end

% tlim = [-0.1 lim];
% slim = [-2 2];
munit = 'w';
scalebar.center=[1 0.05];
scalebar.xdim=0.20;
scalebar.ydim=0.10;
spdstrval=['Leadership - ', num2str(phaseno)];

% keyboard

th = atan2(tset.targets(:,2),tset.targets(:,1));

for trial = 1:tset.Ntrials
    for i = 1:length(th)
%                axes(hva(i));
               trial_set = tset.Ntrials*(phaseno-1);
               if mean(isforce(tset.traj{trial}))>0.5 % force is on
                   [pow1, pow2, tc11, tc22, tc12, tc21] = plot_leadership(tset.traj{trial},col_for);
               else
                   [pow1, pow2, tc11, tc22,tc12, tc21] = plot_leadership(tset.traj{trial},col_nul);
               end
               
%          keyboard
                 pow_matrix1(1:length(pow1),trial+trial_set) = pow1;
                 tc_matrix1(1,trial+trial_set) = tc11;
                 tc_matrix1(2,trial+trial_set) = tc12;
                 i11 = int32(tc11*100);
                 i12 = int32(tc12*100);
                 n_sample = 1;
%                   n_sample = 5;
                 
                 pow_matrix2(1:length(pow2),trial+trial_set) = pow2;
                 tc_matrix2(1,trial+trial_set) = tc21;
                 tc_matrix2(2,trial+trial_set) = tc22;
                 i21 = int32(tc21*100);
                 i22 = int32(tc22*100);
%         
%                  if (i11 < 3)
%                      disp('here here')
% keyboard
%                  end
%           keyboard
                 if((i11==0||i11==1||i11==2)&&(i22~=0||i22~=1||i22~=2))
                     pow1_tcvp = [zeros(1,n_sample*2+1); (pow1(i22-n_sample:i22+n_sample))']; %11 campioni (= 100ms) centrati su tc di vp1 e vp2
                     pow2_tcvp = [zeros(1,n_sample*2+1); (pow2(i22-n_sample:i22+n_sample))'];
                 elseif((i22==0||i22==1||i22==2)&&(i11~=0||i11~=1||i11~=2))
                     pow1_tcvp = [(pow1(i11-n_sample:i11+n_sample))'; zeros(1,n_sample*2+1)]; 
                     pow2_tcvp = [(pow2(i11-n_sample:i11+n_sample))'; zeros(1,n_sample*2+1)];
                 elseif((i22==0||i22==1||i22==2)&&(i11==0||i11==1||i11==2))
                     pow1_tcvp = [zeros(1,n_sample*2+1); zeros(1,n_sample*2+1)]; 
                     pow2_tcvp = [zeros(1,n_sample*2+1); zeros(1,n_sample*2+1)];
                 else
                 pow1_tcvp = [(pow1(i11-n_sample:i11+n_sample))'; (pow1(i22-n_sample:i22+n_sample))'];
                 pow2_tcvp = [(pow2(i11-n_sample:i11+n_sample))'; (pow2(i22-n_sample:i22+n_sample))'];
                 
                 
%                  if((i11==0||i11==1||i11==2)&&(i22~=0||i22~=1||i22~=2))
%                      pow1_tcvp = [zeros(1,n_sample*2); (pow1(i22-n_sample:i22))']; %11 campioni (= 100ms) centrati su tc di vp1 e vp2
%                      pow2_tcvp = [zeros(1,n_sample*2); (pow2(i22-n_sample:i22))'];
%                  elseif((i22==0||i22==1||i22==2)&&(i11~=0||i11~=1||i11~=2))
%                      pow1_tcvp = [(pow1(i11-n_sample:i11))'; zeros(1,n_sample*2)]; 
%                      pow2_tcvp = [(pow2(i11-n_sample:i11))'; zeros(1,n_sample*2)];
%                  elseif((i22==0||i22==1||i22==2)&&(i11==0||i11==1||i11==2))
%                      pow1_tcvp = [zeros(1,n_sample*2); zeros(1,n_sample*2)]; 
%                      pow2_tcvp = [zeros(1,n_sample*2); zeros(1,n_sample*2)];
%                  else
%                  pow1_tcvp = [(pow1(i11-n_sample:i11))'; (pow1(i22-n_sample:i22))'];
%                  pow2_tcvp = [(pow2(i11-n_sample:i11))'; (pow2(i22-n_sample:i22))'];
%                  
%                  end
%                  
                 % leadership percentages (normalizing with trial duration)
                   perc_trial(1,trial+trial_set) = (length(find(pow1<0))/length(pow1)).*100; %leader 1 
                   %percentages(3,trial+trial_set) = (length(find(pow1>=0))/interval1(1,trial+trial_set)).*100; %follower 1
                   perc_trial(2,trial+trial_set) = (length(find(pow2<0))/length(pow2)).*100; %leader 2
                   %percentages(4,trial+trial_set) = (length(find(pow2>=0))/interval2(1,trial+trial_set)).*100; %follower 2
                 
                 % leadership percentages (normalizing with trial duration)  
                   perc_vp(1,trial+trial_set) = (length(find(pow1_tcvp(1,:)<0))/length(pow1_tcvp(1,:))).*100; %perc1 vp1
                   perc_vp(3,trial+trial_set) = (length(find(pow1_tcvp(2,:)<0))/length(pow1_tcvp(2,:))).*100; %perc1 vp2
                   perc_vp(2,trial+trial_set) = (length(find(pow2_tcvp(1,:)<0))/length(pow2_tcvp(1,:))).*100; %perc2 vp1
                   perc_vp(4,trial+trial_set) = (length(find(pow2_tcvp(2,:)<0))/length(pow2_tcvp(2,:))).*100; %perc2 vp2  
                   
                   mean_vp(1,trial+trial_set) = mean(pow1_tcvp(1,:)); %perc1 vp1
                   mean_vp(3,trial+trial_set) = mean(pow1_tcvp(2,:)); %perc1 vp2
                   mean_vp(2,trial+trial_set) = mean(pow2_tcvp(1,:)); %perc2 vp1
                   mean_vp(4,trial+trial_set) = mean(pow2_tcvp(2,:)); %perc2 vp2  
                   
    end
end

end

%% old 
% function [pow_matrix1, pow_matrix2,tc_matrix1,tc_matrix2]=plot_leaderships(tset,phaseno,pow_matrix1,pow_matrix2,tc_matrix1,tc_matrix2)
% bwfigs = 0;
% if bwfigs
%     col_avg='k';
%     col_fig='k';
%     col_tgt='k';
%     col_for = 'k';
%     col_nul = 'k';
% else
%     col_avg='b';
%     col_fig='r';
%     col_tgt='g';
%     col_for = 'm';
%     col_nul = 'k';
% end
% 
% % tlim = [-0.1 lim];
% % slim = [-2 2];
% munit = 'w';
% scalebar.center=[1 0.05];
% scalebar.xdim=0.20;
% scalebar.ydim=0.10;
% spdstrval=['Leadership - ', num2str(phaseno)];
% 
% th = atan2(tset.targets(:,2),tset.targets(:,1));
% 
% for trial = 1:tset.Ntrials
%     for i = 1:length(th)
% %                axes(hva(i));
%                trial_set = tset.Ntrials*(phaseno-1);
%                if mean(isforce(tset.traj{trial}))>0.5 % force is on
%                    [pow1, pow2, tc11, tc22, tc12, tc21] = plot_leadership(tset.traj{trial},col_for);
%                else
%                    [pow1, pow2, tc11, tc22, tc12, tc21] = plot_leadership(tset.traj{trial},col_nul);
%                end
%                
% %                for ha = 1:2
% %                    [tc11, tc12] = get_vp_distances(tset.traj{trial},2); % lato slave
% %                    [tc21, tc22] = get_vp_distances(tset.traj{trial},1);
% %                end
%                
%                  pow_matrix1(1:length(pow1),trial+trial_set) = pow1;
% %                  tc_matrix1(1,trial+trial_set) = tc22;
% %                  tc_matrix1(2,trial+trial_set) = tc21;
%                  tc_matrix1(1,trial+trial_set) = tc11;
%                  tc_matrix1(2,trial+trial_set) = tc12;
%  
%                  pow_matrix2(1:length(pow2),trial+trial_set) = pow2;
% %                  tc_matrix2(1,trial+trial_set) = tc11;
% %                  tc_matrix2(2,trial+trial_set) = tc12;
%                  tc_matrix2(1,trial+trial_set) = tc22;
%                  tc_matrix2(2,trial+trial_set) = tc21;
% 
%     end
% end
% end

%%% from here....
