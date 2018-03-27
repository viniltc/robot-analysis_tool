function inds = get_indicator(tset,indname)

% keyboard

inds = zeros(tset.Ntrials,1);
for trial = 1:tset.Ntrials
trial

switch(indname)
    % dyad indicators:
    case 'ts1',  [inds(trial,1),ts2,te1,te2]=get_dyad_durations(tset.traj{trial}); % partner 1: time of start 
    case 'ts2',  [ts1,inds(trial,1),te1,te2]=get_dyad_durations(tset.traj{trial}); % partner 2: time of start 
    case 'te1',  [ts1,ts2,inds(trial,1),te2]=get_dyad_durations(tset.traj{trial}); % partner 1: time of end
    case 'te2',  [ts1,ts2,te1,inds(trial,1)]=get_dyad_durations(tset.traj{trial}); % partner 2: time of end
    case 'md12', [tc1,tc2,md1, inds(trial,1)]=get_vp_distances(tset.traj{trial},1); % partner 1: minimum distance to VP1 
    case 'md11', [tc1,tc2,inds(trial,1),md2]=get_vp_distances(tset.traj{trial},1); % partner 1: minimum distance to VP2
    case 'md22', [tc1,tc2,md1,inds(trial,1)]=get_vp_distances(tset.traj{trial},2); % partner 2: minimum distance to VP1
    case 'md21', [tc1,tc2,inds(trial,1),md2]=get_vp_distances(tset.traj{trial},2); % partner 2: minimum distance to VP2
    case 'tc12', [tc1,inds(trial,1),md1,md2]=get_vp_distances(tset.traj{trial},1); % partner 1: crossing time to VP1
    case 'tc11', [inds(trial,1),tc2,md1,md2]=get_vp_distances(tset.traj{trial},1); % partner 1: crossing time to VP2
    case 'tc22', [tc1, inds(trial,1),md1,md2]=get_vp_distances(tset.traj{trial},2); % partner 2: crossing time to VP1
    case 'tc21', [inds(trial,1),tc2,md1,md2]=get_vp_distances(tset.traj{trial},2); % partner 2: crossing time to VP2
    case 'eff1', [inds(trial,1),eff2]=get_effort(tset.traj{trial}); % partner 1: effort
    case 'eff2', [eff1,inds(trial,1)]=get_effort(tset.traj{trial}); % partner 2: effort
    case 'pow1', [inds(trial,1),pow2]=get_power(tset.traj{trial}); % partner 1: power
    case 'pow2', [pow1,inds(trial,1)]=get_power(tset.traj{trial}); % partner 2: power
    case 'pow_vp1', [inds(trial,1),pow_vp2]=get_power_vp(tset.traj{trial}); % partner 1: power
    case 'pow_vp2', [pow_vp1,inds(trial,1)]=get_power_vp(tset.traj{trial}); % partner 2: power
%     case 'lead1', [inds(trial,1), lead2]=get_power_fac(tset.traj{trial});
%     case 'lead2', [lead1, inds(trial,1)]=get_power_fac(tset.traj{trial});
    case 'speed1', [inds(trial,1),speed2,R]=get_speed_correlation(tset.traj{trial}); 
    case 'speed2', [speed1, inds(trial,1),R]=get_speed_correlation(tset.traj{trial}); 
    case 'R', [speed1,speed2, inds(trial,1)]=get_speed_correlation(tset.traj{trial}); 
    case 'forcenorm', [inds(trial,1)]=get_intforce(tset.traj{trial}); 
    case 'score1', [inds(trial,1),score2]=get_scores(tset.traj{trial}); 
    case 'score2', [score1,inds(trial,1)]=get_scores(tset.traj{trial}); 
    case 'pkspeed', [inds(trial,1),asp]=get_speed(tset.traj{trial}); 
    case 'avspeed', [psp,inds(trial,1)]=get_speed(tset.traj{trial}); 
    
     case 'LI11', [inds(trial,1),LI12,  LI21, LI22]=get_power_fac(tset.traj{trial}); 
     case 'LI12', [LI11, inds(trial,1), LI21, LI22]=get_power_fac(tset.traj{trial}); 
     case 'LI21', [LI11,LI12,inds(trial,1), LI22]=get_power_fac(tset.traj{trial}); 
     case 'LI22', [LI11,LI12,LI21, inds(trial,1)]=get_power_fac(tset.traj{trial}); 
        
        
    % bimanual indicators:
    case 'hma', [inds(trial,1), mean_deviation]= get_delta_horiz_deviation(tset.traj{trial});
    case 'hme', [max_deviation, inds(trial,1) ]= get_delta_horiz_deviation(tset.traj{trial});
    case 'lme',  [inds(trial,1)] = get_delta_longit_displacement(tset.traj{trial});
    case 'dft', inds(trial,1) = get_final_deltat(tset.traj{trial});
    
    % unimanual indicators:
    case 'dur', [inds(trial,1), adur,ddur,rtime]= get_duration(tset.traj{trial});
    case 'adu', [tdur,inds(trial,1), ddur,rtime]= get_duration(tset.traj{trial});        
    case 'ddu', [tdur,adur, inds(trial,1),rtime]= get_duration(tset.traj{trial});
    case 'sim', [tdur, adur,ddur,rtime]= get_duration(tset.traj{trial});
                inds(trial,1) = adur/ddur;    
    case 'rti', [tdur,adur, ddur, inds(trial,1)]= get_duration(tset.traj{trial});
    case 'ae1', inds(trial,1)=get_aiming_error(tset.traj{trial},0.1,'regr');    
    case 'aa1', inds(trial,1)=abs(get_aiming_error(tset.traj{trial},0.1,'regr'));    
    case 'ae3', inds(trial,1)=get_aiming_error(tset.traj{trial},0.3,'regr');    
    case 'aa3', inds(trial,1)=abs(get_aiming_error(tset.traj{trial},0.3,'regr'));    
    case 'ae2', inds(trial,1)=get_aiming_error(tset.traj{trial},0.2,'regr');
    case 'aef', inds(trial,1)=get_aiming_error(tset.traj{trial},0.1,'final');
    case 'aee', inds(trial,1)=get_aiming_error(tset.traj{trial},'all','final');
    case 'aep', inds(trial,1)=get_aiming_error(tset.traj{trial},'peak','final');
    case 'aap', inds(trial,1)=abs(get_aiming_error(tset.traj{trial},'peak','final'));
    case 'pat', [inds(trial,1),ampl] = get_length(tset.traj{trial});
    case 'len', [pathlength, inds(trial,1)] = get_length(tset.traj{trial});
    case 'lin', [pathlength, amplitude] = get_length(tset.traj{trial});
                inds(trial,1)=100.*(pathlength/amplitude-1);
    case 'pax', [inds(trial,1),peak_accy,tmp]=get_peak_acceleration(tset.traj{trial}); 
    case 'pay', [peak_accx,inds(trial,1),tmp]=get_peak_acceleration(tset.traj{trial}); 
    case 'pae', [tmp1,tmp2,inds(trial,1)]=get_peak_acceleration(tset.traj{trial});     
    case 'jer', [inds(trial,1),logprejerk,logpostjerk]=get_jerk(tset.traj{trial});
    case 'jei', [logtotjerk,inds(trial,1),logpostjerk]=get_jerk(tset.traj{trial});
    case 'jef', [logtotjerk,logprejerk,inds(trial,1)]=get_jerk(tset.traj{trial});
    case 'jra', [teutotjerk,teuprejerk,teupostjerk]=get_jerk(tset.traj{trial});
                %inds(trial,1) = 10.^(logpostjerk-logprejerk);
                inds(trial,1) = teupostjerk./teuprejerk;
    case 'lat', inds(trial,1) = get_lateral_deviation(tset.traj{trial});
    case 'cai',[inds(trial,1),late_cac]=get_centripetal_acceleration(tset.traj{trial});
    case 'caf',[early_cac, inds(trial,1)]=get_centripetal_acceleration(tset.traj{trial});
    case 'fer',inds(trial,1)=get_final_error(tset.traj{trial});
    case 'psp',[inds(trial,1),asp]=get_speed(tset.traj{trial});
    case 'asp',[psp, inds(trial,1)]=get_speed(tset.traj{trial});
        
        
end

end