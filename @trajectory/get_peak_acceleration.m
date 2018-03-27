function [peak_accx,peak_accy,acc_ang_error]=get_peak_acceleration(traj)
% Calculates peak acceleration (X and Y components, angular error)
% (C) V. Sanguineti, 2008

red_fpos= traj.pos(traj.rind:traj.pkind,:);
red_facc= traj.acc(traj.rind:traj.pkind,:);
nfacc = sqrt(red_facc(:,1).^2+red_facc(:,1).^2);

red_fvel=traj.vel(traj.rind:traj.pkind,:);
dir_vers = red_fvel./(sqrt(red_fvel(:,1).^2+red_fvel(:,2).^2)*ones(1,2));

tacc = red_facc(:,1).*dir_vers(:,1)+red_facc(:,2).*dir_vers(:,2); % vector quantity of the tangential accleration
pacc = red_facc - (tacc*ones(1,2)).*dir_vers; % Residual acceleration
cacc = -pacc(:,1).*dir_vers(:,2)+pacc(:,2).*dir_vers(:,1); %Centripetal Acceleration

% [peak_accx,tindx]=max(abs(red_facc(:,1)));
% [peak_accy,tindy]=max(abs(red_facc(:,2)));
% peak_accx = red_facc(tindx,1);
% peak_accy = red_facc(tindy,2);

% this definition is exactly like in Gordon et al (1994)
[peak_acc,tind]=max(tacc); % Find peak of tangential acceleration
peak_dir = atan2(red_fpos(tind,2)-red_fpos(1,2),red_fpos(tind,1)-red_fpos(1,1));
peak_accx = peak_acc*cos(peak_dir);
peak_accy = peak_acc*sin(peak_dir);

actual_dir = traj.angle; 
isneg = find(actual_dir<0);
actual_dir(isneg) = actual_dir(isneg)+2*pi;
actual_dir = rem(actual_dir,2*pi);

ae = (peak_dir-actual_dir);
ae = rem(ae,2*pi);
ae(find(ae > pi)) =  ae(find(ae > pi)) -2*pi;
ae(find(ae < -pi)) =  ae(find(ae < -pi))+2*pi;
acc_ang_error = 180./pi.*ae;
