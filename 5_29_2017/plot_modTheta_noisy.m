


clc;

%% Plot the outcome from diffusion maps

figure;
plot(time, prevtime(:,1));
xlabel('spike time [s]'); 
ylabel('time since the last spike (prevtime) [s]');
title('prevtime vs spiketime, noisestd = eta2*M, tau = 0.1, laps = 4.33');


figure;
scatter(mappedX2(:,1), mappedX2(:,2));
xlabel('eigvec 1'); ylabel('eigvec 2'); 
title('[prev only], noisestd = eta2*M, tau = 0.1, laps = 4.33');



figure;
plot3(time, mappedX2(:,1), mappedX2(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('[prev only], , noisestd = eta2*M, tau = 0.1, laps = 4.33');


figure;
a1=subplot(3,1,1);
plot(time, mod(theta, 2*pi));
xlabel('time [s]');ylabel('mod(theta, 2*pi) in [rad]');
title(a1, 'mod(theta, 2*pi) vs time, noisestd = eta2*M, tau = 0.1, laps = 4.33')


a2 = subplot(3,1,2);
plot(mod(theta, 2*pi), sin(mod(theta, 2*pi)));
xlabel('mod(theta, 2*pi) in [rad]');ylabel('sin(mod(theta, 2*pi)) in [rad]');
title(a2, 'Simulated rat position: mod(theta) Vs Sin(mod(theta))');

a3 = subplot(3,1,3);
plot(time, sin(mod(theta, 2*pi)));
xlabel('time[s]');ylabel('sin(mod(theta, 2*pi) in [rad]');
title(a3, 'Simulated rat position: sin(mod(theta)) Vs Time(s)');




%% plot the projected theta angle in feature space

% new_theta = atan2(mappedX2(:,2), mappedX2(:,1));---this reverses the axes

new_theta = atan2(mappedX2(:,1), mappedX2(:,2));


figure;
a1=subplot(3,1,1);
plot(time, new_theta);
xlabel('time [s]');ylabel('mod(newtheta, 2*pi) in [rad]');
title(a1, 'newtheta vs time, noisestd = eta2*M, tau = 0.1, laps = 4.33')


a2 = subplot(3,1,2);
plot(new_theta, sin(new_theta));
xlabel('newtheta [rad]');ylabel('sin(newtheta) in [rad]');
title(a2, 'Projected rat position: newtheta Vs Sin(newtheta)');

a3 = subplot(3,1,3);
plot(time, sin(new_theta));
xlabel('time[s]');ylabel('sin(newtheta) [rad])');
title(a3, 'Projected rat position: sin(newtheta) Vs Time(s)');


figure;

a1 = subplot(2,1,1);
scatter3(time, mod(theta, 2*pi), sin(mod(theta, 2*pi)));
xlabel('time[s]'); ylabel('mod(theta) [rad]'); zlabel('sin(mod(theta))');
title(a1,'Simulated Rat Position, noisestd = eta2*M, tau=0.1, laps = 4.33');


a2 = subplot(2,1,2);
scatter3(time, new_theta, sin(new_theta));
xlabel('time[s]'); ylabel('newtheta [rad]'); zlabel('sin(newtheta)');
title(a2, 'Projected Rat Position, noisestd = eta2*M, tau = 0.1, laps = 4.33');


% %% Compute the 2-norm of the error
% error = theta - new_theta;
% error = norm(error, 2) %compute the two norm of the error

figure;
a4 = subplot(2,1,1);
plot(time, mod(theta, 2*pi));
xlabel('time[s]');ylabel('mod(theta) [rad])');
title(a4, 'SimulatedRatPosition: mod(theta, 2*pi) Vs Time(s)');
a5 = subplot(2,1,2);
plot(time, sin(new_theta));
xlabel('time[s]');ylabel('sin(newtheta) [rad])');
title(a5, 'Projected rat position: sin(newtheta) Vs Time(s)');

%% calculate the pointwise difference between the two angles 
% wrapped between -pi and pi since that matches atan2 matlab function

%difference = angdiff(theta, new_theta);






