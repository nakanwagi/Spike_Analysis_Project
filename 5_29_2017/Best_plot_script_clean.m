


clc;

%% Plot the outcome from diffusion maps

figure;
plot(time, prevtime(:,1));
xlabel('spike time [s]'); 
ylabel('time since the last spike (prevtime) [s]');
title('prevtime vs spiketime, clean, tau = 0.9, laps = 4.33');


figure;
scatter(mappedX2(:,1), mappedX2(:,2));
xlabel('eigvec 1'); ylabel('eigvec 2'); 
title('[prev only], clean, tau = 0.9, laps = 4.33');



figure;
plot3(time, mappedX2(:,1), mappedX2(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('[prev only], clean, tau = 0.9, laps = 4.33');


figure;
a1=subplot(3,1,1);
plot(time, mod(theta, 2*pi));
xlabel('time [s]');ylabel('mod(theta, 2*pi) in [rad]');
title(a1, 'mod(theta, 2*pi) vs time, clean, tau = 0.9, laps = 4.33')


a2 = subplot(3,1,2);
plot(theta, sin(theta));
xlabel('theta in [rad]');ylabel('sin(theta) in [rad]');
title(a2, 'Simulated rat position: theta Vs Sin(theta)');

a3 = subplot(3,1,3);
plot(time, sin(theta));
xlabel('time[s]');ylabel('sin(theta[rad])');
title(a3, 'Simulated rat position: sin(theta) Vs Time(s)');




%% plot the projected theta angle in feature space

% new_theta = atan2(mappedX2(:,2), mappedX2(:,1));---this reverses the axes

new_theta = atan2(mappedX2(:,1), mappedX2(:,2));


figure;
a1=subplot(3,1,1);
plot(time, mod(new_theta, 2*pi));
xlabel('time [s]');ylabel('mod(newtheta, 2*pi) in [rad]');
title(a1, 'mod(newtheta, 2*pi) vs time, clean, tau = 0.9, laps = 4.33')


a2 = subplot(3,1,2);
plot(new_theta, sin(new_theta));
xlabel('newtheta [rad]');ylabel('sin(newtheta) in [rad]');
title(a2, 'Projected rat position: newtheta Vs Sin(newtheta)');

a3 = subplot(3,1,3);
plot(time, sin(new_theta));
xlabel('time[s]');ylabel('sin(newtheta) in [rad])');
title(a3, 'Projected rat position: sin(newtheta) Vs Time(s)');


figure;

a1 = subplot(2,1,1);
scatter3(time,theta, sin(theta));
xlabel('time[s]'); ylabel('theta [rad]'); zlabel('sin(theta)');
title(a1,'Simulated Rat Position, clean, tau=0.9, laps = 4.33');


a2 = subplot(2,1,2);
scatter3(time, new_theta, sin(new_theta));
xlabel('time[s]'); ylabel('newtheta [rad]'); zlabel('sin(newtheta)');
title(a2, 'Projected Rat Position, clean, tau = 0.9, laps = 4.33');


% %% Compute the 2-norm of the error
% error = theta - new_theta;
% error = norm(error, 2) %compute the two norm of the error




