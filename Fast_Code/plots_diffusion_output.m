
close all;
clc;

%% Plot the outcome from diffusion maps

figure;
plot(time, prevtime(:,1));
xlabel('spike time [s]'); 
ylabel('time since the last spike (prevtime) [s]');
title('plot of previous time vs spike time');



figure;
plot(time, nextime(:,1));
xlabel('spike time [s]'); 
ylabel('time until the next spike (nextime) [s]');
title('plot of nextime vs spike time');



figure;
A = [prevtime nextime];
plot(time, A(:,1));
xlabel('spike time [s]'); 
ylabel('concatenation of [prevtime nextime]');
title('plot of concatenation of prevtime_nextime vs spike time');



figure
scatter(mappedX1(:,1), mappedX1(:,2));
xlabel('eigvec 1 '); ylabel('eigvec 2'); 
title('first 2 eigen vectors [prev next]');


figure;
scatter(mappedX2(:,1), mappedX2(:,2));
xlabel('eigvec 1'); ylabel('eigvec 2'); 
title('first 2 eigen vectors [prev only]');

figure;
scatter(mappedX3(:,1), mappedX3(:,2));
xlabel('eigevec 1'); ylabel('eigvec 2');
title('first 2 eigen vectors [nextime only]');

figure;
plot3(time, mappedX1(:,1), mappedX1(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev next]');


figure;
plot3(time, mappedX2(:,1), mappedX2(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev only]');
