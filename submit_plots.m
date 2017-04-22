
close all;




figure;
plot(time, prevtime(:,1));
xlabel('spike time [s]'); 
ylabel('time since the last spike (prevtime) [s]');
title('plot of previous time vs spike time');
print('prevtime', '-depsc');
print('prevtime', '-dpng');



figure;
plot(time, nextime(:,1));
xlabel('spike time [s]'); 
ylabel('time until the next spike (nextime) [s]');
title('plot of nextime vs spike time');
print('nextime', '-depsc');
print('nextime', '-dpng');



figure;
A = [prevtime nextime];
plot(time, A(:,1));
xlabel('spike time [s]'); 
ylabel('concatenation of [prevtime nextime]');
title('plot of concatenation of prevtime_nextime vs spike time');
print('prevnext', '-depsc');
print('prevnext', '-dpng');



figure
plot3(time, mappedX1(:,1), mappedX1(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev next]');
grid on;
print('nextPrevEigvec', '-depsc')
print('nextPrevEigvec', '-dpng')


figure;
plot3(time, mappedX2(:,1), mappedX2(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev only]');
grid on;
print('prevEigvec', '-depsc')
print('prevEigvec', '-dpng')


figure;
plot3(time, mappedX3(:,1), mappedX3(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [nextime only]');
grid on;
print('nextimeEigV', '-depsc')
print('nextimeEigV', '-dpng')




