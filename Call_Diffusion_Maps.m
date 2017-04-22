
disp('plotting the graphs'); tic;

figure;
scatter(mappedX1(:,1), mappedX1(:,2)); 
title('DiffusionMaps on [exp(-abs(prevtime)) exp(-abs(prevtime))]');
xlabel('eigvec 1'); ylabel('eigvec 2');


figure;
scatter(mappedX2(:,1), mappedX2(:,2)); 
title('DiffusionMaps on [exp(-abs(prevtime))] only');
xlabel('eigvec 1'); ylabel('eigvec 2');


figure;
scatter(mappedX3(:,1), mappedX3(:,2)); 
title('DiffusionMaps on [exp(-abs(nextime))] only');
xlabel('eigvec 1'); ylabel('eigvec 2');





figure
plot3(mod(theta, 2*pi), mappedX1(:,1), mappedX1(:,2));
xlabel('theta'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('theta and first 2 eigen vectors [prev next]');


figure;
plot3(mod(theta, 2*pi), mappedX2(:,1), mappedX2(:,2));
xlabel('theta'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('theta and first 2 eigen vectors [prev only]');

figure;
plot3(mod(theta, 2*pi), mappedX3(:,1), mappedX3(:,2));
xlabel('theta'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('theta and first 2 eigen vectors [nextime only]');

toc;








