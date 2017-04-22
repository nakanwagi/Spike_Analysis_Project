
disp('recompute prev/next to remove nans from prevtime');tic;
%% ==========run this part only if there are nan values==================
%initialize prevtime and nextime
nextime = nan(length(time), ncols);
prevtime = nan(length(time), ncols);
idx = cell(length(time), ncols);
idx2 = cell(length(time), ncols);
tmp = cell(length(time), ncols);
tmp2 = cell(length(time), ncols);
 
 %% Comoute prevtime and nextime for the 32 spike trains
 for i = 1:nrows
     for j = 1:ncols
          idx{i,j} = sd.S{j}.T < time(i);  % get logicals from the jth spike train corresponding to i^th time point.
          idx2{i,j} = sd.S{j}.T > time(i);  %get logicals as above for nexttime.
          
          if sum(idx{i,j})==0  %set prevtime = nan if idx is empty
              prevtime(i,j) = nan;
          else
              tmp{i,j} = sd.S{j}.T(idx{i,j}); %compute prevtime if idx is not empty
          prevtime(i,j) = tmp{i,j}(end) - time(i);
          
          end
          
         if sum(idx2{i,j})==0 %set nextime = nan if idx2 is empty
            nextime(i,j) = nan;
          else
              tmp2{i,j} = sd.S{j}.T(idx2{i,j}); %compute prevtime if idx is not empty
              nextime(i,j) = tmp2{i,j}(1) - time(i);
          
          end        
          
     end
 end

 toc;


%% Run diffusion Maps Alagorithm
disp('run diffusion Maps Algorithm');tic;


X1 = [exp(-abs(prevtime)) exp(-abs(nextime))];
X2 = exp(-abs(prevtime));
X3 = exp(-abs(nextime));

%X = exp(-NextTime);

%=============Use the Dim reduction tool box =============================



% call diffusion maps algorithm on [exp-(abs(pretime)) exp-(abs(nexttime))]

[mappedX1, mapping1] = compute_mapping(X1, 'DiffusionMaps');

% call diffusion maps algorithm on [-exp(prevtime)] only
 
[mappedX2, mapping2] = compute_mapping(X2, 'DiffusionMaps');

 
% call diffusion maps algorithm on [-exp(nextime)] only
 
[mappedX3, mapping3] = compute_mapping(X3, 'DiffusionMaps');


toc;


%% Plot the outcome

disp('plot output from diffusion maps');tic;


% figure;
% scatter(mappedX1(:,1), mappedX1(:,2)); 
% title('DiffusionMaps on [exp(-abs(prevtime)) exp(-abs(nextime))]');
% xlabel('eigvec 1'); ylabel('eigvec 2');
% 
% 
% figure;
% scatter(mappedX2(:,1), mappedX2(:,2)); 
% title('DiffusionMaps on [exp(-abs(prevtime))] only');
% xlabel('eigvec 1'); ylabel('eigvec 2');
% 
% 
% figure;
% scatter(mappedX3(:,1), mappedX3(:,2)); 
% title('DiffusionMaps on [exp(-abs(nextime))] only');
% xlabel('eigvec 1'); ylabel('eigvec 2');
% 

% figure;
% plot(time, prevtime(:,1));
% xlabel('spike time [s]');ylabel('time since the last spike (prevtime) [s]');
% title('plot of previous time vs time');


%% Plot the outcome

% disp('plot output from diffusion maps');tic;


% figure;
% scatter(mappedX1(:,1), mappedX1(:,2)); 
% title('DiffusionMaps on [exp(-abs(prevtime)) exp(-abs(nextime))]');
% xlabel('eigvec 1'); ylabel('eigvec 2');
% 
% 
% figure;
% scatter(mappedX2(:,1), mappedX2(:,2)); 
% title('DiffusionMaps on [exp(-abs(prevtime))] only');
% xlabel('eigvec 1'); ylabel('eigvec 2');
% 
% 
% figure;
% scatter(mappedX3(:,1), mappedX3(:,2)); 
% title('DiffusionMaps on [exp(-abs(nextime))] only');
% xlabel('eigvec 1'); ylabel('eigvec 2');
% 



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
plot3(time, mappedX1(:,1), mappedX1(:,2));
xlabel('time'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev next]');


figure;
plot3(time, mappedX2(:,1), mappedX2(:,2));
xlabel('time'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev only]');

figure;
plot3(time, mappedX3(:,1), mappedX3(:,2));
xlabel('time'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [nextime only]');

figure;
plot3(time, mappedX1(:,1), mappedX1(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev next]');


figure;
plot3(time, mappedX2(:,1), mappedX2(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [prev only]');

figure;
plot3(time, mappedX3(:,1), mappedX3(:,2));
xlabel('spike time [s]'); ylabel('eigvec 1'); zlabel('eigvec 2'); 
title('time and first 2 eigen vectors [nextime only]');










