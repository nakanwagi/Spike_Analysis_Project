%% Script to visualize how the eigenvectors reconstruct the rat position
% load('R232-2011-10-16-vt.mat'); %this loads the rat positions in 2011 data
% whos

%% Convert the eigenvectors to ts objects
x1 = tsd(time, mappedX2(:,1)); %first largest eigenvector
y1 = tsd(time, mappedX2(:,2)); %second largest eigen vector
new_angle = atan2(mappedX2(:,1), mappedX2(:,2));

%%Make a scree plot to see if you can guess the effective dimension
lambda = mapping2.val; % eigenvalues of the kernel matrix sorted in descending order.
ind = mapping2.index; %indices of the eigenvalues
figure;
plot(ind, lambda, 'ko-');
xlabel('index');ylabel('eigenvalue');
title('Scree plot: 2011Data with 32 cells');


%% plot the rat position using only the top eigenvector
figure;
plot3(x.data, y.data, x.data(x.range), 'g');
hold on;
plot3(x.data, y.data,  x1.data(x.range), 'r');
xlabel('x');ylabel('y');zlabel('range of x')
title('Rat position in green, top eigenvector only in  red: 2011Data');
hold off;

%% ???To compute mutual info, we need to do a histogram-based estimation
% method to estimate the joint estimated probability distribution of the 
% the two random variables

%FORMULA: I(X,Y) = H(X) + H(Y) - H(X,Y)
% H(x) = -sum p(x) log2(p(x))---marginal pdf of x.
% why log 2 and not ln?

bound1 = linspace(min(x.data), max(x.data), 80); %range for x
bound2 = linspace(min(y.data), max(y.data), 80); %range for y

jointHist = histcn([x.data y.data], bound1, bound2,...,
    'AccumData', x1.data(x.range), 'FUN', @mean); %estimate joint pdf using eigvec 1

%make log2 defined by setting zero values in the jointHist equal to nan
% jointHistogram(jointHistogram == 0) = nan; 
idx = unique(jointHist);


%plot the histogram without zero values
figure;
jointHist(jointHist==0)= nan; %set zero values in jointhistogram to nan
pcolor(jointHist);
colorbar;
title('Rat Position based on the first largest eigenvector only');


%%  Find the joint probability distribution with zero values removed.
jointProbDist = jointHist./nansum(jointHist); %joint probability distribution
jointEntropy = -nansum(jointHist(:).* log2(jointHist(:)));

% can one assume that X=[x.data y.data] is in rows so that
% marginal = nansum(JointHist, 1); and then repeat the computation?

%most code assumes vectors contain integers but mine are continuous
%variables.
% is your code set up to address that
% are you familair with spiketrain analysis tool?





