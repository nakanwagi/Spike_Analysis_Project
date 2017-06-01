%% Script to visualize how the eigenvectors reconstruct the rat position
% load('R232-2011-10-16-vt.mat'); %this loads the rat positions in 2011 data
% whos

%% Convert the eigenvectors to ts objects
x1 = tsd(time, mappedX2(:,1)); %first eigenvector
y1 = tsd(time, mappedX2(:,2)); %second eigen vector
new_angle = atan2(mappedX2(:,1), mappedX2(:,2));

%load the actual rat positions using load('-vt.mat');
% plot the rat position
figure;
plot3(x.data, y.data, x.data(x.range), 'g');
%xlabel('x');ylabel('y'); zlabel('x.data(x.range)');
%title('Rat position in Real Data');

%plot rat position using the first eigenvector only
hold on;
plot3(x.data, y.data,  x1.data(x.range), 'r');

%plot the rat position using the second eigenvector only
plot3(x.data, y.data,  y1.data(x.range), 'b');

xlabel('x');ylabel('y');zlabel('range of x')

title('Rat pos in green, eigvec 1 in  red and eigvec 2 in blue');
%% Question: How many eigenvectors do we need to reconstruct the 
% actual rat position???  Use the mutual information!

%% To compute mutual info, we need to do a histogram-based estimation
% method to estimate the joint estimated probability distribution of the 
% the two random variables



