% Script by Mary Sylvia Agwang

%% Load data from the 1TB projects space

load('/project/agwan003/Sim_noisyFRpca.mat');

%% Mutual  Information  for PCA on Simulated Firing Rate
% First Convert the eigenvectors to tsd objects stored in cells
% EV stands for eigen vector


FRpcaEV = cell(1, size(FRpca,2)); %FRpca is pca on simulated noisy FireRate
for i = 1:size(FRpca, 2) % i runs over the columns of mappedX2
    FRpcaEV{i} = tsd(time, FRpca(:,i));
end

figure;
a1 = subplot(2, 2, 1);
plot3(theta, FRpca(:,1), FRpca(:,2));
xlabel('theta');ylabel('eigvec 1');zlabel('eigvec 2');
title(a1, 'SimRatPos and top 2 eigvecs from simFRpca');


% Defines position in terms of cosine and sine.
tsd_theta = tsd(time, theta');
x = cos(mod(tsd_theta.data, 2*pi));
y = sin(mod(tsd_theta.data, 2*pi));


%define max and min for input in histcn
xb = linspace(min(x), max(x), 2*size(FRpca,2)); %range for x
yb = linspace(min(y), max(y), 2*size(FRpca,2)); %range for y


% plot the rat position using only the top eigenvector
a2 = subplot(2, 2, 2);
plot3(time, x, y, 'b');
xlabel('time'); ylabel('cos(modtheta)'); zlabel('sin(modtheta)');
title(a2,'Sim Pos and time:SimFRpca'); %SimFRpca stands for pca on simulated noisy FireRate



a3 = subplot(2, 2, 3);
plot3(x, y, FRpcaEV{1}.data, 'r');
xlabel('cos(modtheta)');ylabel('sin(modtheta)');zlabel('1st eigvec');
title(a3, 'Sim Pos and 1st eigvec');
%projected theta using 1st eig vec from pca on Simulated noisy FireRate 


% Compute the Mutual_Information for each of the 32 eigenvectors
disp('computing Mutual Information for SimFRpca'); tic;

H_simFRpca=cell(1,32);
MI_simFRpca = cell(1,32);
Hnew_simFRpca = cell(1,32);

for iV = 1:size(FRpca, 2);
    % I have only put 10 because it's easy to see the contributions of
    % each eigenvector on submatrix plot
    % EV_FRpca stands for eivectors from pca
EV_FRpca = linspace(min(FRpcaEV{iV}.data), max(FRpcaEV{iV}.data), length(time));
    %compute a 3D hsitogram to estimate the joint entropy
    H_simFRpca{iV} = histcn([x, y, FRpcaEV{iV}.data],xb, yb, EV_FRpca);
    % compute the normalizing constant;

% Matricize the  3D tensor (Histogram)  to speed up MATLAB computations

Hnew_simFRpca{iV} = reshape(H_simFRpca{iV}, [], size(H_simFRpca{iV},3)); %convert to 64^2*length(x.data)
%remove zeros so that log2 is defined
Hnew_simFRpca{iV}(Hnew_simFRpca{iV}==0)=nan;  %set all zero values to nans.


Hnew_simFRpca{iV} = Hnew_simFRpca{iV}./nansum(nansum(nansum(Hnew_simFRpca{iV}))); %normalize the probabilities


lh1 = log2(nansum(Hnew_simFRpca{iV},1));
lh2 = log2(nansum(Hnew_simFRpca{iV},2));
MI_simFRpca{iV} = nansum(nansum(Hnew_simFRpca{iV} .* bsxfun(@minus,bsxfun(@minus,log2(Hnew_simFRpca{iV}),lh1),lh2)));

end

% Create a Mutual Information vector by concatenating all the mutual Information
MI_simFRpca = cat(1, MI_simFRpca{:});


% Now plot the Mutual information
a4 = subplot(2, 2, 4);
p = plot(MI_simFRpca);
p.Marker = 'square';
xlabel('Eigen vector index'); ylabel('Mutual Information');
title(a4, 'Mutual Info Vs Eigindx : SimFRpca');
grid on;

print( 'MI_on_SimNoisyFRpca','-dpdf');

toc;


%% Compute the Mutual Information.

% %% estimate joint histogram by accumulating ,means of using eigvec 1 over rat position
% 
% jointHist1 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{1}.data, 'FUN', @mean); 
% 
% jointHist2 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{2}.data, 'FUN', @mean); 
% 
% jointHist3 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{3}.data, 'FUN', @mean); 
% 
% jointHist4 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{4}.data, 'FUN', @mean); 
% 
% 
% jointHist5 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{5}.data, 'FUN', @mean); 
% 
% 
% jointHist6 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{6}.data, 'FUN', @mean); 
% 
% 
% 
% jointHist7 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{7}.data, 'FUN', @mean); 
% 
% 
% jointHist8 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{8}.data, 'FUN', @mean); 
% 
% 
% jointHist9 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{9}.data, 'FUN', @mean); 
% 
% 
% jointHist10 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{10}.data, 'FUN', @mean); 
% 
% 
% jointHist11 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{11}.data, 'FUN', @mean); 
% 
% 
% jointHist12 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{12}.data, 'FUN', @mean); 
% 
% 
% 
% jointHist13 =  histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{13}.data, 'FUN', @mean); 
% 
% 
% jointHist14 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{14}.data, 'FUN', @mean); 
% 
% 
% jointHist15 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{15}.data, 'FUN', @mean); 
% 
% 
% jointHist16 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{16}.data, 'FUN', @mean); 
% 
% 
% jointHist17 =  histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{17}.data, 'FUN', @mean); 
% 
% 
% 
% jointHist18 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{18}.data, 'FUN', @mean); 
% 
% 
% 
% jointHist19 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{19}.data, 'FUN', @mean); 
% 
% 
% jointHist20 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{20}.data, 'FUN', @mean); 
% 
% 
% jointHist21 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{21}.data, 'FUN', @mean); 
% 
% 
% jointHist22 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{22}.data, 'FUN', @mean); 
% 
% 
% jointHist23 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{23}.data, 'FUN', @mean); 
% 
% 
% jointHist24 =  histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{24}.data, 'FUN', @mean); 
% 
% 
% jointHist25 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{25}.data, 'FUN', @mean); 
% 
% 
% jointHist26 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{26}.data, 'FUN', @mean); 
% 
% 
% jointHist27 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{27}.data, 'FUN', @mean); 
% 
% 
% jointHist28 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{28}.data, 'FUN', @mean); 
% 
% 
% jointHist29 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{29}.data, 'FUN', @mean); 
% 
% 
% jointHist30 = histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{30}.data, 'FUN', @mean); 
% 
% 
% jointHist31 =  histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{31}.data, 'FUN', @mean); 
% 
% 
% jointHist32 =  histcn([x, y], xb, yb,...,
%     'AccumData', eigVec{32}.data, 'FUN', @mean); 
% 
% 


% %% plot the histogram without zero values
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b1 = subplot(2, 2, 1);
% jointHist1(jointHist1==0)= nan; 
% pcolor(jointHist1);
% colorbar;
% title(b1, 'Position and Eigev 1');
% 
% b2 = subplot(2, 2, 2);
% jointHist2(jointHist2==0)= nan; 
% pcolor(jointHist2);
% colorbar;
% title(b2, 'Position and Eigev 2');
% 
% b3 = subplot(2, 2, 3);
% jointHist3(jointHist3==0)= nan; 
% pcolor(jointHist3);
% colorbar;
% title(b3, 'Position and Eigev 3');
% 
% 
% b4 = subplot(2, 2, 4);
% jointHist4(jointHist4==0)= nan; 
% pcolor(jointHist4);
% colorbar;
% title(b4, 'Position and Eigev 4');
% 
% print('sublot1','-dpng');
% 
% 
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b5 = subplot(2, 2, 1);
% jointHist5(jointHist5==0)= nan; 
% pcolor(jointHist5);
% colorbar;
% title(b5, 'Position and Eigev 5');
% 
% b6 = subplot(2, 2, 2);
% jointHist6(jointHist6==0)= nan; 
% pcolor(jointHist6);
% colorbar;
% title(b6, 'Position and Eigev 6');
% 
% b7 = subplot(2, 2, 3);
% jointHist7(jointHist7==0)= nan; 
% pcolor(jointHist7);
% colorbar;
% title(b7, 'Position and Eigev 7');
% 
% 
% b8 = subplot(2, 2, 4);
% jointHist8(jointHist8==0)= nan; 
% pcolor(jointHist8);
% colorbar;
% title(b8, 'Position and Eigev 8');
% 
% print('Subplot 2','-dpng');
% 
% 
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b9 = subplot(2, 2, 1);
% jointHist9(jointHist9==0)= nan; 
% pcolor(jointHist9);
% colorbar;
% title(b9, 'Position and Eigev 9');
% 
% b10 = subplot(2, 2, 2);
% jointHist10(jointHist10==0)= nan; 
% pcolor(jointHist10);
% colorbar;
% title(b10, 'Position and Eigev 10');
% 
% b11 = subplot(2, 2, 3);
% jointHist11(jointHist11==0)= nan; 
% pcolor(jointHist11);
% colorbar;
% title(b11, 'Position and Eigev 11');
% 
% 
% b12 = subplot(2, 2, 4);
% jointHist12(jointHist12==0)= nan; 
% pcolor(jointHist12);
% colorbar;
% title(b12, 'Position and Eigev 12');
% 
% print('Subplot 3','-dpng');
% 
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b13 = subplot(2, 2, 1);
% jointHist13(jointHist13==0)= nan; 
% pcolor(jointHist13);
% colorbar;
% title(b13, 'Position and Eigev 13');
% 
% b14 = subplot(2, 2, 2);
% jointHist14(jointHist14==0)= nan; 
% pcolor(jointHist14);
% colorbar;
% title(b14, 'Position and Eigev 14');
% 
% b15 = subplot(2, 2, 3);
% jointHist15(jointHist15==0)= nan; 
% pcolor(jointHist15);
% colorbar;
% title(b15, 'Position and Eigev 15');
% 
% 
% b16 = subplot(2, 2, 4);
% jointHist16(jointHist16==0)= nan; 
% pcolor(jointHist16);
% colorbar;
% title(b16, 'Position and Eigev 16');
% 
% print('Subplot 4','-dpng');
% 
% 
% 
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b17 = subplot(2, 2, 1);
% jointHist17(jointHist17==0)= nan; 
% pcolor(jointHist17);
% colorbar;
% title(b17, 'Position and Eigev 17');
% 
% b18 = subplot(2, 2, 2);
% jointHist18(jointHist18==0)= nan; 
% pcolor(jointHist18);
% colorbar;
% title(b18, 'Position and Eigev 18');
% 
% b19 = subplot(2, 2, 3);
% jointHist19(jointHist19==0)= nan; 
% pcolor(jointHist19);
% colorbar;
% title(b19, 'Position and Eigev 19');
% 
% 
% b20 = subplot(2, 2, 4);
% jointHist20(jointHist20==0)= nan; 
% pcolor(jointHist20);
% colorbar;
% title(b20, 'Position and Eigev 20');
% 
% print('Subplot 5','-dpng');
% 
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b21 = subplot(2, 2, 1);
% jointHist21(jointHist21==0)= nan; 
% pcolor(jointHist21);
% colorbar;
% title(b21, 'Position and Eigev 21');
% 
% b22 = subplot(2, 2, 2);
% jointHist22(jointHist22==0)= nan; 
% pcolor(jointHist22);
% colorbar;
% title(b22, 'Position and Eigev 22');
% 
% b23 = subplot(2, 2, 3);
% jointHist23(jointHist23==0)= nan; 
% pcolor(jointHist23);
% colorbar;
% title(b23, 'Position and Eigev 23');
% 
% 
% b24 = subplot(2, 2, 4);
% jointHist24(jointHist24==0)= nan; 
% pcolor(jointHist24);
% colorbar;
% title(b24, 'Position and Eigev 24');
% 
% print('Subplot 6','-dpng');
% 
% 
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b25 = subplot(2, 2, 1);
% jointHist25(jointHist25==0)= nan; 
% pcolor(jointHist25);
% colorbar;
% title(b25, 'Position and Eigev 25');
% 
% b26 = subplot(2, 2, 2);
% jointHist26(jointHist26==0)= nan; 
% pcolor(jointHist26);
% colorbar;
% title(b26, 'Position and Eigev 26');
% 
% b27 = subplot(2, 2, 3);
% jointHist27(jointHist27==0)= nan; 
% pcolor(jointHist27);
% colorbar;
% title(b27, 'Position and Eigev 27');
% 
% 
% b28 = subplot(2, 2, 4);
% jointHist28(jointHist28==0)= nan; 
% pcolor(jointHist28);
% colorbar;
% title(b28, 'Position and Eigev 28');
% 
% print('Subplot 7','-dpng');
% 
% 
% figure;
% %set zero values in jointhistogram to nan so log2 is defined
% b29 = subplot(2, 2, 1);
% jointHist29(jointHist29==0)= nan; 
% pcolor(jointHist29);
% colorbar;
% title(b29, 'Position and Eigev 29');
% 
% 
% b30 = subplot(2, 2, 2);
% jointHist30(jointHist30==0)= nan; 
% pcolor(jointHist30);
% colorbar;
% title(b30, 'Position and Eigev 30');
% 
% b31 = subplot(2, 2, 3);
% jointHist31(jointHist31==0)= nan; 
% pcolor(jointHist31);
% colorbar;
% title(b31, 'Position and Eigev 31');
% 
% 
% b32 = subplot(2, 2, 4);
% jointHist32(jointHist32==0)= nan; 
% pcolor(jointHist32);
% colorbar;
% title(b32, 'Position and Eigev 32');
% 
% print('Subplot 8','-dpng');
% 



