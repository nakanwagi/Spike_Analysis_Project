% Script by Mary Sylvia Agwang
%% Make a scree plot to see if you can guess the effective dimension
lambda = fireRate_struct.val; % eigenvalues of the kernel matrix sorted in descending order.
ind = fireRate_struct.index; %indices of the eigenvalues

% figure;
% plot(ind, lambda, 'ko-');
% xlabel('index');ylabel('eigenvalue');
% title('Scree plot: 2011Data with 32 cells');
% grid on; 

% %print the plot
% print('Scree Plot','-dpng');

%% find out how much variance the first k eigenvectors explain
% k=2;
% eigvec_no = 2;
% explained_variance = 100*cumsum(lambda(1:eigvec_no))/sum(lambda(1:eigvec_no));
% %%var = mean(var)*100; %percentage of variance explained by first k eigenvectors

%% First Convert the eigenvectors to tsd objects stored in cells
eigVec = cell(1, size(fireRate1,2));
for i = 1:size(fireRate1, 2) % i runs over the columns of mappedX2
    eigVec{i} = tsd(time, fireRate1(:,i));
end


% new_angle = atan2(mappedX2(:,1), mappedX2(:,2));

%% plot the rat position using only the top eigenvector
figure;
plot3(x.data, y.data, x.data(x.range), 'b');
hold on;
plot3(x.data, y.data,  eigVec{1}.data(x.range), 'r');
xlabel('x');ylabel('y');zlabel('First Eigen Vector');
title('Rat position in blue, first eigenvector only in  red: 2011Data on FireR');

% print the plot
print('Rat Position and EigVector 1 on FireR','-dpng');



xb = linspace(min(x.data), max(x.data), 2*size(fireRate1,2)); %range for x
yb = linspace(min(y.data), max(y.data), 2*size(fireRate1,2)); %range for y



%% estimate joint histogram by accumulating ,means of using eigvec 1 over rat position

FjointHist1 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{1}.data(x.range), 'FUN', @mean); 

FjointHist2 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{2}.data(x.range), 'FUN', @mean); 

FjointHist3 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{3}.data(x.range), 'FUN', @mean);

FjointHist4 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{4}.data(x.range), 'FUN', @mean);

FjointHist5 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{5}.data(x.range), 'FUN', @mean); 

FjointHist6 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{6}.data(x.range), 'FUN', @mean); 


FjointHist7 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{7}.data(x.range), 'FUN', @mean); 

FjointHist8 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{8}.data(x.range), 'FUN', @mean); 


FjointHist9 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{9}.data(x.range), 'FUN', @mean); 



FjointHist10 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{10}.data(x.range), 'FUN', @mean); 

FjointHist11 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{11}.data(x.range), 'FUN', @mean); 

FjointHist12 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{12}.data(x.range), 'FUN', @mean); 


FjointHist13 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{13}.data(x.range), 'FUN', @mean); 

FjointHist14 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{14}.data(x.range), 'FUN', @mean); 

FjointHist15 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{15}.data(x.range), 'FUN', @mean); 

FjointHist16 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{16}.data(x.range), 'FUN', @mean); 

FjointHist17 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{17}.data(x.range), 'FUN', @mean); 


FjointHist18 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{18}.data(x.range), 'FUN', @mean); 


FjointHist19 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{19}.data(x.range), 'FUN', @mean); 

FjointHist20 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{20}.data(x.range), 'FUN', @mean); 

FjointHist21 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{21}.data(x.range), 'FUN', @mean); 

FjointHist22 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{22}.data(x.range), 'FUN', @mean); 

FjointHist23 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{23}.data(x.range), 'FUN', @mean); 

FjointHist24 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{24}.data(x.range), 'FUN', @mean); 

FjointHist25 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{25}.data(x.range), 'FUN', @mean); 

FjointHist26 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{26}.data(x.range), 'FUN', @mean); 

FjointHist27 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{27}.data(x.range), 'FUN', @mean); 

FjointHist28 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{28}.data(x.range), 'FUN', @mean); 

FjointHist29 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{29}.data(x.range), 'FUN', @mean); 

FjointHist30 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{30}.data(x.range), 'FUN', @mean); 

FjointHist31 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{31}.data(x.range), 'FUN', @mean); 

FjointHist32 = histcn([x.data y.data], xb, yb,...,
    'AccumData', eigVec{32}.data(x.range), 'FUN', @mean); 



%% plot the histogram without zero values
figure;
%set zero values in jointhistogram to nan so log2 is defined
b1 = subplot(2, 2, 1);
FjointHist1(FjointHist1==0)= nan; 
pcolor(FjointHist1);
colorbar;
title(b1, 'Position and Eigev 1 on FireR');

b2 = subplot(2, 2, 2);
FjointHist2(FjointHist2==0)= nan; 
pcolor(FjointHist2);
colorbar;
title(b2, 'Position and Eigev 2 on FireR');

b3 = subplot(2, 2, 3);
FjointHist3(FjointHist3==0)= nan; 
pcolor(FjointHist3);
colorbar;
title(b3, 'Position and Eigev 3 on FireR');


b4 = subplot(2, 2, 4);
FjointHist4(FjointHist4==0)= nan; 
pcolor(FjointHist4);
colorbar;
title(b4, 'Position and Eigev 4 on FireR');

print('sublot1 on FireR','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b5 = subplot(2, 2, 1);
FjointHist5(FjointHist5==0)= nan; 
pcolor(FjointHist5);
colorbar;
title(b5, 'Position and Eigev 5 on FireR');

b6 = subplot(2, 2, 2);
FjointHist6(FjointHist6==0)= nan; 
pcolor(FjointHist6);
colorbar;
title(b6, 'Position and Eigev 6 on FireR');

b7 = subplot(2, 2, 3);
FjointHist7(FjointHist7==0)= nan; 
pcolor(FjointHist7);
colorbar;
title(b7, 'Position and Eigev 7 on FireR');


b8 = subplot(2, 2, 4);
FjointHist8(FjointHist8==0)= nan; 
pcolor(FjointHist8);
colorbar;
title(b8, 'Position and Eigev 8 on FireR');

print('Subplot 2 on FireR','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b9 = subplot(2, 2, 1);
FjointHist9(FjointHist9==0)= nan; 
pcolor(FjointHist9);
colorbar;
title(b9, 'Position and Eigev 9 on FireR');

b10 = subplot(2, 2, 2);
FjointHist10(FjointHist10==0)= nan; 
pcolor(FjointHist10);
colorbar;
title(b10, 'Position and Eigev 10 on FireR');

b11 = subplot(2, 2, 3);
FjointHist11(FjointHist11==0)= nan; 
pcolor(FjointHist11);
colorbar;
title(b11, 'Position and Eigev 11 on FireR');


b12 = subplot(2, 2, 4);
FjointHist12(FjointHist12==0)= nan; 
pcolor(FjointHist12);
colorbar;
title(b12, 'Position and Eigev 12 on FireR');

print('Subplot 3 on FireR','-dpng');

figure;
%set zero values in jointhistogram to nan so log2 is defined
b13 = subplot(2, 2, 1);
FjointHist13(FjointHist13==0)= nan; 
pcolor(FjointHist13);
colorbar;
title(b13, 'Position and Eigev 13 on FireR');

b14 = subplot(2, 2, 2);
FjointHist14(FjointHist14==0)= nan; 
pcolor(FjointHist14);
colorbar;
title(b14, 'Position and Eigev 14 on FireR');

b15 = subplot(2, 2, 3);
FjointHist15(FjointHist15==0)= nan; 
pcolor(FjointHist15);
colorbar;
title(b15, 'Position and Eigev 15 on FireR');


b16 = subplot(2, 2, 4);
FjointHist16(FjointHist16==0)= nan; 
pcolor(FjointHist16);
colorbar;
title(b16, 'Position and Eigev 16 on FireR');

print('Subplot 4 on FireR','-dpng');



figure;
%set zero values in jointhistogram to nan so log2 is defined
b17 = subplot(2, 2, 1);
FjointHist17(FjointHist17==0)= nan; 
pcolor(FjointHist17);
colorbar;
title(b17, 'Position and Eigev 17 on FireR');

b18 = subplot(2, 2, 2);
FjointHist18(FjointHist18==0)= nan; 
pcolor(FjointHist18);
colorbar;
title(b18, 'Position and Eigev 18 on FireR');

b19 = subplot(2, 2, 3);
FjointHist19(FjointHist19==0)= nan; 
pcolor(FjointHist19);
colorbar;
title(b19, 'Position and Eigev 19 on FireR');


b20 = subplot(2, 2, 4);
FjointHist20(FjointHist20==0)= nan; 
pcolor(FjointHist20);
colorbar;
title(b20, 'Position and Eigev 20 on FireR');

print('Subplot 5 on FireR','-dpng');

figure;
%set zero values in jointhistogram to nan so log2 is defined
b21 = subplot(2, 2, 1);
FjointHist21(FjointHist21==0)= nan; 
pcolor(FjointHist21);
colorbar;
title(b21, 'Position and Eigev 21 on FireR');

b22 = subplot(2, 2, 2);
FjointHist22(FjointHist22==0)= nan; 
pcolor(FjointHist22);
colorbar;
title(b22, 'Position and Eigev 22 on FireR');

b23 = subplot(2, 2, 3);
FjointHist23(FjointHist23==0)= nan; 
pcolor(FjointHist23);
colorbar;
title(b23, 'Position and Eigev 23 on FireR');


b24 = subplot(2, 2, 4);
FjointHist24(FjointHist24==0)= nan; 
pcolor(FjointHist24);
colorbar;
title(b24, 'Position and Eigev 24 on FireR');

print('Subplot 6 on FireR','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b25 = subplot(2, 2, 1);
FjointHist25(FjointHist25==0)= nan; 
pcolor(FjointHist25);
colorbar;
title(b25, 'Position and Eigev 25 on FireR');

b26 = subplot(2, 2, 2);
FjointHist26(FjointHist26==0)= nan; 
pcolor(FjointHist26);
colorbar;
title(b26, 'Position and Eigev 26 on FireR');

b27 = subplot(2, 2, 3);
FjointHist27(FjointHist27==0)= nan; 
pcolor(FjointHist27);
colorbar;
title(b27, 'Position and Eigev 27 on FireR');


b28 = subplot(2, 2, 4);
FjointHist28(FjointHist28==0)= nan; 
pcolor(FjointHist28);
colorbar;
title(b28, 'Position and Eigev 28 on FireR');

print('Subplot 7 on FireR','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b29 = subplot(2, 2, 1);
FjointHist29(FjointHist29==0)= nan; 
pcolor(FjointHist29);
colorbar;
title(b29, 'Position and Eigev 29 on FireR');


b30 = subplot(2, 2, 2);
FjointHist30(FjointHist30==0)= nan; 
pcolor(FjointHist30);
colorbar;
title(b30, 'Position and Eigev 30 on FireR');

b31 = subplot(2, 2, 3);
FjointHist31(FjointHist31==0)= nan; 
pcolor(FjointHist31);
colorbar;
title(b31, 'Position and Eigev 31 on FireR');


b32 = subplot(2, 2, 4);
FjointHist32(FjointHist32==0)= nan; 
pcolor(FjointHist32);
colorbar;
title(b32, 'Position and Eigev 32  on FireR');

print('Subplot 8 FireR','-dpng');

% %% Recover the Rat position using the first eigenvector only.
% 
% xb = linspace(min(x.data), max(x.data), 2*size(mappedX2,2)); %range for x
% yb = linspace(min(y.data), max(y.data), 2*size(mappedX2,2)); %range for y




%% Compute the Mutual_Information for each of the 32 eigenvectors

iV = 32;
% I have only put 10 because it's easy to see the contributions of
% each eigenvector on submatrix plot
eigVecb = linspace(min(eigVec{iV}.data), max(eigVec{iV}.data), length(x.range));
%compute a 3D hsitogram to estimate the joint entropy
H = histcn([x.data, y.data, eigVec{iV}.data(x.range)],xb, yb, eigVecb);
% compute the normalizing constant;

% Matricize the tensor to speed up MATLAB computations

Hnew = reshape(H, [], size(H,3)); %convert to 64^2*length(x.data)
% remove zeros so that log2 is defined
Hnew(Hnew==0)=nan;


Hnew = Hnew./nansum(nansum(nansum(Hnew)));


lh1 = log2(nansum(Hnew,1));
lh2 = log2(nansum(Hnew,2));
MI_Fire32 = nansum(nansum(Hnew .* bsxfun(@minus,bsxfun(@minus,log2(Hnew),lh1),lh2)));

 %% Create a Mutual Information vector.

MI_FireRate = [MI_Fire1, MI_Fire2,  MI_Fire3,  MI_Fire4,...
    MI_Fire5, MI_Fire6,  MI_Fire7, MI_Fire8, MI_Fire9,...
    MI_Fire10, MI_Fire11, MI_Fire12, MI_Fire13, MI_Fire14,...
    MI_Fire15,  MI_Fire16,  MI_Fire17, MI_Fire18,  MI_Fire19,...
    MI_Fire20, MI_Fire21,  MI_Fire22, MI_Fire23,  MI_Fire24,...
   MI_Fire25, MI_Fire26,  MI_Fire27,  MI_Fire28,  MI_Fire29,...
   MI_Fire30, MI_Fire31,  MI_Fire32];


figure;
p = plot(MI_FireRate);
p.Marker = 'square';
xlabel('Eigen vector index'); ylabel('Mutual Information in bits');
title('Mutual Information Vs EigenVector index on FireRate');
grid on;

print('Mutual Information on FireRate','-dpng');





%use I(X,Y) = -sum(sum(p(x,y)log(p(x,y)/(p(x)*p(y)))))
% I(X,Y) = -sum(sum(p(x,y)[log(p(x,y)-log(p(x))-log(p(y))]).






