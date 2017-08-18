% Script by Mary Sylvia Agwang

%% First Convert the eigenvectors to tsd objects stored in cells
FRdiffmapeigVec = cell(1, size(FRdiffmap,2));
for i = 1:size(FRdiffmap, 2) % i runs over the columns of mappedX2
    FRdiffmapeigVec{i} = tsd(time, FRdiffmap(:,i));
end


% new_angle = atan2(mappedX2(:,1), mappedX2(:,2));

%% plot the rat position using only the top eigenvector
figure;
plot3(x.data, y.data, x.data(x.range), 'b');
hold on;
plot3(x.data, y.data,  FRdiffmapeigVec{1}.data(x.range), 'r');
xlabel('x');ylabel('y');zlabel('First Eigen Vector');
title('RatPos in blue, first eigenvector only in  red: diffmaps on FR');

% print the plot
print('Rat Position and EigVector 1 on FRdiffmaps','-dpng');



xb = linspace(min(x.data), max(x.data), 2*size(FRdiffmap,2)); %range for x
yb = linspace(min(y.data), max(y.data), 2*size(FRdiffmap,2)); %range for y



%% estimate joint histogram by accumulating ,means of using eigvec 1 over rat position

FjointHist1 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{1}.data(x.range), 'FUN', @mean); 

FjointHist2 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{2}.data(x.range), 'FUN', @mean); 

FjointHist3 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{3}.data(x.range), 'FUN', @mean);

FjointHist4 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{4}.data(x.range), 'FUN', @mean);

FjointHist5 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{5}.data(x.range), 'FUN', @mean); 

FjointHist6 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{6}.data(x.range), 'FUN', @mean); 


FjointHist7 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{7}.data(x.range), 'FUN', @mean); 

FjointHist8 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{8}.data(x.range), 'FUN', @mean); 


FjointHist9 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{9}.data(x.range), 'FUN', @mean); 



FjointHist10 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{10}.data(x.range), 'FUN', @mean); 

FjointHist11 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{11}.data(x.range), 'FUN', @mean); 

FjointHist12 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{12}.data(x.range), 'FUN', @mean); 


FjointHist13 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{13}.data(x.range), 'FUN', @mean); 

FjointHist14 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{14}.data(x.range), 'FUN', @mean); 

FjointHist15 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{15}.data(x.range), 'FUN', @mean); 

FjointHist16 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{16}.data(x.range), 'FUN', @mean); 

FjointHist17 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{17}.data(x.range), 'FUN', @mean); 


FjointHist18 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{18}.data(x.range), 'FUN', @mean); 


FjointHist19 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{19}.data(x.range), 'FUN', @mean); 

FjointHist20 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{20}.data(x.range), 'FUN', @mean); 

FjointHist21 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{21}.data(x.range), 'FUN', @mean); 

FjointHist22 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{22}.data(x.range), 'FUN', @mean); 

FjointHist23 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{23}.data(x.range), 'FUN', @mean); 

FjointHist24 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{24}.data(x.range), 'FUN', @mean); 

FjointHist25 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{25}.data(x.range), 'FUN', @mean); 

FjointHist26 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{26}.data(x.range), 'FUN', @mean); 

FjointHist27 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{27}.data(x.range), 'FUN', @mean); 

FjointHist28 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{28}.data(x.range), 'FUN', @mean); 

FjointHist29 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{29}.data(x.range), 'FUN', @mean); 

FjointHist30 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{30}.data(x.range), 'FUN', @mean); 

FjointHist31 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{31}.data(x.range), 'FUN', @mean); 

FjointHist32 = histcn([x.data y.data], xb, yb,...,
    'AccumData', FRdiffmapeigVec{32}.data(x.range), 'FUN', @mean); 



%% plot the histogram without zero values
figure;
%set zero values in jointhistogram to nan so log2 is defined
b1 = subplot(2, 2, 1);
FjointHist1(FjointHist1==0)= nan; 
pcolor(FjointHist1);
colorbar;
title(b1, 'Position and Eigev 1 on FRdiffmap');

b2 = subplot(2, 2, 2);
FjointHist2(FjointHist2==0)= nan; 
pcolor(FjointHist2);
colorbar;
title(b2, 'Position and Eigev 2 on FRdiffmap');

b3 = subplot(2, 2, 3);
FjointHist3(FjointHist3==0)= nan; 
pcolor(FjointHist3);
colorbar;
title(b3, 'Position and Eigev 3 on FRdiffmap');


b4 = subplot(2, 2, 4);
FjointHist4(FjointHist4==0)= nan; 
pcolor(FjointHist4);
colorbar;
title(b4, 'Position and Eigev 4 on FRdiffmap');

print('sublot1 on FRdiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b5 = subplot(2, 2, 1);
FjointHist5(FjointHist5==0)= nan; 
pcolor(FjointHist5);
colorbar;
title(b5, 'Position and Eigev 5 on FRdiffmap');

b6 = subplot(2, 2, 2);
FjointHist6(FjointHist6==0)= nan; 
pcolor(FjointHist6);
colorbar;
title(b6, 'Position and Eigev 6 on FRdiffmap');

b7 = subplot(2, 2, 3);
FjointHist7(FjointHist7==0)= nan; 
pcolor(FjointHist7);
colorbar;
title(b7, 'Position and Eigev 7 on FRdiffmap');


b8 = subplot(2, 2, 4);
FjointHist8(FjointHist8==0)= nan; 
pcolor(FjointHist8);
colorbar;
title(b8, 'Position and Eigev 8 on FRdiffmap');

print('Subplot 2 on FRdiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b9 = subplot(2, 2, 1);
FjointHist9(FjointHist9==0)= nan; 
pcolor(FjointHist9);
colorbar;
title(b9, 'Position and Eigev 9 on FRdiffmap');

b10 = subplot(2, 2, 2);
FjointHist10(FjointHist10==0)= nan; 
pcolor(FjointHist10);
colorbar;
title(b10, 'Position and Eigev 10 on FRdiffmap');

b11 = subplot(2, 2, 3);
FjointHist11(FjointHist11==0)= nan; 
pcolor(FjointHist11);
colorbar;
title(b11, 'Position and Eigev 11 on FRdiffmap');


b12 = subplot(2, 2, 4);
FjointHist12(FjointHist12==0)= nan; 
pcolor(FjointHist12);
colorbar;
title(b12, 'Position and Eigev 12 on FRdiffmap');

print('Subplot 3 on FRdiffmap','-dpng');

figure;
%set zero values in jointhistogram to nan so log2 is defined
b13 = subplot(2, 2, 1);
FjointHist13(FjointHist13==0)= nan; 
pcolor(FjointHist13);
colorbar;
title(b13, 'Position and Eigev 13 on FRdiffmap');

b14 = subplot(2, 2, 2);
FjointHist14(FjointHist14==0)= nan; 
pcolor(FjointHist14);
colorbar;
title(b14, 'Position and Eigev 14 on FRdiffmap');

b15 = subplot(2, 2, 3);
FjointHist15(FjointHist15==0)= nan; 
pcolor(FjointHist15);
colorbar;
title(b15, 'Position and Eigev 15 on FRdiffmap');


b16 = subplot(2, 2, 4);
FjointHist16(FjointHist16==0)= nan; 
pcolor(FjointHist16);
colorbar;
title(b16, 'Position and Eigev 16 on FRdiffmap');

print('Subplot 4 on FRdiffmap','-dpng');



figure;
%set zero values in jointhistogram to nan so log2 is defined
b17 = subplot(2, 2, 1);
FjointHist17(FjointHist17==0)= nan; 
pcolor(FjointHist17);
colorbar;
title(b17, 'Position and Eigev 17 on FRdiffmap');

b18 = subplot(2, 2, 2);
FjointHist18(FjointHist18==0)= nan; 
pcolor(FjointHist18);
colorbar;
title(b18, 'Position and Eigev 18 on FRdiffmap');

b19 = subplot(2, 2, 3);
FjointHist19(FjointHist19==0)= nan; 
pcolor(FjointHist19);
colorbar;
title(b19, 'Position and Eigev 19 on FRdiffmap');


b20 = subplot(2, 2, 4);
FjointHist20(FjointHist20==0)= nan; 
pcolor(FjointHist20);
colorbar;
title(b20, 'Position and Eigev 20 on FRdiffmap');

print('Subplot 5 on FRdiffmap','-dpng');

figure;
%set zero values in jointhistogram to nan so log2 is defined
b21 = subplot(2, 2, 1);
FjointHist21(FjointHist21==0)= nan; 
pcolor(FjointHist21);
colorbar;
title(b21, 'Position and Eigev 21 on FRdiffmap');

b22 = subplot(2, 2, 2);
FjointHist22(FjointHist22==0)= nan; 
pcolor(FjointHist22);
colorbar;
title(b22, 'Position and Eigev 22 on FRdiffmap');

b23 = subplot(2, 2, 3);
FjointHist23(FjointHist23==0)= nan; 
pcolor(FjointHist23);
colorbar;
title(b23, 'Position and Eigev 23 on FRdiffmap');


b24 = subplot(2, 2, 4);
FjointHist24(FjointHist24==0)= nan; 
pcolor(FjointHist24);
colorbar;
title(b24, 'Position and Eigev 24 on FRdiffmap');

print('Subplot 6 on FRdiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b25 = subplot(2, 2, 1);
FjointHist25(FjointHist25==0)= nan; 
pcolor(FjointHist25);
colorbar;
title(b25, 'Position and Eigev 25 on FRdiffmap');

b26 = subplot(2, 2, 2);
FjointHist26(FjointHist26==0)= nan; 
pcolor(FjointHist26);
colorbar;
title(b26, 'Position and Eigev 26 on FRdiffmap');

b27 = subplot(2, 2, 3);
FjointHist27(FjointHist27==0)= nan; 
pcolor(FjointHist27);
colorbar;
title(b27, 'Position and Eigev 27 on FRdiffmap');


b28 = subplot(2, 2, 4);
FjointHist28(FjointHist28==0)= nan; 
pcolor(FjointHist28);
colorbar;
title(b28, 'Position and Eigev 28 on FRdiffmap');

print('Subplot 7 on FRdiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b29 = subplot(2, 2, 1);
FjointHist29(FjointHist29==0)= nan; 
pcolor(FjointHist29);
colorbar;
title(b29, 'Position and Eigev 29 on FRdiffmap');


b30 = subplot(2, 2, 2);
FjointHist30(FjointHist30==0)= nan; 
pcolor(FjointHist30);
colorbar;
title(b30, 'Position and Eigev 30 on FRdiffmap');

b31 = subplot(2, 2, 3);
FjointHist31(FjointHist31==0)= nan; 
pcolor(FjointHist31);
colorbar;
title(b31, 'Position and Eigev 31 on FRdiffmap');


b32 = subplot(2, 2, 4);
FjointHist32(FjointHist32==0)= nan; 
pcolor(FjointHist32);
colorbar;
title(b32, 'Position and Eigev 32  on FRdiffmap');

print('Subplot 8 FRdiffmap','-dpng');



%% Compute the Mutual_Information for each of the 32 eigenvectors

iV = 32;
% I have only put 10 because it's easy to see the contributions of
% each eigenvector on submatrix plot
FRdiffmapeigVecb = linspace(min(FRdiffmapeigVec{iV}.data), max(FRdiffmapeigVec{iV}.data), length(x.range));
%compute a 3D hsitogram to estimate the joint entropy
H = histcn([x.data, y.data, FRdiffmapeigVec{iV}.data(x.range)],xb, yb, FRdiffmapeigVecb);
% compute the normalizing constant;

% Matricize the tensor to speed up MATLAB computations

Hnew = reshape(H, [], size(H,3)); %convert to 64^2*length(x.data)
% remove zeros so that log2 is defined
Hnew(Hnew==0)=nan;


Hnew = Hnew./nansum(nansum(nansum(Hnew)));


lh1 = log2(nansum(Hnew,1));
lh2 = log2(nansum(Hnew,2));
MI_FRdiffmap32 = nansum(nansum(Hnew .* bsxfun(@minus,bsxfun(@minus,log2(Hnew),lh1),lh2)));

 %% Create a Mutual Information vector.

MI_FRdiffmap = [MI_FRdiffmap1, MI_FRdiffmap2,  MI_FRdiffmap3,  MI_FRdiffmap4,...
    MI_FRdiffmap5, MI_FRdiffmap6,  MI_FRdiffmap7, MI_FRdiffmap8, MI_FRdiffmap9,...
    MI_FRdiffmap10, MI_FRdiffmap11, MI_FRdiffmap12, MI_FRdiffmap13, MI_FRdiffmap14,...
    MI_FRdiffmap15,  MI_FRdiffmap16,  MI_FRdiffmap17, MI_FRdiffmap18,  MI_FRdiffmap19,...
    MI_FRdiffmap20, MI_FRdiffmap21,  MI_FRdiffmap22, MI_FRdiffmap23,  MI_FRdiffmap24,...
   MI_FRdiffmap25, MI_FRdiffmap26,  MI_FRdiffmap27,  MI_FRdiffmap28,  MI_FRdiffmap29,...
   MI_FRdiffmap30, MI_FRdiffmap31,  MI_FRdiffmap32];


figure;
p = plot(MI_FRdiffmap);
p.Marker = 'square';
xlabel('Eigen vector index'); ylabel('Mutual Information in bits');
title('Mutual info using diffMaps on Gauss-smoothed Spiketimes(FR)');
grid on;

print('Mutual Information on FRdiffusionMaps','-dpng');




%%
%use I(X,Y) = -sum(sum(p(x,y)log(p(x,y)/(p(x)*p(y)))))
% I(X,Y) = -sum(sum(p(x,y)[log(p(x,y)-log(p(x))-log(p(y))]).






