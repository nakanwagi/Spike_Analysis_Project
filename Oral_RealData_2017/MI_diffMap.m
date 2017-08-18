% Script by Mary Sylvia Agwang

%% First Convert the eigenvectors to tsd objects stored in cells
PrevDiffmaPeigVec = cell(1, size(diffmap,2));
for i = 1:size(diffmap, 2) % i runs over the columns of mappedX2
    PrevDiffmaPeigVec{i} = tsd(time, diffmap(:,i));
end


% new_angle = atan2(mappedX2(:,1), mappedX2(:,2));

%% plot the rat position using only the top eigenvector
figure;
plot3(x.data, y.data, x.data(x.range), 'b');
hold on;
plot3(x.data, y.data,  PrevDiffmaPeigVec{1}.data(x.range), 'r');
xlabel('x');ylabel('y');zlabel('First Eigen Vector');
title('Rat position in blue, first eigenvector only in  red: 2011Data on prevDiffmap');

% print the plot
print('Rat Position and EigVector 1 and 5 on prevDiffmap','-dpng');



xb = linspace(min(x.data), max(x.data), 2*size(diffmap,2)); %range for x
yb = linspace(min(y.data), max(y.data), 2*size(diffmap,2)); %range for y





%% estimate joint histogram by accumulating ,means of using eigvec 1 over rat position

jointHist1 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{1}.data(x.range), 'FUN', @mean); 

jointHist2 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{2}.data(x.range), 'FUN', @mean); 

jointHist3 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{3}.data(x.range), 'FUN', @mean);

jointHist4 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{4}.data(x.range), 'FUN', @mean);

jointHist5 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{5}.data(x.range), 'FUN', @mean); 

jointHist6 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{6}.data(x.range), 'FUN', @mean); 


jointHist7 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{7}.data(x.range), 'FUN', @mean); 

jointHist8 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{8}.data(x.range), 'FUN', @mean); 


jointHist9 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{9}.data(x.range), 'FUN', @mean); 



jointHist10 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{10}.data(x.range), 'FUN', @mean); 

jointHist11 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{11}.data(x.range), 'FUN', @mean); 

jointHist12 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{12}.data(x.range), 'FUN', @mean); 


jointHist13 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{13}.data(x.range), 'FUN', @mean); 

jointHist14 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{14}.data(x.range), 'FUN', @mean); 

jointHist15 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{15}.data(x.range), 'FUN', @mean); 

jointHist16 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{16}.data(x.range), 'FUN', @mean); 

jointHist17 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{17}.data(x.range), 'FUN', @mean); 


jointHist18 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{18}.data(x.range), 'FUN', @mean); 


jointHist19 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{19}.data(x.range), 'FUN', @mean); 

jointHist20 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{20}.data(x.range), 'FUN', @mean); 

jointHist21 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{21}.data(x.range), 'FUN', @mean); 

jointHist22 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{22}.data(x.range), 'FUN', @mean); 

jointHist23 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{23}.data(x.range), 'FUN', @mean); 

jointHist24 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{24}.data(x.range), 'FUN', @mean); 

jointHist25 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{25}.data(x.range), 'FUN', @mean); 

jointHist26 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{26}.data(x.range), 'FUN', @mean); 

jointHist27 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{27}.data(x.range), 'FUN', @mean); 

jointHist28 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{28}.data(x.range), 'FUN', @mean); 

jointHist29 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{29}.data(x.range), 'FUN', @mean); 

jointHist30 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{30}.data(x.range), 'FUN', @mean); 

jointHist31 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{31}.data(x.range), 'FUN', @mean); 

jointHist32 = histcn([x.data y.data], xb, yb,...,
    'AccumData', PrevDiffmaPeigVec{32}.data(x.range), 'FUN', @mean); 



%% plot the histogram without zero values
figure;
%set zero values in jointhistogram to nan so log2 is defined
b1 = subplot(2, 2, 1);
jointHist1(jointHist1==0)= nan; 
pcolor(jointHist1);
colorbar;
title(b1, 'Position and Eigev 1 on prevDiffmap');

b2 = subplot(2, 2, 2);
jointHist2(jointHist2==0)= nan; 
pcolor(jointHist2);
colorbar;
title(b2, 'Position and Eigev 2 on prevDiffmap');

b3 = subplot(2, 2, 3);
jointHist3(jointHist3==0)= nan; 
pcolor(jointHist3);
colorbar;
title(b3, 'Position and Eigev 3 on prevDiffmap');


b4 = subplot(2, 2, 4);
jointHist4(jointHist4==0)= nan; 
pcolor(jointHist4);
colorbar;
title(b4, 'Position and Eigev 4 on prevDiffmap');

print('sublot1 on prevDiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b5 = subplot(2, 2, 1);
jointHist5(jointHist5==0)= nan; 
pcolor(jointHist5);
colorbar;
title(b5, 'Position and Eigev 5 on prevDiffmap');

b6 = subplot(2, 2, 2);
jointHist6(jointHist6==0)= nan; 
pcolor(jointHist6);
colorbar;
title(b6, 'Position and Eigev 6 on prevDiffmap');

b7 = subplot(2, 2, 3);
jointHist7(jointHist7==0)= nan; 
pcolor(jointHist7);
colorbar;
title(b7, 'Position and Eigev 7 on prevDiffmap');


b8 = subplot(2, 2, 4);
jointHist8(jointHist8==0)= nan; 
pcolor(jointHist8);
colorbar;
title(b8, 'Position and Eigev 8 on prevDiffmap');

print('Subplot 2 on prevDiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b9 = subplot(2, 2, 1);
jointHist9(jointHist9==0)= nan; 
pcolor(jointHist9);
colorbar;
title(b9, 'Position and Eigev 9 on prevDiffmap');

b10 = subplot(2, 2, 2);
jointHist10(jointHist10==0)= nan; 
pcolor(jointHist10);
colorbar;
title(b10, 'Position and Eigev 10 on prevDiffmap');

b11 = subplot(2, 2, 3);
jointHist11(jointHist11==0)= nan; 
pcolor(jointHist11);
colorbar;
title(b11, 'Position and Eigev 11 on prevDiffmap');


b12 = subplot(2, 2, 4);
jointHist12(jointHist12==0)= nan; 
pcolor(jointHist12);
colorbar;
title(b12, 'Position and Eigev 12 on prevDiffmap');

print('Subplot 3 on prevDiffmap','-dpng');

figure;
%set zero values in jointhistogram to nan so log2 is defined
b13 = subplot(2, 2, 1);
jointHist13(jointHist13==0)= nan; 
pcolor(jointHist13);
colorbar;
title(b13, 'Position and Eigev 13 on prevDiffmap');

b14 = subplot(2, 2, 2);
jointHist14(jointHist14==0)= nan; 
pcolor(jointHist14);
colorbar;
title(b14, 'Position and Eigev 14 on prevDiffmap');

b15 = subplot(2, 2, 3);
jointHist15(jointHist15==0)= nan; 
pcolor(jointHist15);
colorbar;
title(b15, 'Position and Eigev 15 on prevDiffmap');


b16 = subplot(2, 2, 4);
jointHist16(jointHist16==0)= nan; 
pcolor(jointHist16);
colorbar;
title(b16, 'Position and Eigev 16 on prevDiffmap');

print('Subplot 4 on prevDiffmap','-dpng');



figure;
%set zero values in jointhistogram to nan so log2 is defined
b17 = subplot(2, 2, 1);
jointHist17(jointHist17==0)= nan; 
pcolor(jointHist17);
colorbar;
title(b17, 'Position and Eigev 17 on prevDiffmap');

b18 = subplot(2, 2, 2);
jointHist18(jointHist18==0)= nan; 
pcolor(jointHist18);
colorbar;
title(b18, 'Position and Eigev 18 on prevDiffmap');

b19 = subplot(2, 2, 3);
jointHist19(jointHist19==0)= nan; 
pcolor(jointHist19);
colorbar;
title(b19, 'Position and Eigev 19 on prevDiffmap');


b20 = subplot(2, 2, 4);
jointHist20(jointHist20==0)= nan; 
pcolor(jointHist20);
colorbar;
title(b20, 'Position and Eigev 20 on prevDiffmap');

print('Subplot 5 on prevDiffmap','-dpng');

figure;
%set zero values in jointhistogram to nan so log2 is defined
b21 = subplot(2, 2, 1);
jointHist21(jointHist21==0)= nan; 
pcolor(jointHist21);
colorbar;
title(b21, 'Position and Eigev 21 on prevDiffmap');

b22 = subplot(2, 2, 2);
jointHist22(jointHist22==0)= nan; 
pcolor(jointHist22);
colorbar;
title(b22, 'Position and Eigev 22 on prevDiffmap');

b23 = subplot(2, 2, 3);
jointHist23(jointHist23==0)= nan; 
pcolor(jointHist23);
colorbar;
title(b23, 'Position and Eigev 23 on prevDiffmap');


b24 = subplot(2, 2, 4);
jointHist24(jointHist24==0)= nan; 
pcolor(jointHist24);
colorbar;
title(b24, 'Position and Eigev 24 on prevDiffmap');

print('Subplot 6 on prevDiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b25 = subplot(2, 2, 1);
jointHist25(jointHist25==0)= nan; 
pcolor(jointHist25);
colorbar;
title(b25, 'Position and Eigev 25 on prevDiffmap');

b26 = subplot(2, 2, 2);
jointHist26(jointHist26==0)= nan; 
pcolor(jointHist26);
colorbar;
title(b26, 'Position and Eigev 26 on prevDiffmap');

b27 = subplot(2, 2, 3);
jointHist27(jointHist27==0)= nan; 
pcolor(jointHist27);
colorbar;
title(b27, 'Position and Eigev 27 on prevDiffmap');


b28 = subplot(2, 2, 4);
jointHist28(jointHist28==0)= nan; 
pcolor(jointHist28);
colorbar;
title(b28, 'Position and Eigev 28 on prevDiffmap');

print('Subplot 7 on prevDiffmap','-dpng');


figure;
%set zero values in jointhistogram to nan so log2 is defined
b29 = subplot(2, 2, 1);
jointHist29(jointHist29==0)= nan; 
pcolor(jointHist29);
colorbar;
title(b29, 'Position and Eigev 29 on prevDiffmap');


b30 = subplot(2, 2, 2);
jointHist30(jointHist30==0)= nan; 
pcolor(jointHist30);
colorbar;
title(b30, 'Position and Eigev 30 on prevDiffmap');

b31 = subplot(2, 2, 3);
jointHist31(jointHist31==0)= nan; 
pcolor(jointHist31);
colorbar;
title(b31, 'Position and Eigev 31 on prevDiffmap');


b32 = subplot(2, 2, 4);
jointHist32(jointHist32==0)= nan; 
pcolor(jointHist32);
colorbar;
title(b32, 'Position and Eigev 32  on prevDiffmap');

print('Subplot 8 prevDiffmap','-dpng');

% %% Recover the Rat position using the first eigenvector only.
% 
% xb = linspace(min(x.data), max(x.data), 2*size(mappedX2,2)); %range for x
% yb = linspace(min(y.data), max(y.data), 2*size(mappedX2,2)); %range for y




%% Compute the Mutual_Information for each of the 32 eigenvectors

iV = 32;
% I have only put 10 because it's easy to see the contributions of
% each eigenvector on submatrix plot
eigVecPrevDiffmap = linspace(min(PrevDiffmaPeigVec{iV}.data), max(PrevDiffmaPeigVec{iV}.data), length(x.range));
%compute a 3D hsitogram to estimate the joint entropy
H = histcn([x.data, y.data, PrevDiffmaPeigVec{iV}.data(x.range)],xb, yb, eigVecPrevDiffmap);
% compute the normalizing constant;

% Matricize the tensor to speed up MATLAB computations

Hnew = reshape(H, [], size(H,3)); %convert to 64^2*length(x.data)
% remove zeros so that log2 is defined
Hnew(Hnew==0)=nan;

Hnew = Hnew./nansum(nansum(nansum(Hnew)));

lh1 = log2(nansum(Hnew,1));
lh2 = log2(nansum(Hnew,2));

MI_prevDiffmap32 = nansum(nansum(Hnew .* bsxfun(@minus,bsxfun(@minus,log2(Hnew),lh1),lh2)));


%% Create a Mutual Information vector.

MI_prevDiffmap = [MI_prevDiffmap1, MI_prevDiffmap2,  MI_prevDiffmap3,  MI_prevDiffmap4,...
    MI_prevDiffmap5, MI_prevDiffmap6,  MI_prevDiffmap7, MI_prevDiffmap8, MI_prevDiffmap9,...
    MI_prevDiffmap10, MI_prevDiffmap11, MI_prevDiffmap12, MI_prevDiffmap13, MI_prevDiffmap14,...
    MI_prevDiffmap15,  MI_prevDiffmap16,  MI_prevDiffmap17, MI_prevDiffmap18,  MI_prevDiffmap19,...
    MI_prevDiffmap20, MI_prevDiffmap21,  MI_prevDiffmap22, MI_prevDiffmap23,  MI_prevDiffmap24,...
   MI_prevDiffmap25, MI_prevDiffmap26,  MI_prevDiffmap27,  MI_prevDiffmap28,  MI_prevDiffmap29,...
   MI_prevDiffmap30, MI_prevDiffmap31,  MI_prevDiffmap32];


figure;
p = plot(MI_prevDiffmap);
p.Marker = 'square';
xlabel('Eigen vector index'); ylabel('Mutual Information in bits');
title('Mutual info Vs idx using DiffMaps on Previous time');
grid on;

print('Mutual Information on prevDiffMaps','-dpng');

%%
% %use I(X,Y) = -sum(sum(p(x,y)log(p(x,y)/(p(x)*p(y)))))
% % I(X,Y) = -sum(sum(p(x,y)[log(p(x,y)-log(p(x))-log(p(y))]).
% 
% 
% 



