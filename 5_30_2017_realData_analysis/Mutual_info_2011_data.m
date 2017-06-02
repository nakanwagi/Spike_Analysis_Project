%% First get an estimate of how many dimensions
% are most significant.   
lambda = mapping2.val;
ind = mapping2.index;
    
%% Make a scree plot to guess the most important dimensions
figure;
plot(ind, lambda, 'ko-');
xlabel('index');ylabel('eigenvalue');
title('Scree plot: 2011RedishLabData');

% %% find out how much variance the first k eigenvectors explain
 eigvec_no = 2;
explained_variance = 100*cumsum(lambda(1:eigvec_no))/sum(lambda(1:eigvec_no));
% %var = mean(var)*100; %percentage of variance explained by first k eigenvectors












    
    