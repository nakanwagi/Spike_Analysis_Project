
%clear all;
close all;
%clc;


%%

B = cell(5,1);
B{1} = [2;3;3];
B{2} = [4;5;6;7];
B{3} = [3;4;6;8;10];
B{4} = [1;2;6;8;10;12];
B{5} = [6;8;10];


%%
disp('computing length of longest spike train');tic;
nrows = zeros(length(B),1);

for i = 1:length(B)
    nrows(i) = length(B{i});   %find the length of each spike train
end
Nrows = max(nrows); %get the length of longest spike train
toc;

%% Convert from cell to matrix array

disp('converting real data from cells to a matrix');tic;


Data = nan(Nrows, length(B)); %initialize the real data array
ncols = length(B);
for i = 1:ncols
    Data(1:length(B{i}), i) = B{i}; %imitating my real data set up
end

toc;

disp('computing prevtime/nextime');tic;

dt =2;
time = 0:dt:10;

%%  Initialize prevtime and nextime matrices
idx = zeros(length(time), ncols);
idx2 = zeros(length(time), ncols);
tmp =  zeros(length(time), ncols);
tmp2 = nan(length(time), ncols);
prevtime = zeros(length(time), ncols);
nextime = nan(length(time), ncols);

for i = 1:Nrows
  for j = 1:ncols
  
      % look for elements </> time respectively
      
        idx(:,j) =  Data(:,j)  <  time(i);
        idx2(:,j) =  Data(:,j) >  time(i);
    
    if (sum(idx(:,j))~=0)  % exclude empty arrays

   %save the indices corresponding to non-empty arrays from from idx into tmp     
   
   tmp = max(Data(find(idx(1:nnz(~isnan(Data(:,j))),j)),j));
   % compute prevtime 
   prevtime(i,j) = tmp-time(i);%max(tmp(:,j)) - time(i);
 
        else
          prevtime(i,j) = nan;
    end
    
    %repeat the same exact proceedure above for nextime
    if (sum(idx2(:,j))~=0)  
        
   tmp2=min(Data(find(idx2(1:nnz(~isnan(Data(:,j))),j)),j));
   nextime(i,j) = tmp2-time(i);%min(tmp2(:,j)) - time(i);
 
        else
          nextime(i,j) = nan;
    end
    
    
  end
end
toc;
