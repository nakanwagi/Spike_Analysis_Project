%% Find the length of the longest spike train to convert
% all logical arrays into matrices
disp('find the length of the longest spike train'); tic;
nrows = zeros(32, 1);
for i = 1:32
    nrows(i) = length(Spikes{i}.T);
end
Nrows = max(nrows);
toc;

%% Store all the real spike data into a matrix/2D array
disp('converting real data from cells to a matrix');tic;
RealData = nan(Nrows, length(Spikes));
for i = 1:length(Spikes)
    RealData(1:length(Spikes{i}.T), i) =  Spikes{i}.T;
end
toc;

%%
disp('initializing tstart , dt = 0.2 [s]');tic;

% Limit spiking to time on and off track
tstart = sd.ExpKeys.TimeOnTrack; %start time in [s] of the real experiment
tend = sd.ExpKeys.TimeOffTrack; %end time in [s] of the real experiment
ncols = size(RealData, 2); %number of neurons or columns in the data

MaxTime = (tend - tstart); %Duration of the experiment.
dt = 0.2; %time step in [s];

%disp(min(unique(time) == time));
%-1*ones(4,5);
%data(1,1:length(x{2}))) = x{2};


%initialize a time vector to index rows of prevtime/nextime
time = tstart:dt:tend-dt;

%caution; the rows of real data increase by one due to two inserted spikes

idx = zeros(size(RealData, 1), ncols); %matrix for prevtime logicals
idx2 = zeros(size(RealData, 1), ncols); %matrix for next time logicals
prevtime = zeros(length(time), ncols); %initialize the prevtime vector
nextime = nan(length(time), ncols);%initialize the nextime vector
tmp = zeros(length(time), ncols); % store indices of idx in tmp matrix
tmp2 = nan(length(time), ncols); % store indices of idx2 in tmp2 matrixx
toc;

%%

disp('Recomputing prevtime/nextime');tic;

for i = 1:length(time)
    for j = 1:ncols
        
        % look for elements </> time respectively
        
        idx(:,j) =  RealData(:,j)  <  time(i);
        idx2(:,j) =  RealData(:,j) >  time(i);
        
        if (sum(idx(:,j))~=0)  % exclude empty arrays
            
            %save the indices corresponding to non-empty arrays from from idx into tmp
            
            count = nnz(~isnan(RealData(:,j))); %count the number of non_NAN values in each column
            tmp =   max(RealData(find(idx(1:count, j)), j)); %choose the maximum number satisfying the idx condition
            
            % compute prevtime
            prevtime(i,j) = tmp-time(i);
            
        else
            prevtime(i,j) = nan;
        end
        
        %repeat the same exact proceedure above for nextime
        if (sum(idx2(:,j))~=0)
            
            count2 = nnz(~isnan(RealData(:,j))); %count the number of non_NAN values in each column
            tmp2 =   min(RealData(find(idx2(1:count2, j)), j)); %choose the minimum number satisfying the idx2 condition
            
            %compute nextime
            nextime(i,j) = tmp2-time(i);
            
        else
            nextime(i,j) = nan;
        end
        
        
        
    end
end
toc;