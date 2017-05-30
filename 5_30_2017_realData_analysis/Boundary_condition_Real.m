function NewSpikes = Boundary_condition_Real(Spikes, Dt, tstart, tend, time)
% BoundaryCondition_NANs  generates a new set of spike times
% called NewSpikes based on SpikeTimes generated by SimHeterogeniousSpikes.

% NewSpikes is a one-by-Ncols cell array  generated 
% by inserting one spike at the beginning and end of the experiment.
% The inserted spikes act as  a boundary condition used to remove NAN
% values (i.e. clean the data obtained) from the previous time (prevtime)
% and  next time functions (nextime).
% Prevtime and nextime functions are designed to overcome the
% coverage problem encoutered while using the firing rate  to do
% dimensionality reduction

% Usage: NewSpikes = BoundaryCondition_NANs(Spikes)

% Input arguments: 
% Spikes-- one-by-ncols cell array obtained from the RedishLab
%  Dt -- time step or bin width in [s] 
%  time -- length(time)-by-one vector describing
% the time stamps for each row of the data
% length(time) = total number of time points.

%Output arguments
% NewSpikes -- one-by-Ncols cell array with an extra spike time inserted
% at the beginning and end of the experiment.
% each cell in NewSpikes has 2 more spike times than generated SpikeTimes

%Default parameters Dt and time are the same as in SimHeterogeneousSpikes.m

%% 
%if nargin == 1
%    Dt = 0.1; % time step in [s]
%    % Limit spiking to time on and off track
%    tstart = sd.ExpKeys.TimeOnTrack; %start time in [s] of the real experiment
%    tend = sd.ExpKeys.TimeOffTrack; %end time in [s] of the real experiment
%    time = tstart:Dt:tend-Dt; %initialize a time vector to index rows of prevtime/nextime
%    %initialize a time vector to index rows of prevtime/nextime;  
%end


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
disp('initializing variables, dt in [s]');tic;


ncols = size(RealData, 2); %number of neurons or columns in the data

%MaxTime = (tend - tstart); %Duration of the experiment.
% dt = 0.15; %MaxTime/Nrows; %time step in [s];

%disp(min(unique(time) == time));
%-1*ones(4,5);
%data(1,1:length(x{2}))) = x{2};


idx = zeros(size(RealData,1), ncols); %matrix for prevtime logicals
idx2 = zeros(size(RealData,1), ncols); %matrix for next time logicals
prevtime = zeros(length(time), ncols); %initialize the prevtime vector
nextime = nan(length(time), ncols);%initialize the nextime vector
%tmp = zeros(length(time), ncols); % store indices of idx in tmp matrix
%tmp2 = nan(length(time), ncols); % store indices of idx2 in tmp2 matrix
toc;

%%

disp('computing prevtime/nextime');tic;

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

%%

disp('removing nans from prevtime/nextime');tic;
%% Compute the mean of PrevTime and NextTime along the time dimension
% and add it to the spiketime vector cell wise.

prevtime(isnan(prevtime)) = []; %delete all nan values from PrevTime---indeed the sizes change!
nextime(isnan(nextime)) = []; %delete all nan values from NextTime---indeed the sizes change!
% before i was using the inbuilt function omitnan which is even too slow.

mean_PrevTime = mean(prevtime, 2);  %compute mean of PrevTime without nans
mean_NextTime = mean(nextime, 2);  %compute mean of NextTime without nans

% Now add a spike before and after spiketime into each individual cell
% do this by adding SpikeTimes - mean(PrevTime) at the beginning of each
% cell array and then statetime(end) + mean(NextTime) at the end of each cell
% array
NewSpikes = cell(1,ncols);

for j = 1:ncols
    NewSpikes{j}.T = [time(1)+mean_PrevTime  Spikes{j}.T' time(end)+mean_NextTime]; 
    %insert a spike at the beginning to remove nans from prevtime.
end

toc;
%% Now, re-calculate the PrevTime and NextTime as before;







              