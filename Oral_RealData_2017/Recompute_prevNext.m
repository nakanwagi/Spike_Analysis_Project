function [prevtime, nextime] = Recompute_prevNext(NewSpikes, Dt, tstart,tend,time)
% RecomputePrevtimeNextime  generates  clean previous time (prevtime)
% and next time functions (nextime) by  recomputing the prevtime
% nextime functions.

% NewSpikes is a one-by-Ncols cell array  generated 
% by Boundary_condition_Real.m file
% Prevtime and nextime functions are designed to overcome the
% coverage problem encoutered while using the firing rate  to do
% dimensionality reduction

% Usage: [prevtime, nextime] = RecomputePrevtimeNextime(NewSpikes, Dt, tstart, tend, time)

% Input arguments: 
% NewSpikes-- one-by-ncols cell array obtained using the
%  BoundaryCondition_NANs  function.
%  Dt -- time step or bin width in [s] 
%  times -- length(times)-by-one vector describing
% the time stamps for each row of the data
% length(time) = total number of time points.

%Output arguments
% prevtime -- length(times)-by-ncols matrix containing the time
% since the last spike encoding past firing activity of the neurons
% nextime - length(time)-by-ncols matrix  containing the time
% until the next spike encoding future firing activity of the neurons.


% %Default parameters Dt and time are the same as in SimHeterogeneousSpikes.m
% 
% %% 
% if nargin == 1
%    Dt = 0.0015; % time step in [s]
%    tstart = % start time of experiment in [s]
%    tend = % end time of the experiment in [s]
%    times = tstart:Dt:tend-Dt; 
%    %initialize a time vector to index rows of prevtime/nextime;  
% end


%% Find the length of the longest spike train to convert 
% all logical arrays into matrices
disp('find the length of the longest spike train'); tic;
nrows = zeros(length(NewSpikes), 1);
for i = 1:length(NewSpikes)
    nrows(i) = length(NewSpikes{i}.T);
end
Nrows = max(nrows);
toc;

%% Store all the real spike data into a matrix/2D array
disp('converting real data from cells to a matrix');tic;
RealData = nan(Nrows, length(NewSpikes));
for i = 1:length(NewSpikes)
    RealData(1:length(NewSpikes{i}.T), i) =  NewSpikes{i}.T; 
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









              
