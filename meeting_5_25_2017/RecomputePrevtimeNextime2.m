function [prevtime, nextime] = RecomputePrevtimeNextime2(NewSpikes, ~, times, ~, ~)
% RecomputePrevtimeNextime  generates  clean previous time (prevtime)
% and next time functions (nextime) by  recomputing the prevtime
% nextime functions.

% NewSpikes is a one-by-Ncols cell array  generated 
% by BoundaryCondition_NANs m-file
% Prevtime and nextime functions are designed to overcome the
% coverage problem encoutered while using the firing rate  to do
% dimensionality reduction

% Usage: [prevtime, nextime] = RecomputePrevtimeNextime(NewSpikes)

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


%Default parameters Dt and time are the same as in SimHeterogeneousSpikes.m

%% 
if nargin == 1
   Dt = 0.0015; % time step in [s]
   tstart = 0; % start time of experiment in [s]
   tend = 6.5; % end time of the experiment in [s]
   times = tstart:Dt:tend-Dt; 
   %initialize a time vector to index rows of prevtime/nextime;  
end



%% Find the length of the longest spike train to convert 
% all logical arrays into matrices

nrows = zeros(length(NewSpikes), 1);
for i = 1:length(NewSpikes)
    nrows(i) = length(NewSpikes{i});
end
Nrows = max(nrows);


%% Store all the real spike data into a matrix/2D array

SyntheticData = nan(Nrows, length(NewSpikes));
for i = 1:length(NewSpikes)
    SyntheticData(1:length(NewSpikes{i}), i) =  NewSpikes{i}; 
end



%%  Limit spiking to time on and off track

Ncols = size(SyntheticData, 2); %number of neurons or columns in the data

idx = zeros(size(SyntheticData,1), Ncols); %matrix for prevtime logicals
idx2 = zeros(size(SyntheticData,1), Ncols); %matrix for next time logicals
prevtime = zeros(length(times), Ncols); %initialize the prevtime vector
nextime = nan(length(times), Ncols);%initialize the nextime vector
%tmp = zeros(length(times), Ncols); % store indices of idx in tmp matrix
%tmp2 = nan(length(times), Ncols); % store indices of idx2 in tmp2 matrix



%%


for i = 1:length(times)
  for j = 1:Ncols
  
      % look for elements </> time respectively
      
        idx(:,j) =  SyntheticData(:,j)  <  times(i);
        idx2(:,j) =  SyntheticData(:,j) >  times(i);
    
    if (sum(idx(:,j))~=0)  % exclude empty arrays

   %save the indices corresponding to non-empty arrays from from idx into tmp
   count = nnz(~isnan(SyntheticData(:,j))); %count the number of non_NAN values in each column
    tmp =   max(SyntheticData(find(idx(1:count, j)), j)); %choose the maximum number satisfying the idx condition
  
   % compute prevtime 
   prevtime(i,j) = tmp-times(i);
 
        else
          prevtime(i,j) = nan;
    end
    
    %repeat the same exact proceedure above for nextime
    if (sum(idx2(:,j))~=0)  

   count2 = nnz(~isnan(SyntheticData(:,j))); 
   %count the number of non_NAN values in each column
  tmp2 =   min(SyntheticData(find(idx2(1:count2, j)), j)); 
  %choose the minimum number satisfying the idx2 condition
  
  %compute nextime
  nextime(i,j) = tmp2-times(i);
 
        else
          nextime(i,j) = nan;
    end
    
    
    
  end
end
