

%% Clear workspace, command window and close figures

close all;
clc;

%% Add TaskInit directory
%% cd /home/agwan003/Sylvia_Analysis
addpath(pwd)

disp('reading raw data');
tic;
%sd = TaskInit('/home/tesylvia/Documents/MATLAB_working_dir/David_Hippocampal_Data/R232-2011-10-16'); %load the spike data
sd = TaskInit('/home/agwan003/Sylvia_Analysis/David_Data/R232-2011-10-16'); %load the spike data
Spikes = sd.S; %get spikes cells from the struct sd
toc;

% %% Plot the Spike Data to visualize
% %  Q = MakeQfromS(Spikes, 1); %bin width dt = 0.5  %0.1 ms
%  DT = 0.05;
%  sigma = 1;
%  
%  Q = MakeQfromS_spreadSpikes(Spikes, DT, sigma); %bin width dt = 0.5  %0.1 ms
%  
%  %Q = MakeQfromS(Spikes, 1); %bin width dt = 0.5  %0.1 ms
%  
%  Data = Q.D;
%  
%  Data = full(Data);
%  
%  %visualize the data
%  figure(1);
%  imagesc(Q.range, 1:32, Q.data');
%  title('image of firing rates data');
%  
%  %plot the data with the bursting neuron 8 in red.
%  figure(2);
%   plot(sd.x.data, sd.y.data, '.', sd.x.data(sd.S{12}), sd.y.data(sd.S{12}), 'ro');
%  title('bursting neurons in red and the rest of the data');
%  
% %  tsd(Q); %gives the data with rows as spike times and columns as neurons
% %  
% %  sd.S{8}; %refers to the 8th spike (T=5361 means this neuron spikes this frequently)
% 

%% Initialize the main variables for analysis.

disp('initializing tsart_2660_2750, dt = 0.02 [s]');tic;

% Limit spiking to time on and off track
tstart = sd.ExpKeys.TimeOnTrack; %start time in [s] of the real experiment
tend = sd.ExpKeys.TimeOffTrack; %end time in [s] of the real experiment
ncols = length(Spikes); %number of neurons or columns in the data
dt = 0.002; % timestep in [s] (due to memory problem of diffusion maps)
Start = 2660; %start time based on frequency of firing in [s]
Stop = 2750; % end time based on frequency of firing in [s]
maxTime = Stop - Start; % Duration based on frequency of firing in [s]
nrows = ceil(maxTime/dt); %number of rows in data matrix based on time vec
dt = maxTime/nrows; %redefine dt based on nrows
time = Start:dt:Stop-dt;%time vector in [s] indexing the rows of the matrix


%% =========================================================================
%initialize prevtime and nextime
nextime = nan(length(time), ncols);
prevtime = nan(length(time), ncols);
idx = cell(length(time), ncols);
idx2 = cell(length(time), ncols);
tmp = cell(length(time), ncols);
tmp2 = cell(length(time), ncols);

toc;
 
 %% Comoute prevtime and nextime for the 32 spike trains.
 
 disp('computing prevtime and nextime');tic;
 
 
 for i = 1:nrows
     for j = 1:ncols
          idx{i,j} = sd.S{j}.T < time(i);  % get logicals from the jth spike train corresponding to i^th time point.
          idx2{i,j} = sd.S{j}.T > time(i);  %get logicals as above for nexttime.
          
          if sum(idx{i,j})==0  %set prevtime = nan if idx is empty
              prevtime(i,j) = nan;
          else
              tmp{i,j} = sd.S{j}.T(idx{i,j}); %compute prevtime if idx is not empty
          prevtime(i,j) = tmp{i,j}(end) - time(i);
          
          end
          
         if sum(idx2{i,j})==0 %set nextime = nan if idx2 is empty
            nextime(i,j) = nan;
          else
              tmp2{i,j} = sd.S{j}.T(idx2{i,j}); %compute prevtime if idx is not empty
              nextime(i,j) = tmp2{i,j}(1) - time(i);
          
          end        
          
     end
 end

 toc;
 
 
%% Compute the mean of PrevTime and NextTime along the time dimension
% and add it to the spiketime vector cell wise.

disp('remove NANs from prevtime by inserting spike at tstart');tic;


prevtime(isnan(prevtime)) = []; %delete all nan values from PrevTime---indeed the sizes change!
nextime(isnan(nextime)) = []; %delete all nan values from NextTime---indeed the sizes change!
% before i was using the inbuilt function omitnan which is even too slow.

mean_PrevTime = mean(prevtime, 2);  %compute mean of PrevTime without nans
mean_NextTime = mean(nextime, 2);  %compute mean of NextTime without nans

% Now add a spike before and after spiketime into each individual cell
% do this by adding SpikeTimes - mean(PrevTime) at the beginning of each
% cell array and then statetime(end) + mean(NextTime) at the end of each cell
% array
for j = 1:ncols
    sd.S{j}.T = [mean_PrevTime  sd.S{j}.T']; %insert a spike at the beginning to remove nans from prevtime.
end

toc;
%% Now, re-calculate the PrevTime and NextTime as before;
