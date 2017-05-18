%%Matlab Script by Mary Sylvia Agwang


disp('computing firing rate'); tic;
%% Clear the command window and close all open figure


clc;


%% Facts about the implementation of Laplacian eigen Maps
% default settings in order for varargin
% k = number of nearest neigbors in a neighborhood graph and default = 12.
% t = 2*sigma^2 where default sigma = 1.
%'MATLAB' is the default setting when data points are less than 10,000.
% mappedX are no_dims eigenvectors e.g 2 eigenvectors if no_dims = 2;

%Fact: default noise or std of the noise added to Laplacian synthetic data is 0.05.

%% Define the required parameters
Duration = 6.5; %maximum time taken to carry out the experiment in [s]
T_lap  = 1.5; % time taken to make one lap around the circle in (period) [s]
number_of_laps = Duration/T_lap;
%dist = 2*pi--- the distance around a circle of radius 1
speed = 2*pi/T_lap; %speed of the rat (~~ 4.2 rads/sec)
% angle  = speed*time;

%%Question:   ========= Is it possible for dt in FireRate to be different 
%%============= from the DT in isSpike?==================================

X2 =  bsxfun(@times, heaviside(time)',((1/tau)*exp(-abs(prevtime)/tau)));
dt = 0.0015;         %0.00001;   % in [s] 
%divide the time vector into very small intervals in [sec]
% so that the poisson process is approximately bernouli since dt<<<1.
% dt is the step size or bin width depending on the purpose.


nrows = Duration/dt; % number of time points 
ncols = 32; % number of neurons or columns of the data
N = ncols;
%c = 2*pi; %speed of the rat
t = 0:dt:Duration-dt;  % tvec is a vector describing the time stamps for each row of the data
% i.e. length(tvec) = Timepoints;----useful for making plots.
theta =  speed*t; %*t; %c*t; % position of the rat at time t



 %time(find(time == 0.0001)) = [];



%==no heterogeneity in sigma, all receptive fields have the same width=====
%sig = 0.07;
sig = 1; %perfect case
%sig = 0.01; %bad case

%Case 1:
% when sigma is the same through out 


 f=@(x, maxR) maxR*exp(-((mod(x+pi, 2*pi)-pi)).^2./(2.*sig.^2));


%=================== add heterogeneity in sigma so that other receptive
%===================fields are wider than others=========================

 % sig_max = 1; %0.1;
 % sig_min = 0.01;
 % sigma = sig_min + rand(1, N)*(sig_max - sig_min);
%=======================================================================

% we don't want sig_min to be smaller than the smallest
%sigma for which the algorithm breaks down


%=======================================================
% Case 2:
% when sigma is heterogeneous.
%f=@(x, sig) exp(-((mod(x+pi, 2*pi)-pi)).^2./(2.*sig.^2));



%==========================================================

%% Using uniform sampling
% for each rat position, all the 32 neurons fire according to the
% function FireRate.
%==========================================================================
%case 1: center = (2*pi*i)/N %centers evenly spread out in space.

%=========================================================================
%case 2: centers uniformly distributed in space.
%center = rand(1, N)*2*pi; % centers must be the same size as N

%========================================================================
FireRate = zeros(length(t), ncols);
maxRate = 20; %the maximum firing rate helps to ensure that we have many spiketimes for each neuron
% the larger the max rate, the more spike times we get in each cell
% be aware that some cells may not spike in which case we get an empty
% spike  time matrix.
for i = 1:ncols
    
    % heterogeneous sigma
    %FireRate(:,i) = f(theta - center(i), sigma(i)); %Uniform sampling of neuron positions
    
    % same sigma for all receptive fields
    FireRate(:,i) = f(theta - (2*pi*i)/ncols, maxRate); %Uniform sampling of neuron positions
end

%==========================================================================

toc;

%%

disp('generating spikes from firing rate'); tic;



time = 0:dt:Duration-dt;
%check if uniform rand number in (0, 1) is less than prob of 1 spike
% in interval Dt and store any spikes in the matrix isSpike.
isSpike = rand(length(time), ncols) < FireRate*dt;

%extract the time when isSpikes had the value 1 and
%store it in the vector SpikeTime(i) for the i^th cell/neuron.

SpikeTimes = cell(1, ncols);
for ii = 1:ncols
   SpikeTimes{ii} = time(isSpike(:,ii)); %SpikeTimes contains the spike times for the i^th neuron
   
end

toc;
%========Chose state times uniformly in space=============================

%================changed this part yesterday 3/13/2017 ========================================%

%  mean_FireRate = 1./mean(FireRate); %how do I make sure I have the right units?
% 
% statetime = time;
% 
% extra_spike = mean_FireRate + statetime(end);
%  SpikeTimes = [num2cell(-mean_FireRate) SpikeTime num2cell(extra_spike)];


%===========end change of 3/13/2017============================================%

% %Compute the time of previous spike of i^{th} neuron and store it
% % in the vector called PrevTime

%=======Initialize with nan values=========================================
%PrevTime = zeros(length(time), nCols);  % loop over all brain states

%=======Initialize with nan values=========================================

disp('computing prevtime/nextime'); tic;

prevtime = nan(length(time), ncols);  % loop over all brain states
%Prevtime = (time of previous spike before time t) - t

%=========================================================================
%PrevTime = nan(1, nCols);  % to check what happens at a single brain state
%=========================================================================


%time since the last spike (encodes the past of the neural activity)
% Based on PrevTime, the neuron which fired most recently given a specific
% statetime is one which fired abs(x) seconds ago where abs(x) satisfies
% these conditions: abs(x) < statetime and has 
% smallest |statetime - abs(x)|.
 
%=============initialize with zeros=====================================
%NextTime = zeros(length(time), nCols); % to check what happens at a single brain state

%================initialize with nan values=============================
nextime = nan(length(time), ncols); % to check what happens at a single brain state
%NextTime = (time of next spike after time t) - t.

%=========================================================================X2 =  bsxfun(@times, heaviside(time)',((1/tau)*exp(-abs(prevtime)/tau)));
%NextTime = nan(1, nCols);  %loop over all brain states
%=========================================================================


%time after the last spike (encodes the past of the neural activity)
%Based on PrevTime, the neuron for which it's been a very long while since
%it last fired is the one which has x seconds where x  statisfies two
%conditions: x > statetime and has the smallest |x - statetime|.


%========================================================================
% Uncomment this if you want to see only one brain state
%nn = find(statetime==0);  
%nn = find(statetime==2.2); %at this values all neurons are active


%=========================================================================
idx = cell(length(time), ncols);
 idx2 = cell(length(time), ncols);
 tmp = cell(length(time), ncols);
 tmp2 = cell(length(time), ncols);
 
 %%
 for ll = 1:nrows
     for kk = 1:ncols
          idx{ll,kk} = SpikeTimes{kk} < time(ll);  % get logicals from the jth spike train corresponding to i^th time point.
          idx2{ll,kk} = SpikeTimes{kk} > time(ll);  %get logicals as above for nexttime.
          
          if sum(idx{ll,kk})==0
              prevtime(ll,kk) = nan;
          else
              tmp{ll,kk} = SpikeTimes{kk}(idx{ll,kk});
          prevtime(ll,kk) = tmp{ll,kk}(end) - time(ll);
          
          end
          
         if sum(idx2{ll,kk})==0
            nextime(ll,kk) = nan;
          else
              tmp2{ll,kk} = SpikeTimes{kk}(idx2{ll,kk});
              nextime(ll,kk) = tmp2{ll,kk}(1) - time(ll);
          
          end        
          
     end
 end

 
 toc;
%% Compute the mean of PrevTime and NextTime along the time dimension
% and add it to the spiketime vector cell wise.

disp('removing NANs from prevtime'); tic;

prevtime(isnan(prevtime)) = []; %delete all nan values from PrevTime---indeed the sizes change!
nextime(isnan(nextime)) = []; %delete all nan values from NextTime---indeed the sizes change!
% before i was using the inbuilt function omitnan which is even too slow.

mean_PrevTime = mean(prevtime, 2);  %compute mean of PrevTime without nans
mean_NextTime = mean(nextime, 2);  %compute mean of NextTime without nans

%% Now add a spike before and after spiketime into each individual cell
% do this by adding SpikeTimes - mean(PrevTime) at the beginning of each
% cell array and then statetime(end) + mean(NextTime) at the end of each cell
% array
for mm = 1:numel(SpikeTimes)
    SpikeTimes{mm} = [mean_PrevTime  SpikeTimes{mm}  time(end)+mean_NextTime];
end

toc;
%% Now, re-calculate the PrevTime and NextTime as before;
