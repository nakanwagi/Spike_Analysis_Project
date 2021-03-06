%% Clear workspace, command window and close figures

close all;
clc;

%% Redish Lab spike data is sd.S

disp('reading raw data');tic;

%2011 data
  sd = TaskInit('/home/agwan003/Sylvia_Analysis/David_Data/R232-2011-10-16'); %load the spike data ----2011 is what is in my workspace

%2012 data
% sd = TaskInit('/home/agwan003/Sylvia_Analysis/David_Data/R236-2012-12-24'); %load the spike data ----2011 is what is in my workspace'); %load the spike data ----2011 is what is in my workspace

%2010  data
%sd = TaskInit('/home/agwan003/Sylvia_Analysis/David_Data/R199-2010-06-11'); %load the spike data ----2011 is what is in my workspace



Spikes = sd.S; %get spikes cells from the struct sd

toc;

%% Plot the Spike Data to visualize
%  Q = MakeQfromS(Spikes, 1); %bin width dt = 0.5  %0.1 ms
 DT = 0.1;
 sigma = 1;
 
 Q = MakeQfromS_spreadSpikes(Spikes, DT, sigma); %bin width dt = 0.5  %0.1 ms
 
 %Q = MakeQfromS(Spikes, 1); %bin width dt = 0.5  %0.1 ms
 
 Data = Q.D;
 
 Data = full(Data);
 
 %visualize the data
 figure(1);
 imagesc(Q.range, 1:32, Q.data');
 %title('image of firing rates data');
 xlabel('firing rate');ylabel('neuron label');
 print('Firing Rates data','-dpng');
 
 %plot the data with the bursting neuron 8 in red.
 figure(2);
  plot(sd.x.data, sd.y.data, '.', sd.x.data(sd.S{12}), sd.y.data(sd.S{12}), 'ro');
 %title('bursting neurons in red and the rest of the data');
 print('Bursting neuron and the rest of the data','-dpng');
 
 
%  tsd(Q); %gives the data with rows as spike times and columns as neurons
%  
%  sd.S{8}; %refers to the 8th spike (T=5361 means this neuron spikes this frequently)

%% Initialize the input variables.
 Dt = 0.1; % time step in [s]
   % Limit spiking to time on and off track
   tstart = sd.ExpKeys.TimeOnTrack; %start time in [s] of the real experiment
   tend = sd.ExpKeys.TimeOffTrack; %end time in [s] of the real experiment
   time = tstart:Dt:tend-Dt; %initialize a time vector to index rows of prevtime/nextime
   %initialize a time vector to index rows of prevtime/nextime;


% %%  Insert a boundary condition to remove NaN values
% NewSpikes = Boundary_condition_Real(Spikes, Dt, tstart, tend, time);
% 
% 
% %%  Recompute prevtime and nextime without NaN values.
% [prevtime, nextime] = Recompute_prevNext(NewSpikes, Dt, tstart,tend,time);
% 
%  
% 
% %% Generate the distance matrix using the exponential kernel
% 
% disp('run Diffusion Maps Algorithm');tic;
% 
% tau = 12; %-----with noise; %0.9----no noise; 
% 
% 
% X2 = (1/tau)*exp(-abs(prevtime)/tau);  
% 
% 
% 
% %% call diffusion maps algorithm on [exp-(abs(pretime)) exp-(abs(nexttime))]
% 
% %[mappedX1, mapping1] = compute_mapping(X1, 'DiffusionMaps');
% 
% % call diffusion maps algorithm on [-exp(prevtime)] only
% 
% sig_d = 1; %variance of gaussian in diffusion maps (default = 1)
% 
% alpha = 1; %parameter in diffusion maps (default = 1)
% 
% no_dims = length(Spikes); %number of dimensions. (default = 2)
% 
% %% Run diffusion Maps on  prevtime/nextime data
%  
% [mappedX2, mapping2] = compute_mapping_sylvia(X2, 'DiffusionMaps', no_dims, alpha, sig_d);
% 
%  
% 
% %% Run diffusion Maps on Firing rate data
% 
% [fireRate1, fireRate_struct] = compute_mapping_sylvia(Data, 'DiffusionMaps', no_dims, alpha, sig_d);
% 
% toc;
% 
% 
% 
% 







