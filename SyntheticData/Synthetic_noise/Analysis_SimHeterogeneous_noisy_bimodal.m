%% Script by Mary Sylvia Agwang to analyse synthetic data
% Oral Project: Summer 2017


clc;

%% Fix the random number generator for MATLAB
% ---why does change in the way I wrote the function cause an error????
% function SpikeTimes = SimHeterogeneousSpikes(tstart,...,
%     tend,dt, theta, time, ~, ~, ncols,sigma)


%% Define the default input parameters.

ncols = 32; % number of neurons or columns of the data
dt = 0.0015;  % time step or bin width in [s] 
% divide the time vector into very small intervals in [s]
% so that the poisson process is approximately bernouli since dt<<<1.
% dt is the step size or bin width depending on the purpose.
tend = 6.5; % end time of the experiment in [s]
tstart = 0; % start time  of the experiment in [s]
Duration = tend - tstart; % Maximum time of the experiment in [s]
% nrows = round(Duration/dt); % total number of time points/rows in the 

%% Default parameters are given below

% ncols = 32; % number of neurons or columns of the data
 T_lap  = 1.5; % period or time to make one lap in [s] around circle.
%  sig_max = 1; 
%  sig_min = 0.01;
%  sigma = sig_min + rand(1, ncols)*(sig_max - sig_min);
% % we don't want sig_min to be smaller than the smallest
% %sigma for which Laplacian eigenMaps breaks down.

time = tstart:dt:tend-dt;  % length(time)-by-one vector describing
% the time stamps for each row of the data
% length(time) = total number of time points;----useful for making plots.

speed = 2*pi/T_lap; %speed of the rat in [rad/sec]

% % This speed is the angular speed or how fast we're  revolving around the circle
% % not the same as linear speed in terms of arc length unless radius is the
% % same. But we have no definition of radius.
% 
% % Infact, our simulation even without an explicit radius shows that
% % mod(theta, 2*pi) is a circle centered at the origin with radius one.

number_of_laps = Duration/T_lap; %number of laps made by the rat around

% the circle.
% %dist = 2*pi--- the distance around a circle of radius 1

theta =  speed*time;  % simulated position of the rat at t = time in [s] 



%% call SimHeterogeneousSpikes to generate SpikeTimes, theta and time.
%Output arguments:
% SpikeTimes -- one-by-ncols cell array.
% the i^th cell in SpikeTimes corresponds to spikes times 
% of the i^th neuron
%  theta -- simulated position of the rat to be recovered by 
% the diffusionMaps dim reduction algorithm
% FireRate -- a length(time)-by-ncols matrix containing the firing rates 
% of the neurons.
% each row of FireRate represents the firing behaviour of all the
% ncols neurons at a given rat position, theta

% time -- length(time)-by-one vector indexing the rows of the FireRate data 
% matrix. The units of time are in [s]

%% Generate Synthetic Data using Prevtime and Nextime functions.

disp('Simulating rat position, firing rate and SpikeTimes');tic;
[SpikeTimes, FireRate] = SimHetSpikes_noisy_bimodal(tstart, tend, dt);
toc;

disp('Computing mean of prevtime and nextime to generate NewSpikes');tic;
NewSpikes=BoundaryCondition_NANs_noisy_bimodal(SpikeTimes, dt, time, tstart, tend);
toc;

disp('Computing clean prevtime and nextime functions for synthetic data');
tic;

[prevtime, nextime] = RecomputePrevtimeNextime_noisy_bimodal(NewSpikes, dt, time, tstart, tend);
toc;


%% Generate the distance matrix using the exponential kernel

disp('run Diffusion Maps Algorithm');tic;

tau = 0.1; % tuning parameter for the exponential kernel

prevData = (1/tau)*exp(-abs(prevtime)/tau); %previous time with noisy firing rate



%% call diffusion maps algorithm on [exp-(abs(pretime)) exp-(abs(nexttime))]

%[mappedX1, mapping1] = compute_mapping(X1, 'DiffusionMaps');

% call diffusion maps algorithm on [-exp(prevtime)] only

sig_d = 1; %variance of gaussian in diffusion maps (default = 1)

alpha = 0.5; %parameter in diffusion maps (default = 1)

no_dims = length(NewSpikes); %number of dimensions. (default = 2)

%% Run diffusion Maps on noisy  data


% Run diffusion Maps on  prevtime/nextime data
 
[preVdiffmapL2dist, diffL2struct] = compute_mapping_sylvia(prevData, 'DiffusionMaps', no_dims, alpha, sig_d);

[preVdiffmapL1dist, diffL1struct] = compute_mapping_sylaplacekernel(prevtime, 'DiffusionMaps', no_dims, alpha, sig_d);

% Run diffusion Maps on Firing rate data

[FRdiffmapL2dist, FRdiffL2struct] = compute_mapping_sylvia(FireRate, 'DiffusionMaps', no_dims, alpha, sig_d);

[FRdiffmapL1dist, FRdiffL1struct] = compute_mapping_sylaplacekernel(FireRate, 'DiffusionMaps', no_dims, alpha, sig_d, tau);


% Run Principle Component  Analysis (PCA) on  prevtime/nextime data

[prevPcaL2dist, pcaL2struct] = compute_mapping_sylvia(prevData, 'PCA',  no_dims);

% Run Principle Component Analysis (PCA) on Firing rate data

[FRpcaL2dist, FRpcaL2struct] = compute_mapping_sylvia(FireRate,'PCA',no_dims);



% future ========================================================================================
% % Run Kernel PCA with the Gaussian Kernel on  prevtime/nextime data
%  
% [preVgKpca, gKpcaStruct] = compute_mapping_sylvia(prevData, 'KernelPCA', no_dims, 'gauss', sig_d);
% 
% % Run Kernel PCA  with gaussian kernel on Firing rate data
% 
% [FRgKpca, FRgKpcaStruct] = compute_mapping_sylvia(FireRate, 'KernelPCA', no_dims, 'gauss', sig_d);
% 
% % Run Kernel PCA with the linear Kernel on  prevtime/nextime data
%  
% [preVlKpca, lKpcaStruct] = compute_mapping_sylvia(prevData, 'KernelPCA', no_dims, 'linear', sig_d);
% 
% % Run Kernel PCA  with gaussian kernel on Firing rate data
% 
% [FRlKpca, FRlKpcaStruct] = compute_mapping_sylvia(FireRate, 'KernelPCA', no_dims, 'linear', sig_d);
% 
% % Run Gaussian Process latent variable model on  prevtime/nextime data
%  
% [preVgplvm, gplvmStruct] = compute_mapping_sylvia(prevData, 'GPLVM', no_dims, sig_d);
% 
% % Run Gaussian Process latent variable model on Firing rate data
% 
% [FRgplvm, FRgplvmStruct] = compute_mapping_sylvia(FireRate,'GPLVM', no_dims, sig_d);


 












 








