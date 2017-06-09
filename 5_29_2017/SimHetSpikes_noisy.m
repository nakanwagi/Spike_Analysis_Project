function SpikeTimes = SimHetSpikes_noisy(tstart,...,
    tend,dt, theta, time, ~, ~, ncols,sigma)
%SimHeterogeneousSpikes generates SpikeTimes given the input arguments below:

%Usage: SpikeTimes = SimHeterogeneousSpikes(tstart, tend, dt)

% Input arguments:
% tstart -- the start time of the experiment in [s]
% tend -- the end time of the experiment in [s]
% dt -- step size/bin width at which the data is sampled in [s] 
% small dt ensures  that the poisson process is approximately
% bernouli since dt<<<1.
% ncols -- number of neurons/cells in the data set from which data is
% recorded
% sigma -- standard deviation of the Gaussian used encode the sensitivity
% of the receptive fields while simulating the firing rate.
% The variable sigma should be a vector of size one-by-ncols 
% We add heterogeneity in sigma so that some receptive fields are
% wider than others.
% (default: sigma = 0.01 + rand(1, ncols)*(1-0.01))

%Output arguments:
% SpikeTimes -- one-by-ncols cell array.
% the i^th cell in SpikeTimes corresponds to spikes times 
% of the i^th neuron
%  theta -- simulated position of the rat to be recovered by 
% the diffusionMaps dim reduction algorithm
% FireRate -- an length(time)-by-ncols matrix containing the firing rates 
% of the neurons.
% each row of FireRate represents the firing behaviour of all the
% ncols neurons at a given rat position theta

% time -- length(time)-by-one vector indexing the rows of the FireRate data 
% matrix. The units of time are in [s]
  
%% Define some default parameters
if nargin == 3
   ncols = 32; % number of neurons or columns of the data
   T_lap  = 1.5; % period or time to make one lap in [s] around circle.
   sig_max = 1; 
   sig_min = 0.01;
   sigma = sig_min + rand(1, ncols)*(sig_max - sig_min);
% we don't want sig_min to be smaller than the smallest
%sigma for which Laplacian eigenMaps breaks down.

%---------intialize the time vector-----------------------------
time = tstart:dt:tend-dt;  % tvec is a vector describing the time 
% stamps for each row of the data
% length(time) = total number of time points;----useful for making plots.

speed = 2*pi/T_lap; %speed of the rat in [rad/sec]
% This speed is the angular speed or how fast we're  revolving around the circle
% not the same as linear speed in terms of arc length unless radius is the
% same. But we have no definition of radius.

% Infact, our simulation even without an explicit radius shows that
% mod(theta, 2*pi) is a circle centered at the origin with radius one.

% number_of_laps = Duration/T_lap; %number of laps made by the rat around
% the circle.
%dist = 2*pi--- the distance around a circle of radius 1

theta =  speed*time;  % simulated position of the rat at t = time in [s] 

end


%% Using uniform sampling
% for each rat position, all the 32 neurons fire according to the
% function FireRate.
%==========================================================================
%case 1: center = (2*pi*i)/N %centers evenly spread out in space.

%=========================================================================
%case 2: centers uniformly distributed in space.
center = rand(1, ncols)*2*pi; % centers must be the same size as N
eta = 0.05; % parameter for std of noise added to each center
center = center + (eta*2*pi)*randn(size(center)); 
%add Gaussian noise with std = eta*2*pi to the centers

%========================================================================
% initialize the firing rate matrix

FireRate = zeros(length(time), ncols);

% Parameters for the maximum and minimum firing rate
M = 20; % maximum firing rate added as a parameter.
S = 400; % parameter for the minimum firing rate.
maxRate = rand(1, ncols)*M; %the maximum firing rate helps to ensure that we have many spiketimes for each neuron
% the larger the max rate, the more spike times we get in each cell
% be aware that some cells may not spike in which case we get an empty
% spike  time matrix.

minRate = maxRate./S;



%%  define an anonymous function to define the firing rate function f 

f=@(x,sig,maxR,minR) ((maxR-minR).*exp(-((mod(x+pi, 2*pi)-pi)).^2./(2.*sig.^2))) + minR;



%%  Add  noise component-wise in the firing rate
eta2 = 0.05; %parameter for std of noise for the firing rate
noise = (eta2*M)*randn(size(FireRate)); % std of gaussian noise = (eta2*M)



%% Compute the noisy firing rate

for i = 1:ncols
    
    % heterogeneous sigma
    FireRate(:,i) = f(theta - center(i), sigma(i), maxRate(i), minRate(i)) + noise(i); 
    %Uniform sampling of neuron positions/centers
    
   
end
%==========================================================================

%check if uniform rand number in (0, 1) is less than prob of 1 spike
% in interval Dt and store any spikes in the matrix isSpike.
isSpike = rand(length(time), ncols) < FireRate*dt;

%extract the time when isSpikes had the value 1 and
%store it in the vector SpikeTime(i) for the i^th cell/neuron.

SpikeTimes = cell(1, ncols);
for ii = 1:ncols
   SpikeTimes{ii} = time(isSpike(:,ii)); 
   %SpikeTimes contains the spike times for the i^th neuron
   
end




end


