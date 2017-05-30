%% Generate the distance matrix using the exponential kernel

disp('run Diffusion Maps Algorithm');tic;

tau = 15; %-----with noise; %0.9----no noise; 

X2 = (1/tau)*exp(-abs(prevtime)/tau);  



%% call diffusion maps algorithm on [exp-(abs(pretime)) exp-(abs(nexttime))]

%[mappedX1, mapping1] = compute_mapping(X1, 'DiffusionMaps');

% call diffusion maps algorithm on [-exp(prevtime)] only

sig_d = 1; %variance of gaussian in diffusion maps (default = 1)

alpha = 1; %parameter in diffusion maps (default = 1)

no_dims = 2; %number of dimensions. (default = 2)

%% Run diffusion Maps on clean data
 
[mappedX2, mapping2] = compute_mapping(X2, 'DiffusionMaps', no_dims, alpha, sig_d);

 
toc;



