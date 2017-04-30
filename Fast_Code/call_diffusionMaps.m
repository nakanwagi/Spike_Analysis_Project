%% Run diffusion Maps Alagorithm
disp('run diffusion Maps Algorithm');tic;

tau = 1; %our second parameter.

X1 = [(1/tau)*exp(-abs(prevtime)/tau) (1/tau)*exp(-abs(nextime)/tau)];
X2 = (1/tau)*exp(-abs(prevtime)/tau);
X3 = (1/tau)*exp(-abs(nextime)/tau);

%X = exp(-NextTime);

%=============Use the Dim reduction tool box =============================



% call diffusion maps algorithm on [exp-(abs(pretime)) exp-(abs(nexttime))]

[mappedX1, mapping1] = compute_mapping(X1, 'DiffusionMaps');

% call diffusion maps algorithm on [-exp(prevtime)] only
 
[mappedX2, mapping2] = compute_mapping(X2, 'DiffusionMaps');

 
% call diffusion maps algorithm on [-exp(nextime)] only
 
[mappedX3, mapping3] = compute_mapping(X3, 'DiffusionMaps');


toc;
