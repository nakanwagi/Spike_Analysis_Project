%Recompute PrevTime and NextTime again

disp('recomputing prevtime/nextime'); tic;


prevtime = nan(length(time), ncols);
nextime = nan(length(time), ncols);
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

%% Run diffusion Maps Alagorithm

disp('running diffusion maps algorithm'); tic;


tau = 1;
X1 = [(1/tau)*exp(-abs(prevtime)/tau)  (1/tau)*exp(-abs(nextime)/tau)];
X2 = (1/tau)*exp(-abs(prevtime)/tau);
X3 = (1/tau)*exp(-abs(nextime)/tau);

% at tau = 0.5, the wiggly tail at the end of the experiment becomes
% shorter in [next prev] but the trajectories get more rugged and less
% smooth in prev/next seperately and start to cross each other.
% at tau = 0.1; the different laps start becoming indistinguishable
% tau = 0.05 all the structure is lost and trajectories are so jagged and
% indistiguishable
% tau = 0.25, the wiggly tails in [prev next] become shorter while
% these tails get longer in next/prev time seperately.
% also the circle becomes more smooth and laps are less distinguishable.


%=============Use the Dim reduction tool box =============================



% call diffusion maps algorithm on [exp-(abs(pretime)) exp-(abs(nexttime))]

[mappedX1, mapping1] = compute_mapping(X1, 'DiffusionMaps');

% call diffusion maps algorithm on [-exp(prevtime)] only
 
[mappedX2, mapping2] = compute_mapping(X2, 'DiffusionMaps');

 
% call diffusion maps algorithm on [-exp(nextime)] only
 
[mappedX3, mapping3] = compute_mapping(X3, 'DiffusionMaps');

toc;
