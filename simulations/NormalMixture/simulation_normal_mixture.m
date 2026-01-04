

%%
% Start a pool once (optional; function will start it if needed)
parpool('local', 6);


%%
tic

% Run 5000 parallel replications
T = demo_normmix_epa_sim(5000, 6);

toc

% takes 22742.777573 seconds with 6 cores of cpu. 



