%% no outlier

tic
% Example truth for PG:
truth_pois.pi    = [0.4; 0.6];
truth_pois.theta = [.5; 10];  % a; b

out = mix_mc_identifyK_perMethodFast( ...
    1, 20000, 'poiss', truth_pois, ...
    'methods', {'em','hd','vned'}, ...
    'Kmax', 3,  ...
    'penalty','bic','gh_n',12,'seed',333);

toc
