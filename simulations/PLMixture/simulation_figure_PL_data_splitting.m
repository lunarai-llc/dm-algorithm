%%
tic;
truth = struct('pi',[0.3;0.7], 'theta',[3 1; 0.5 0.5]);  % PL: [mu;sigma]
out = mix_mc_identifyK_perMethodFast_v2( ...
    1, 20000, 'pl', truth, ...          % B=20 gives meaningful percentages
    'Kmax', 5, 'n_sub', 5000, ...        % stronger subsample
    'methods', {'em','hd','vned'}, ...
    'sel_C', 3, ...                      % a few selection splits per dataset
    'sel_nStarts', 5, 'fit_nStarts', 5, ...  % multistarts to avoid bad local optima
    'gh_n', 20, ...                      % more accurate GH quadrature
    'UseParallel', false, 'Verbose', true);
out.khat
out.prop_correct
out.avg_overall
out.avg_when_correct
toc



