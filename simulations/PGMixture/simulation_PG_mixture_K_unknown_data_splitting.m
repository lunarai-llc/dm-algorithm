

%% no outlier

tic
% Example truth for PG:
truth_pg.pi    = [0.3; 0.7];
truth_pg.theta = [10 1; 1 2];  % a; b

out = mix_mc_identifyK_perMethodFast( ...
    1, 20000, 'pg', truth_pg, ...
    'methods', {'em','hd','vned'}, ...
    'Kmax', 3,  ...
    'penalty','bic','gh_n',12,'seed',333);

toc


%% summarize data no outlier

T = summarize_pg_results_data_splitting(out);   % creates table T and prints it
% Optionally save to CSV:
% writetable(T, 'pg_summary_table.csv');


%% model misspecification: 

tic
% Generator (Poisson truth)
gen.poiss.pi    = [0.3; 0.7];
gen.poiss.theta = [10; 20];   % Poisson rates

out_model_mis = mix_mc_identifyK_misspecFast( ...
    48, 5000, 'pg', gen.poiss, ...
    'genModel','poiss', ...
    'methods',{'em','hd','vned'}, ...
    'Kmax',5, ...
    'n_sub',1000, 'sel_C',1, 'sel_nStarts',1, 'sel_MaxIter',60, ...
    'fit_C',5,    'fit_nStarts',2, 'fit_MaxIter',100, ...
    'penalty','bic', 'sel_rule','stop', 'init','quantile', ...
    'gh_n',12, 'seed',7, 'UseParallel',true);

toc


