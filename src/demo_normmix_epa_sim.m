function T = demo_normmix_epa_sim(Nrep, NumWorkers)
% DEMO_NORMMIX_EPA_SIM_PARALLEL
% Parallel (parfor) simulation for Normal mixtures with Epanechnikov KDE.
% Methods: EM, HMIX (HD), VNEDMIX (VNED). Selection via stopping rule.
%
% Usage:
%   T = demo_normmix_epa_sim_parallel(5000, 8);

if nargin<1, Nrep = 5000; end
if nargin<2, NumWorkers = 8; end

% ----- Spin up a pool (reuse if already up) -----
p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= NumWorkers
    if ~isempty(p), delete(p); end
    parpool('local', NumWorkers);
end

% Make sure workers can see selectK_normmix.m
w = which('selectK_normmix.m');
if ~isempty(w)
    addAttachedFiles(gcp, {w});
end

% ----- True parameters -----
pi1 = 0.3;
mu1 = 10;  sig1 = 1;
mu2 = 0;   sig2 = 1;
n   = 200;

% ----- Settings -----
Kmax    = 5;
nStarts = 3;        % increase for more robustness (slower)
MaxIter = 400;
Tol     = 1e-6;
penalty = 'bic';
kernel  = 'epa';
grid_n  = 512;
methods = {'em','hd','vned'}; M = numel(methods);

% Precompute independent seeds per repetition (reproducible)
rng(1);                             % master seed
seeds = randi(2^31-1, Nrep, 1);

% ----- Preallocate outputs (Nrep x M), NaN where K≠2 -----
K2        = false(Nrep, M);         % indicator Khat_stop==2
pi10_mat  = nan(Nrep, M);           % larger-mean component (≈10)
mu10_mat  = nan(Nrep, M);
sig10_mat = nan(Nrep, M);
pi0_mat   = nan(Nrep, M);           % smaller-mean component (≈0)
mu0_mat   = nan(Nrep, M);
sig0_mat  = nan(Nrep, M);

% ----- PARFOR over repetitions -----
parfor r = 1:Nrep
    % Per-iteration RNG
    rng(seeds(r), 'CombRecursive');

    % Generate data
    n1 = sum(rand(n,1) < pi1);
    Y  = [mu1 + sig1*randn(n1,1); mu2 + sig2*randn(n-n1,1)];

    for m = 1:M
        meth = methods{m};
        switch meth
            case 'em'
                out = selectK_normmix(Y, 'method','em','em_div','hd', ...
                    'kernel',kernel,'Kmax',Kmax, ...
                    'nStarts',nStarts,'MaxIter',MaxIter,'Tol',Tol, ...
                    'penalty',penalty,'select_rule','stop', ...
                    'grid_n',grid_n,'min_sigma',1e-3,'seed',seeds(r)+101);
            case 'hd'
                out = selectK_normmix(Y, 'method','hd', ...
                    'kernel',kernel,'Kmax',Kmax, ...
                    'nStarts',nStarts,'MaxIter',MaxIter,'Tol',Tol, ...
                    'penalty',penalty,'select_rule','stop', ...
                    'grid_n',grid_n,'min_sigma',1e-3,'seed',seeds(r)+202);
            case 'vned'
                out = selectK_normmix(Y, 'method','vned', ...
                    'kernel',kernel,'Kmax',Kmax, ...
                    'nStarts',nStarts,'MaxIter',MaxIter,'Tol',Tol, ...
                    'penalty',penalty,'select_rule','stop', ...
                    'grid_n',grid_n,'min_sigma',1e-3,'seed',seeds(r)+303);
        end

        if out.Khat_stop == 2
            K2(r,m) = true;
            fit = out.fits(2);
            pi_hat    = fit.pi(:);
            mu_hat    = fit.mu(:);
            sigma_hat = fit.sigma(:);

            % Map larger mean to "≈10" component
            [~, idx_max] = max(mu_hat);
            idx_min = 3 - idx_max;

            pi10_mat(r,m)  = pi_hat(idx_max);
            mu10_mat(r,m)  = mu_hat(idx_max);
            sig10_mat(r,m) = sigma_hat(idx_max);

            pi0_mat(r,m)   = pi_hat(idx_min);
            mu0_mat(r,m)   = mu_hat(idx_min);
            sig0_mat(r,m)  = sigma_hat(idx_min);
        end
    end
end

% ----- Aggregate -----
PctK2  = 100 * mean(K2, 1).';           % Mx1
cntEst = sum(K2, 1).';                  % Mx1

% Helper to get nan-aware mean/std along rows
nanmean_col = @(A) mean(A, 1, 'omitnan').';
nanstd_col  = @(A) std(A, 0, 'omitnan').';

pi10_mean  = nanmean_col(pi10_mat);
pi10_std   = nanstd_col( pi10_mat);
mu10_mean  = nanmean_col(mu10_mat);
mu10_std   = nanstd_col( mu10_mat);
sig10_mean = nanmean_col(sig10_mat);
sig10_std  = nanstd_col( sig10_mat);

pi0_mean   = nanmean_col(pi0_mat);
pi0_std    = nanstd_col( pi0_mat);
mu0_mean   = nanmean_col(mu0_mat);
mu0_std    = nanstd_col(  mu0_mat);
sig0_mean  = nanmean_col(sig0_mat);
sig0_std   = nanstd_col(  sig0_mat);

Method = upper(methods(:));

T = table( Method, ...
    PctK2, cntEst, ...
    pi10_mean, pi10_std, mu10_mean, mu10_std, sig10_mean, sig10_std, ...
    pi0_mean,  pi0_std,  mu0_mean,  mu0_std,  sig0_mean,  sig0_std );

disp('Summary over runs with Khat\_stop = 2 (Epanechnikov KDE):');
disp(T);
end
