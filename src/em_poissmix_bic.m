function out = em_poissmix_bic(Y, varargin)
% em_poissmix_bic  EM for Poisson mixtures with BIC/AIC model selection
%   out = em_poissmix_bic(Y, 'Kmax',6, 'nStarts',10, 'seed',0, ...
%                             'maxIter',1000, 'tol',1e-6, 'init','quantile')
%
% Inputs:
%   Y        : nonnegative integer vector (n x 1 or 1 x n)
%   Kmax     : maximum number of components to fit (default 6)
%   nStarts  : multistart EM runs per K (default 10)
%   seed     : RNG seed (default 0)
%   maxIter  : EM max iterations (default 1000)
%   tol      : EM stopping tolerance on loglik (default 1e-6)
%   init     : 'quantile' (default) | 'random'
%
% Outputs (struct):
%   out.Kgrid      : (1..Kmax)'
%   out.loglik     : best log-likelihood per K
%   out.BIC        : -2*loglik + p*log(n), with p = 2K-1
%   out.AIC        :  2*p - 2*loglik
%   out.Khat_BIC   : argmin_K BIC(K)
%   out.fits(K)    : best fit for each K:
%                    .pi (Kx1), .lambda (Kx1), .iters, .converged (bool)
%
% Notes:
%   - Poisson mixture parameters per K: K lambdas + (K-1) weights => p = 2K-1
%   - Uses log-sum-exp for stability. Lambdas are floored at 1e-8, pis at eps.
%   - Labels are not ordered; you can sort components by lambda if desired.

% ---------- parse inputs ----------
p = inputParser;
p.addRequired('Y', @(x) isnumeric(x) && isvector(x) && all(x>=0) && all(mod(x,1)==0));
p.addParameter('Kmax', 6, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('nStarts', 10, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('seed', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('maxIter', 1000, @(x) isnumeric(x) && isscalar(x) && x>=10);
p.addParameter('tol', 1e-6, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('init', 'quantile', @(s) ischar(s) && ismember(lower(s),{'quantile','random'}));
p.parse(Y, varargin{:});
S = p.Results;

Y = Y(:);
n = numel(Y);
rng(S.seed);
Kgrid = (1:S.Kmax).';

loglik = -inf(S.Kmax,1);
BIC    = inf(S.Kmax,1);
AIC    = inf(S.Kmax,1);
fits   = repmat(struct('pi',[],'lambda',[],'iters',0,'converged',false), S.Kmax, 1);

% Precompute const term -gammaln(y+1) (doesn't depend on lambda)
log_y_fact = -gammaln(Y+1);

for K = 1:S.Kmax
    best = struct('ll', -inf, 'pi', [], 'lambda', [], 'iters', 0, 'converged', false);
    for s0 = 1:S.nStarts
        [pi0, lam0] = init_pl(Y, K, S.init);
        [ll, piK, lamK, iters, conv] = em_one(Y, log_y_fact, pi0, lam0, S.maxIter, S.tol);
        if ll > best.ll
            best.ll = ll; best.pi = piK; best.lambda = lamK; best.iters = iters; best.converged = conv;
        end
    end

    pK       = 2*K - 1;                   % free-parameter count
    loglik(K)= best.ll;
    BIC(K)   = -2*best.ll + pK*log(n);
    AIC(K)   =  2*pK - 2*best.ll;

    fits(K).pi       = best.pi;
    fits(K).lambda   = best.lambda;
    fits(K).iters    = best.iters;
    fits(K).converged= best.converged;
end

[~, idx] = min(BIC);
Khat_BIC = Kgrid(idx);

out = struct();
out.Kgrid    = Kgrid;
out.loglik   = loglik;
out.BIC      = BIC;
out.AIC      = AIC;
out.Khat_BIC = Khat_BIC;
out.fits     = fits;
out.settings = S;
end

% ==================== helpers ====================

function [pi0, lam0] = init_pl(Y, K, how)
% initialization: mixing weights & Poisson rates
Y = Y(:);
n = numel(Y);
switch lower(how)
    case 'quantile'
        pi0  = ones(K,1)/K;
        qs   = linspace(5,95,K);
        lam0 = max(1e-3, prctile(Y, qs).'); % column
        if any(~isfinite(lam0))
            mY   = max(1e-3, mean(Y));
            lam0 = mY * linspace(0.5,1.5,K).';
        end
    case 'random'
        a    = rand(K,1) + 0.1;
        pi0  = a/sum(a);
        mY   = max(1e-3, mean(Y));
        sY   = max(1e-3, std(double(Y)));
        lam0 = max(1e-3, mY + sY * randn(K,1));
    otherwise
        error('Unknown init method.');
end
end

function [ll, piK, lamK, iters, converged] = em_one(Y, log_y_fact, piK, lamK, maxIter, tol)
% One EM run (Poisson mixture) with log-sum-exp numerics
n   = numel(Y);
K   = numel(piK);
eps_w = eps;        % weight floor
eps_l = 1e-8;       % lambda floor
converged = false;
prev_ll = -inf;

for it = 1:maxIter
    % ----- E-step: responsibilities via log-sum-exp -----
    % log p(y_i | k) = y_i*log(lam_k) - lam_k - log(y_i!)
    log_lam = log(max(eps_l, lamK(:))).';         % 1 x K
    % build N x K matrix of log pmf terms
    log_pyk = bsxfun(@times, double(Y), log_lam) ...
            - repmat(lamK(:).', n, 1) ...
            + repmat(log_y_fact, 1, K);

    log_pi  = log(max(eps_w, piK(:))).';          % 1 x K
    log_rik = bsxfun(@plus, log_pyk, log_pi);     % N x K
    % log-sum-exp across components
    lse     = logsumexp(log_rik, 2);              % N x 1
    ll      = sum(lse);                           % log-likelihood

    % responsibilities tau = exp(log_rik - lse)
    tau = exp(bsxfun(@minus, log_rik, lse));

    % ----- M-step -----
    Nk    = sum(tau, 1)';                         % K x 1
    piK   = max(eps_w, Nk / n);                   % avoid zeros
    piK   = piK / sum(piK);                       % renormalize exactly
    lamK  = max(eps_l, (tau.' * double(Y)) ./ max(eps_w, Nk));  % K x 1

    % ----- check convergence -----
    if it > 1 && abs(ll - prev_ll) < tol*max(1,abs(prev_ll))
        converged = true;
        break;
    end
    prev_ll = ll;
end

iters = it;
end

function s = logsumexp(A, dim)
% numerically stable log(sum(exp(A), dim))
if nargin < 2, dim = 1; end
amax = max(A, [], dim);
% subtract amax, sum, then add back
s = amax + log( sum( exp(bsxfun(@minus, A, amax)), dim) );
end
