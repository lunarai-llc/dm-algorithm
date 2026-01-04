
function out = selectK_normmix(Y, varargin)
% selectK_normmix
% Normal-mixture selection & estimation with EM / HMIX (HD) / VNEDMIX.
% g_n is a KDE (Epanechnikov or Gaussian). Selection uses:
%   HD:   H^2 = 2*(1 - ∫ sqrt(g f) dx)
%   VNED: ∫ g * exp(-f/g) dx
% and the stopping rule:  Khat_stop = min { K : D(K) <= D(K+1)+eta(K) }.
%
% EM fits by likelihood; K selection for EM uses divergence to g_n (choose via 'em_div').
%
% USAGE:
%   out = selectK_normmix(Y, 'method','hd', 'kernel','epa', 'Kmax',5, ...
%       'penalty','bic','select_rule','stop', 'nStarts',10, ...
%       'MaxIter',400,'Tol',1e-6,'grid_n',512,'min_sigma',1e-3,'seed',1);
%
%   % EM but choose K by HD against KDE:
%   out = selectK_normmix(Y,'method','em','em_div','hd','kernel','gauss');
%
% OUTPUT (struct):
%   out.Kgrid, out.D, out.IC, out.eta, out.Khat_stop, out.Khat_min
%   out.fits(K): .pi, .mu, .sigma, .iters
%   out.grid: .x (grid), .g (KDE), .dx
%   out.settings: echo of options

% -------- Parse inputs --------
ip = inputParser;
ip.addRequired('Y', @(x) isnumeric(x) && isvector(x));
ip.addParameter('method','hd', @(s) ischar(s) && ismember(lower(s),{'em','hd','vned'}));
ip.addParameter('em_div','hd', @(s) ischar(s) && ismember(lower(s),{'hd','vned'})); % divergence to pick K for EM
ip.addParameter('kernel','epa', @(s) ischar(s) && ismember(lower(s),{'epa','gauss'}));
ip.addParameter('Kmax',5, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('penalty','bic', @(s) ischar(s)&&ismember(lower(s),{'bic','aic'}));
ip.addParameter('select_rule','stop', @(s) ischar(s)&&ismember(lower(s),{'stop','min'}));
ip.addParameter('nStarts',5, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('MaxIter',400, @(x) isnumeric(x)&&isscalar(x)&&x>=50);
ip.addParameter('Tol',1e-6, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('grid_n',512, @(x) isnumeric(x)&&isscalar(x)&&x>=128);
ip.addParameter('min_sigma',1e-3, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('seed',0, @(x) isnumeric(x)&&isscalar(x));
ip.parse(Y, varargin{:});
S = ip.Results;

rng(S.seed);
Y = Y(:);
n = numel(Y);

% penalty b(n)
switch lower(S.penalty)
    case 'aic', b_n = 1;
    otherwise,  b_n = log(n)/2; % 'bic'
end

% parameter count v(K) = (K-1) + 2K = 3K-1
vK = @(K) (3*K - 1);

% -------- KDE g_n on a common grid --------
[xx, gkde, dx] = kde1d(Y, S.kernel, S.grid_n);

% -------- Loop over K --------
Kgrid = (1:S.Kmax).';
Dvals = nan(S.Kmax,1);
IC    = nan(S.Kmax,1);
fits  = repmat(struct('pi',[],'mu',[],'sigma',[],'iters',0), S.Kmax,1);

for K = 1:S.Kmax
    best = struct('val',Inf,'pi',[],'mu',[],'sigma',[],'iters',0);

    for s0 = 1:S.nStarts
        [pi0, mu0, sig0] = init_norm_params(Y, K, S.min_sigma);

        switch lower(S.method)
            case 'hd'
                fit = hmix_fit(xx, gkde, dx, pi0, mu0, sig0, S.MaxIter, S.Tol, S.min_sigma);
                D   = hd_divergence(xx, gkde, dx, fit.pi, fit.mu, fit.sigma);
                val = D;

            case 'vned'
                fit = vnedmix_fit(xx, gkde, dx, pi0, mu0, sig0, S.MaxIter, S.Tol, S.min_sigma);
                D   = vned_divergence(xx, gkde, dx, fit.pi, fit.mu, fit.sigma);
                val = D;

            otherwise % 'em' (fit on Y), then compute divergence to KDE for selection
                fit = em_gauss_fit(Y, K, pi0, mu0, sig0, S.MaxIter, S.Tol, S.min_sigma);
                switch lower(S.em_div)
                    case 'hd'
                        D = hd_divergence(xx, gkde, dx, fit.pi, fit.mu, fit.sigma);
                    otherwise
                        D = vned_divergence(xx, gkde, dx, fit.pi, fit.mu, fit.sigma);
                end
                val = D;
        end

        if val < best.val
            best.val   = val;
            best.pi    = fit.pi;
            best.mu    = fit.mu;
            best.sigma = fit.sigma;
            best.iters = fit.iters;
        end
    end

    % record best for this K
    fits(K).pi    = best.pi;
    fits(K).mu    = best.mu;
    fits(K).sigma = best.sigma;
    fits(K).iters = best.iters;

    Dvals(K) = best.val;
    IC(K)    = Dvals(K) + (b_n/n) * vK(K);
end

% -------- Stopping rule eta(K) & Khat --------
eta = nan(S.Kmax,1);
for K = 1:S.Kmax-1
    eta(K) = (b_n/n) * (vK(K+1) - vK(K));
end

Khat_stop = S.Kmax;
for K = 1:S.Kmax-1
    if Dvals(K) <= Dvals(K+1) + eta(K)
        Khat_stop = K; break;
    end
end
[~,iMin] = min(IC);
Khat_min = Kgrid(iMin);

% -------- Pack output --------
out = struct();
out.Kgrid     = Kgrid;
out.D         = Dvals;
out.IC        = IC;
out.eta       = eta;
out.Khat_stop = Khat_stop;
out.Khat_min  = Khat_min;
out.fits      = fits;
out.grid      = struct('x',xx,'g',gkde,'dx',dx);
out.settings  = S;
end

% =======================================================================
% =========================== KDE (g_n) =================================
% =======================================================================
function [xgrid, g, dx] = kde1d(Y, kernel, grid_n)
Y = Y(:);
n = numel(Y);
s = std(Y);
if s<=0, s = 1; end
% Silverman-like bandwidth (works for both kernels)
h = 1.06 * s * n^(-1/5);
if h < 1e-3, h = 1e-3; end

% Grid
lo = min(Y) - 3*s;
hi = max(Y) + 3*s;
xgrid = linspace(lo, hi, grid_n).';
dx = xgrid(2)-xgrid(1);

Z = (xgrid - Y.')/h;  % grid_n x n

switch lower(kernel)
    case 'epa'
        K = max(0, 0.75*(1 - Z.^2));           % 0 outside |Z|>1
        K(abs(Z) > 1) = 0;
        g = sum(K, 2) / (n*h);
    otherwise % 'gauss'
        K = exp(-0.5*Z.^2) / sqrt(2*pi);
        g = sum(K, 2) / (n*h);
end
g = max(g, 1e-12);
end

% =======================================================================
% ===================== Mixture PDF on grid =============================
% =======================================================================
function f = mixpdf_norm(xx, piK, mu, sigma)
% all column vectors
piK   = piK(:);
mu    = mu(:);
sigma = sigma(:);
K = numel(piK);

% Evaluate component PDFs: each column is a component
X = xx(:);
comp = zeros(numel(X), K);
for k = 1:K
    comp(:,k) = normal_pdf(X, mu(k), max(1e-9, sigma(k)));
end
f = max(comp * piK, 1e-12);
end

function y = normal_pdf(x, mu, sigma)
z = (x - mu) ./ sigma;
y = exp(-0.5*z.^2) ./ (sqrt(2*pi) * sigma);
end

% =======================================================================
% ====================== Divergence Functionals =========================
% =======================================================================
function D = hd_divergence(xx, g, dx, piK, mu, sigma)
% H^2 = 2*(1 - ∫ sqrt(g f) dx)   (Riemann sum)
f  = mixpdf_norm(xx, piK, mu, sigma);
BC = sum( sqrt(g .* f) ) * dx;
D  = 2 * max(0, 1 - BC);
end

function D = vned_divergence(xx, g, dx, piK, mu, sigma)
% ∫ g * exp(-f/g) dx
f = mixpdf_norm(xx, piK, mu, sigma);
D = sum( g .* exp(-f./g) ) * dx;
end

% =======================================================================
% ========================= Initialisation ==============================
% =======================================================================
function [pi0, mu0, sig0] = init_norm_params(Y, K, min_sigma)
Y = Y(:);
pi0 = ones(K,1)/K;

% k-means++ like seeds via quantiles (robust)
qs = linspace(5,95,K);
mu0 = prctile(Y, qs).';
if any(~isfinite(mu0))
    mu0 = linspace(min(Y), max(Y), K).';
end

% common sigma start
sY = std(Y);
if sY<=0, sY = 1; end
sig0 = max(min_sigma, 0.5*sY) * ones(K,1);
end

% =======================================================================
% ============================== HMIX ===================================
% =======================================================================
function fit = hmix_fit(xx, g, dx, pi0, mu0, sig0, MaxIter, Tol, min_sigma)
% Alternating scheme on grid: responsibilities wrt current mix f, then
% per-component (mu,sigma) via fminunc to maximize component-wise BC proxy,
% weights updated via squared objective values (as in your HMIX).
piK   = pi0(:); mu = mu0(:); sigma = max(min_sigma, sig0(:));
prev  = [piK; mu; sigma];
iters = 0;

opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                    'MaxIterations',400,'MaxFunctionEvaluations',5e4);

while iters < MaxIter
    iters = iters + 1;

    % E-step on grid
    K = numel(piK);
    comp = zeros(numel(xx), K);
    for k = 1:K, comp(:,k) = normal_pdf(xx, mu(k), sigma(k)); end
    num = bsxfun(@times, comp, piK(:).');
    den = max(sum(num,2), 1e-15);
    R   = bsxfun(@rdivide, num, den);              % responsibilities on grid

    % Per-component update (mu,sigma)
    fvals = zeros(K,1);
    mu_new = mu; sig_new = sigma;
    for k = 1:K
        rk = R(:,k);
        % objective (negative BC_k proxy); dx cancels in pi update ratio
        obj = @(v) -sum( sqrt( max(1e-12, rk).*g .* normal_pdf(xx, v(1), max(min_sigma, exp(v(2)))) ) );
        v0 = [mu(k); log(sigma(k))];
        [vopt, fneg] = fminunc(obj, v0, opts);
        mu_new(k)  = vopt(1);
        sig_new(k) = max(min_sigma, exp(vopt(2)));
        fvals(k)   = fneg;   % negative value
    end

    % weight update: pi_k ∝ f_k^2
    w = fvals.^2;
    if ~all(isfinite(w)), w = ones(K,1); end
    pi_new = max(eps, w/sum(w));

    theta_new = [pi_new; mu_new; sig_new];
    if max(abs(theta_new - prev)) <= Tol
        piK = pi_new; mu = mu_new; sigma = sig_new; break;
    end
    piK = pi_new; mu = mu_new; sigma = sig_new; prev = theta_new;
end

fit = struct('pi',piK,'mu',mu,'sigma',sigma,'iters',iters);
end

% =======================================================================
% ============================= VNEDMIX =================================
% =======================================================================
function fit = vnedmix_fit(xx, g, dx, pi0, mu0, sig0, MaxIter, Tol, min_sigma)
% Alternating scheme on grid for VNED
piK   = pi0(:); mu = mu0(:); sigma = max(min_sigma, sig0(:));
prev  = [piK; mu; sigma];
iters = 0;

opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                    'MaxIterations',400,'MaxFunctionEvaluations',5e4);

while iters < MaxIter
    iters = iters + 1;

    % Responsibilities on grid
    K = numel(piK);
    comp = zeros(numel(xx), K);
    for k = 1:K, comp(:,k) = normal_pdf(xx, mu(k), sigma(k)); end
    num = bsxfun(@times, comp, piK(:).');
    den = max(sum(num,2), 1e-15);
    R   = bsxfun(@rdivide, num, den);
    R_safe = max(R, 1e-12);

    % Per-component update to minimize NED_k, then pi update via vned_k
    vned_k = zeros(K,1);
    mu_new = mu; sig_new = sigma;

    for k = 1:K
        pk = piK(k);
        rk = R(:,k); rk_safe = R_safe(:,k);

        % minimize sum exp(- pk * N_k / (g*rk)) * g * rk (dx cancels in ratio)
        obj = @(v) sum( exp( - pk .* normal_pdf(xx, v(1), max(min_sigma,exp(v(2)))) ./ (g .* rk_safe) ) .* g .* rk );
        v0  = [mu(k); log(sigma(k))];
        vopt = fminunc(obj, v0, opts);
        mu_new(k)  = vopt(1);
        sig_new(k) = max(min_sigma, exp(vopt(2)));

        % compute vned_k for pi-update
        Nk = normal_pdf(xx, mu_new(k), sig_new(k));
        vned_k(k) = sum( exp( - pk .* Nk ./ (g .* rk_safe) ) .* (pk .* Nk) );
        if ~isfinite(vned_k(k)), vned_k(k) = 1e-12; end
    end

    pi_new = max(eps, vned_k / sum(vned_k));

    theta_new = [pi_new; mu_new; sig_new];
    if max(abs(theta_new - prev)) <= Tol
        piK = pi_new; mu = mu_new; sigma = sig_new; break;
    end
    piK = pi_new; mu = mu_new; sigma = sig_new; prev = theta_new;
end

fit = struct('pi',piK,'mu',mu,'sigma',sigma,'iters',iters);
end

% =======================================================================
% =============================== EM ====================================
% =======================================================================
function fit = em_gauss_fit(Y, K, pi0, mu0, sig0, MaxIter, Tol, min_sigma)
% Standard Gaussian-mixture EM on sample Y
Y = Y(:); n = numel(Y);
piK = pi0(:); mu = mu0(:); sigma = max(min_sigma, sig0(:));
prev_ll = -Inf;
iters = 0;

while iters < MaxIter
    iters = iters + 1;

    % E-step
    comp = zeros(n, K);
    for k = 1:K
        comp(:,k) = normal_pdf(Y, mu(k), sigma(k));
    end
    mix = max(comp * piK, 1e-300);
    tau = bsxfun(@rdivide, bsxfun(@times, comp, piK(:).'), mix);

    % M-step
    Nk = sum(tau, 1)';       % K x 1
    piK = max(eps, Nk / n);
    mu  = (tau' * Y) ./ max(Nk, eps);
    % var update
    sig2 = zeros(K,1);
    for k = 1:K
        sig2(k) = sum( tau(:,k) .* (Y - mu(k)).^2 ) / max(Nk(k), eps);
    end
    sigma = sqrt(max(min_sigma^2, sig2));

    % loglik for convergence
    ll = sum(log(mix));
    if isfinite(prev_ll) && abs(ll - prev_ll) < Tol*max(1,abs(prev_ll))
        break;
    end
    prev_ll = ll;
end

fit = struct('pi',piK,'mu',mu,'sigma',sigma,'iters',iters);
end




