function out = selectK_DIC(Y, varargin)
% selectK_DIC  Model selection for Poisson mixtures (Disparity IC + EM/BIC)
% Uses:
%   - HD (Hellinger) via HMIX (your algorithm, generalized to K)
%     H^2 = 2*(1 - sum_x sqrt(g(x) f(x)))
%   - VNED via VNEDMIX (your algorithm, generalized to K)
%     D_VNED = sum_x g(x) * exp( - f(x)/g(x) )
%   - EM/BIC (likelihood) in same call
%
% Example:
% out = selectK_DIC(Y,'method','hd','penalty','bic','Kmax',6,...
%                   'nStarts',10,'seed',1,'dic_solver','auto',...
%                   'em_nStarts',10,'do_dic',true,'do_em',true);

% --------- Parse inputs ---------
p = inputParser;
p.addRequired('Y', @(x) isnumeric(x) && isvector(x) && all(x>=0) && all(mod(x,1)==0));
p.addParameter('method','hd', @(s) ischar(s)&&ismember(lower(s),{'hd','vned'}));
p.addParameter('penalty','aic', @(s) ischar(s)&&ismember(lower(s),{'aic','bic'}));
p.addParameter('Kmax', 6, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('nStarts', 10, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('seed', 0, @(x) isnumeric(x)&&isscalar(x));
% DIC solver: auto ? HD:HMIX, VNED:VNEDMIX
p.addParameter('dic_solver','auto', @(s) ischar(s)&&ismember(lower(s),{'auto','hmix','vnedmix','direct'}));
p.addParameter('dic_maxIter', 500, @(x) isnumeric(x)&&isscalar(x)&&x>=10);
p.addParameter('dic_tol', 1e-7, @(x) isnumeric(x)&&isscalar(x)&&x>0);
% EM settings
p.addParameter('em_nStarts', 10, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('em_maxIter', 1000, @(x) isnumeric(x)&&isscalar(x)&&x>=10);
p.addParameter('em_tol', 1e-6, @(x) isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('em_init','quantile', @(s) ischar(s)&&ismember(lower(s),{'quantile','random'}));
% Toggles
p.addParameter('do_dic', true, @(x) islogical(x)&&isscalar(x));
p.addParameter('do_em',  true, @(x) islogical(x)&&isscalar(x));
p.parse(Y, varargin{:});
S = p.Results;

% ---- Resolve/validate solver pairing ----
if strcmpi(S.dic_solver,'auto')
    S.dic_solver = iff(strcmpi(S.method,'hd'),'hmix','vnedmix');
end
if strcmpi(S.method,'hd')   && ~ismember(lower(S.dic_solver),{'hmix','direct'})
    error('For method="hd", dic_solver must be "hmix" or "direct".');
end
if strcmpi(S.method,'vned') && ~ismember(lower(S.dic_solver),{'vnedmix','direct'})
    error('For method="vned", dic_solver must be "vnedmix" or "direct".');
end

Y = Y(:);
n = numel(Y);
rng(S.seed);

% --------- Empirical pmf g_n on observed support ---------
a  = tabulate(Y);
YY = a(:,1);
g_n = a(:,3)/100;           % strictly >0 on observed support
support = YY;

% --------- b(n) for DIC penalty ---------
switch lower(S.penalty)
    case 'aic', b_n = 1;
    case 'bic', b_n = log(n)/2;
end

% Parameter count v(K) for Poisson mixture: K lambdas + (K-1) weights
vK  = @(K) 2*K - 1;   % for DIC/HIC
pKf = @(K) 2*K - 1;   % for BIC/AIC (likelihood)
Kgrid = (1:S.Kmax).';

% =========================================
% Part 1: DIC (HD?HMIX, VNED?VNEDMIX; or 'direct' if explicit)
% =========================================
D = []; DIC = []; eta = []; fits = []; Khat_DIC = []; Khat_stop = [];
if S.do_dic
    D    = nan(S.Kmax,1);
    DIC  = nan(S.Kmax,1);
    fits = repmat(struct('pi',[],'theta',[],'BC',[],'HD2',[],'iters',0), S.Kmax, 1);

    for K = 1:S.Kmax
        best = struct('obj', Inf, 'pi', [], 'lambda', [], 'BC', [], 'HD2', [], 'iters', 0);
        for s0 = 1:S.nStarts
            [pi0, lam0] = init_pl_for_dic(Y, K);

            switch lower(S.dic_solver)
                case 'hmix'     % HD only
                    [piK, lamK, HD2, BC, iters] = hmix_K(support, g_n, pi0, lam0, S.dic_maxIter, S.dic_tol);
                    obj_here = HD2;

                case 'vnedmix'  % VNED only
                    [piK, lamK, Dvned, iters] = vnedmix_K(support, g_n, pi0, lam0, S.dic_maxIter, S.dic_tol);
                    BC = []; HD2 = [];
                    obj_here = Dvned;

                otherwise       % 'direct' (explicit opt-in)
                    init.pi = pi0; init.lambda = lam0;
                    [obj_here, piK, lamK, BC, HD2] = dic_direct_opt(support, g_n, K, S.method, init);
                    iters = NaN;
            end

            if obj_here < best.obj
                best.obj = obj_here; best.pi = piK; best.lambda = lamK;
                best.BC = BC; best.HD2 = HD2; best.iters = iters;
            end
        end

        fits(K).pi    = best.pi;
        fits(K).theta = struct('lambda', best.lambda);
        fits(K).BC    = best.BC;
        fits(K).HD2   = best.HD2;
        fits(K).iters = best.iters;

        if strcmpi(S.method,'hd')
            if isempty(best.HD2)
                f  = max(pois_mix_pmf(support, best.pi, best.lambda), 1e-15);
                BC = sum( sqrt(g_n .* f) ); D(K) = 2*(1-BC);
            else
                D(K) = best.HD2;
            end
        else % VNED
            if strcmpi(lower(S.dic_solver),'vnedmix')
                D(K) = best.obj;
            else
                f  = max(pois_mix_pmf(support, best.pi, best.lambda), 1e-15);
                D(K) = sum( exp(-f./g_n) .* g_n );
            end
        end

        DIC(K) = D(K) + (b_n / n) * vK(K);
    end

    % Stopping rule
    eta = nan(S.Kmax,1);
    for K = 1:S.Kmax-1
        eta(K) = (b_n / n) * (vK(K+1) - vK(K));
    end
    Khat_stop = S.Kmax;
    for K = 1:S.Kmax-1
        if D(K) <= D(K+1) + eta(K), Khat_stop = K; break; end
    end
    [~, idxMin] = min(DIC); Khat_DIC = Kgrid(idxMin);
end

% =========================================
% Part 2: EM for MLE and traditional BIC/AIC
% =========================================
loglik = []; BIC = []; AIC = []; fits_em = []; Khat_BIC = [];
if S.do_em
    loglik  = -inf(S.Kmax,1);
    BIC     =  inf(S.Kmax,1);
    AIC     =  inf(S.Kmax,1);
    fits_em = repmat(struct('pi',[],'lambda',[],'iters',0,'converged',false), S.Kmax, 1);

    log_y_fact = -gammaln(Y+1);  % constant in log pmf

    for K = 1:S.Kmax
        best = struct('ll', -inf, 'pi', [], 'lambda', [], 'iters', 0, 'converged', false);
        for s0 = 1:S.em_nStarts
            [pi0, lam0] = init_pl(Y, K, S.em_init);
            [ll, piK, lamK, iters, conv] = em_one(Y, log_y_fact, pi0, lam0, S.em_maxIter, S.em_tol);
            if ll > best.ll
                best.ll = ll; best.pi = piK; best.lambda = lamK; best.iters = iters; best.converged = conv;
            end
        end
        pK        = pKf(K);
        loglik(K) = best.ll;
        BIC(K)    = -2*best.ll + pK*log(n);
        AIC(K)    =  2*pK     - 2*best.ll;

        fits_em(K).pi        = best.pi;
        fits_em(K).lambda    = best.lambda;
        fits_em(K).iters     = best.iters;
        fits_em(K).converged = best.converged;
    end

    [~, idxB] = min(BIC); Khat_BIC = Kgrid(idxB);
end

% --------- Pack output ---------
out = struct();
out.Kgrid     = Kgrid;
out.support   = support;
out.g_n       = g_n;
% DIC
out.D         = D;
out.DIC       = DIC;
out.Khat_DIC  = Khat_DIC;
out.Khat_stop = Khat_stop;
out.eta       = eta;
out.fits      = fits;
% EM/BIC
out.loglik    = loglik;
out.BIC       = BIC;
out.AIC       = AIC;
out.Khat_BIC  = Khat_BIC;
out.fits_em   = fits_em;
out.settings  = S;
end

% ===================== DIC helpers =====================

function [pi0, lam0] = init_pl_for_dic(Y, K)
pi0  = ones(K,1)/K;
qs   = linspace(5,95,K);
lam0 = max(1e-3, prctile(Y, qs).');
if any(~isfinite(lam0))
    mY   = max(1e-3, mean(Y));
    lam0 = mY * linspace(0.5,1.5,K).';
end
end

function [piK, lamK, HD2, BC, iters] = hmix_K(x, g, pi0, lam0, maxIter, tol)
% HMIX (generalized from your 2-component code)
piK  = pi0(:);
lamK = max(1e-6, lam0(:));
prev = [piK; lamK]; iters = 0;

for it = 1:maxIter
    iters = it;
    % responsibilities on support
    Pk  = poisspdf(repmat(x,1,numel(lamK)), repmat(lamK',numel(x),1));  % |x| x K
    num = bsxfun(@times, Pk, piK(:)');                                   % |x| x K
    den = max(sum(num,2), 1e-15);
    R   = bsxfun(@rdivide, num, den);                                    % |x| x K

    % per-component minimization of Hk
    Kc = numel(lamK);
    fvals = zeros(Kc,1);
    lam_new = lamK;
    opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton',...
                        'MaxIterations',500,'MaxFunctionEvaluations',5e4);
    for k = 0:Kc-1
        rk = R(:,k+1);
        Hk = @(t) -sum( sqrt( max(1e-12,rk).*g .* poisspdf(x, max(1e-12,exp(t))) ) );
        [tk, fk] = fminunc(Hk, log(lamK(k+1)), opts);
        lam_new(k+1) = max(1e-8, exp(tk));
        fvals(k+1)   = fk;  % negative
    end

    % weight update: pi_k ? f_k^2
    w = fvals.^2; if ~all(isfinite(w)), w = ones(Kc,1); end
    pi_new = max(eps, w / sum(w));

    theta_new = [pi_new; lam_new];
    if max(abs(theta_new - prev)) <= tol
        piK = pi_new; lamK = lam_new; break;
    end
    piK = pi_new; lamK = lam_new; prev = theta_new;
end

f  = max(pois_mix_pmf(x, piK, lamK), 1e-15);
BC = sum( sqrt(g .* f) );
HD2 = 2*(1-BC);
end

function [piK, lamK, Dvned, iters] = vnedmix_K(x, g, pi0, lam0, maxIter, tol)
% VNEDMIX (generalized from your 2-component code)
piK  = pi0(:);
lamK = max(1e-6, lam0(:));
prev = [piK; lamK]; iters = 0;

for it = 1:maxIter
    iters = it;
    % responsibilities
    Pk  = poisspdf(repmat(x,1,numel(lamK)), repmat(lamK',numel(x),1));  % |x| x K
    num = bsxfun(@times, Pk, piK(:)');                                   % |x| x K
    den = max(sum(num,2), 1e-15);
    R   = bsxfun(@rdivide, num, den);                                    % |x| x K
    R_safe = max(R, 1e-12);

    % per-component minimization & vned_k
    Kc = numel(lamK);
    vned_k = zeros(Kc,1);
    lam_new = lamK;
    opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton',...
                        'MaxIterations',500,'MaxFunctionEvaluations',5e4);

    for k = 0:Kc-1
        pk = piK(k+1);
        rk = R(:,k+1); rk_safe = R_safe(:,k+1);
        NEDk = @(t) sum( exp( - pk .* poisspdf(x, max(1e-12,exp(t))) ./ (g .* rk_safe) ) .* g .* rk );
        [tk, ~] = fminunc(NEDk, log(lamK(k+1)), opts);
        lam_new(k+1) = max(1e-8, exp(tk));

        pkphi = poisspdf(x, lam_new(k+1));
        vned_k(k+1) = sum( exp( - pk .* pkphi ./ (g .* rk_safe) ) .* (pk .* pkphi) );
        if ~isfinite(vned_k(k+1)), vned_k(k+1) = 1e-12; end
    end

    pi_new = max(eps, vned_k / sum(vned_k));
    theta_new = [pi_new; lam_new];
    if max(abs(theta_new - prev)) <= tol
        piK = pi_new; lamK = lam_new; break;
    end
    piK = pi_new; lamK = lam_new; prev = theta_new;
end

f  = max(pois_mix_pmf(x, piK, lamK), 1e-15);
Dvned = sum( exp(-f./g) .* g );
end

function [obj, piK, lamK, BCval, HD2val] = dic_direct_opt(x, g, K, method, init)
opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton',...
                    'MaxIterations',500,'MaxFunctionEvaluations',5e4);
z0 = pack_params(init, K);
objfun = @(z) dic_objective(z, x, g, K, method);
try
    zhat = fminunc(objfun, z0, opts);
catch
    zhat = z0;
end
[obj, piK, lamK, BCval, HD2val] = dic_objective(zhat, x, g, K, method);
end

function z = pack_params(init, K)
alpha0 = zeros(K,1);
beta0  = log(max(1e-6, init.lambda(:)));
z = [alpha0; beta0] + 0.05*randn(numel(alpha0)+numel(beta0),1);
end

function [obj, piK, lamK, BCval, HD2val] = dic_objective(z, x, g, K, method)
alpha = z(1:K); beta = z(K+1:2*K);
e = exp(alpha - max(alpha)); piK = e / sum(e);
lamK = max(1e-6, exp(beta));
f = max(pois_mix_pmf(x, piK, lamK), 1e-15);
switch lower(method)
    case 'hd'
        BCval  = sum( sqrt(g .* f) );
        HD2val = 2 * (1 - BCval);
        obj    = HD2val;
    case 'vned'
        BCval = []; HD2val = [];
        obj = sum( exp(-f./g) .* g );
    otherwise, error('Unknown method.');
end
end

function f = pois_mix_pmf(x, piK, lambda)
K = numel(piK);
X = repmat(x(:)', K, 1);
Lam = repmat(lambda(:), 1, numel(x));
Pk = poisspdf(X, Lam);
f = (piK(:)' * Pk).';
end

% ===================== EM/BIC helpers =====================

function [pi0, lam0] = init_pl(Y, K, how)
Y = Y(:);
switch lower(how)
    case 'quantile'
        pi0 = ones(K,1)/K;
        qs = linspace(5,95,K);
        lam0 = max(1e-3, prctile(Y, qs).');
        if any(~isfinite(lam0))
            mY = max(1e-3, mean(Y));
            lam0 = mY * linspace(0.5,1.5,K).';
        end
    case 'random'
        a = rand(K,1) + 0.1; pi0 = a/sum(a);
        mY = max(1e-3, mean(Y)); sY = max(1e-3, std(double(Y)));
        lam0 = max(1e-3, mY + sY*randn(K,1));
    otherwise, error('Unknown init method.');
end
end

function [ll, piK, lamK, iters, converged] = em_one(Y, log_y_fact, piK, lamK, maxIter, tol)
n = numel(Y); K = numel(piK); eps_w = eps; eps_l = 1e-8;
converged = false; prev_ll = -inf;
for it = 1:maxIter
    log_lam = log(max(eps_l, lamK(:))).';
    log_pyk = bsxfun(@times, double(Y), log_lam) ...
            - repmat(lamK(:).', n, 1) ...
            + repmat(log_y_fact, 1, K);
    log_pi  = log(max(eps_w, piK(:))).';
    log_rik = bsxfun(@plus, log_pyk, log_pi);
    lse     = logsumexp(log_rik, 2);
    ll      = sum(lse);
    tau     = exp(bsxfun(@minus, log_rik, lse));

    Nk    = sum(tau, 1)';                   
    piK   = max(eps_w, Nk / n);             
    piK   = piK / sum(piK);                 
    lamK  = max(eps_l, (tau.' * double(Y)) ./ max(eps_w, Nk));

    if it > 1 && abs(ll - prev_ll) < tol*max(1,abs(prev_ll)), converged = true; break; end
    prev_ll = ll;
end
iters = it;
end

function s = logsumexp(A, dim)
if nargin < 2, dim = 1; end
amax = max(A, [], dim);
s = amax + log( sum( exp(bsxfun(@minus, A, amax)), dim) );
end

function y = iff(cond, a, b); if cond, y=a; else, y=b; end; end
