function [avg, est, meta] = demo_pl_outlier_sim()
% DEMO_PL_OUTLIER_SIM
% Monte Carlo for Poisson–Lognormal (PL) mixture with known K=2 and
% ε-fraction outliers (replacing values by a large count).
% Methods: EM (MLE), HMIX (Hellinger divergence), VNEDMIX (VNED).
%
% Returns:
%   avg.(EM|HD|VNED)  -> [E x P] NaN-robust averages across reps
%   est.(EM|HD|VNED)  -> [R x E x P] all estimates (each rep & epsilon)
%   meta              -> struct with cont_grid, truth, etc., and creates a plot.
%
% Parameter order P = [pi1, mu1, sig1, mu2, sig2],
% where PL component is:  lambda ~ LogNormal(mu, sig),  Y|lambda ~ Poisson(lambda).
%
% Truth used here (component 1 has the larger mean for stable alignment):
%   pi1 = 0.3, mu1 = 3, sigma1 = 0.5
%   pi2 = 0.7, mu2 = 1, sigma2 = 0.5
% Mean(lambda_k) = exp(mu_k + 0.5*sigma_k^2).

% ------------------ Simulation settings ------------------
n            = 2000;                            % sample size per dataset
R            = 5000;                            % Monte Carlo reps
cont_grid    = [0 0.05 0.10 0.15 0.20 0.30];    % ε levels (outlier pct)
outlier_val  = 50;                              % value used as outlier
nStarts0     = 6;                               % multistarts at ε=0
nStartsWarm  = 1;                               % warm-starts for next ε
MaxIter      = 150;                             % inner EM/optimizer iters
Tol          = 1e-6;                            % convergence tolerance
init         = 'quantile';                      % 'quantile' | 'random'
nGH          = 20;                              % Gauss–Hermite nodes (accuracy/speed)
rng(1);                                         % base RNG

% Truth (component 1 has the HIGHER PL mean for stable alignment)
pi_true = [0.3; 0.7];
mu_true = [3; 1];
sg_true = [0.5; 0.5];
K       = 2;

% ------------------ Storage ------------------
E = numel(cont_grid);
P = 5;  % [pi1, mu1, sig1, mu2, sig2]
est_EM   = NaN(R,E,P);
est_HD   = NaN(R,E,P);
est_VNED = NaN(R,E,P);

% ------------------ Parallel pool ------------------
usePar = license('test','Distrib_Computing_Toolbox');
if usePar
    try
        p = gcp('nocreate'); if isempty(p), parpool('local'); end %#ok<PPOOL>
    catch
        usePar = false; % fall back to serial
    end
end

% ------------------ Monte Carlo (parallel over reps) ------------------
if usePar
parfor r = 1:R
    rng(1000 + r);             % independent stream
    GH = gh_rule(nGH);         % <-- build GH INSIDE each worker (prevents broadcast bug)

    local_EM   = NaN(E,P);
    local_HD   = NaN(E,P);
    local_VNED = NaN(E,P);

    % warm starts carried across epsilon levels for this replication
    prev_em   = [];
    prev_hd   = [];
    prev_vned = [];

    for e = 1:E
        cont = cont_grid(e);

        % --- Generate clean PL mixture, then inject outliers ---
        Y = sim_pl_mix(n, pi_true, mu_true, sg_true);
        m = round(cont*n);
        if m>0
            idx = randperm(n, m);
            Y(idx) = outlier_val;
        end

        % --- EM (known K), with warm start after first ε ---
        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_em = em_pl_knownK(Y, K, nStarts, MaxIter, Tol, init, prev_em, GH);
        if fit_em.converged
            fit_em = align_pl_by_mean(fit_em);
            local_EM(e,:) = [fit_em.pi(1), fit_em.theta(1,1), fit_em.theta(2,1), ...
                                         fit_em.theta(1,2), fit_em.theta(2,2)];
            prev_em = fit_em; % warm start forward
        else
            prev_em = [];     % drop warm start if it failed
        end

        % --- HMIX (HD) ---
        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_hd = dic_pl_knownK(Y, K, 'hd', nStarts, MaxIter, Tol, init, prev_hd, GH);
        fit_hd = align_pl_by_mean(fit_hd);
        local_HD(e,:) = [fit_hd.pi(1), fit_hd.theta(1,1), fit_hd.theta(2,1), ...
                                      fit_hd.theta(1,2), fit_hd.theta(2,2)];
        prev_hd = fit_hd;

        % --- VNEDMIX (VNED) ---
        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_vn = dic_pl_knownK(Y, K, 'vned', nStarts, MaxIter, Tol, init, prev_vned, GH);
        fit_vn = align_pl_by_mean(fit_vn);
        local_VNED(e,:) = [fit_vn.pi(1), fit_vn.theta(1,1), fit_vn.theta(2,1), ...
                                        fit_vn.theta(1,2), fit_vn.theta(2,2)];
        prev_vned = fit_vn;
    end

    est_EM(r,:,:)   = local_EM;
    est_HD(r,:,:)   = local_HD;
    est_VNED(r,:,:) = local_VNED;
end
else
for r = 1:R
    rng(1000 + r);
    GH = gh_rule(nGH);         % build GH in serial path as well

    local_EM   = NaN(E,P);
    local_HD   = NaN(E,P);
    local_VNED = NaN(E,P);

    prev_em   = [];
    prev_hd   = [];
    prev_vned = [];

    for e = 1:E
        cont = cont_grid(e);

        Y = sim_pl_mix(n, pi_true, mu_true, sg_true);
        m = round(cont*n);
        if m>0
            idx = randperm(n, m);
            Y(idx) = outlier_val;
        end

        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_em = em_pl_knownK(Y, K, nStarts, MaxIter, Tol, init, prev_em, GH);
        if fit_em.converged
            fit_em = align_pl_by_mean(fit_em);
            local_EM(e,:) = [fit_em.pi(1), fit_em.theta(1,1), fit_em.theta(2,1), ...
                                         fit_em.theta(1,2), fit_em.theta(2,2)];
            prev_em = fit_em;
        else
            prev_em = [];
        end

        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_hd = dic_pl_knownK(Y, K, 'hd', nStarts, MaxIter, Tol, init, prev_hd, GH);
        fit_hd = align_pl_by_mean(fit_hd);
        local_HD(e,:) = [fit_hd.pi(1), fit_hd.theta(1,1), fit_hd.theta(2,1), ...
                                      fit_hd.theta(1,2), fit_hd.theta(2,2)];
        prev_hd = fit_hd;

        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_vn = dic_pl_knownK(Y, K, 'vned', nStarts, MaxIter, Tol, init, prev_vned, GH);
        fit_vn = align_pl_by_mean(fit_vn);
        local_VNED(e,:) = [fit_vn.pi(1), fit_vn.theta(1,1), fit_vn.theta(2,1), ...
                                        fit_vn.theta(1,2), fit_vn.theta(2,2)];
        prev_vned = fit_vn;
    end

    est_EM(r,:,:)   = local_EM;
    est_HD(r,:,:)   = local_HD;
    est_VNED(r,:,:) = local_VNED;
end
end

% ------------------ Averages (ignore NaNs) ------------------
avg_EM   = squeeze(nanmean(est_EM,1));   % [E x P]
avg_HD   = squeeze(nanmean(est_HD,1));
avg_VNED = squeeze(nanmean(est_VNED,1));

% ------------------ Package outputs ------------------
est  = struct('EM',est_EM, 'HD',est_HD, 'VNED',est_VNED);
avg  = struct('EM',avg_EM, 'HD',avg_HD, 'VNED',avg_VNED);
meta = struct();
meta.cont_grid   = cont_grid;
meta.param_names = {'\pi_1','\mu_1','\sigma_1','\mu_2','\sigma_2'};
meta.truth       = [pi_true(1), mu_true(1), sg_true(1), mu_true(2), sg_true(2)];
meta.n = n; meta.R = R; meta.outlier_val = outlier_val;

% ------------------ Plot (robust y-limits + dotted truth) ------------------
figure('Color','w','Position',[60 60 1100 740]);
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

for p = 1:P
    nexttile;
    plot(cont_grid, avg.EM(:,p),   '-o','LineWidth',1.6,'DisplayName','EM');   hold on;
    plot(cont_grid, avg.HD(:,p),   '-s','LineWidth',1.6,'DisplayName','HMIX');
    plot(cont_grid, avg.VNED(:,p), '-^','LineWidth',1.6,'DisplayName','VNEDMIX');
    yline(meta.truth(p),'k:','LineWidth',1.5,'DisplayName','Truth');

    xlabel('\epsilon (outlier fraction)'); ylabel('Average estimate');
    title(meta.param_names{p},'Interpreter','tex'); grid on; box on;

    % robust axis limits to avoid blown-up non-convergent runs
    vals = [avg_EM(:,p); avg_HD(:,p); avg_VNED(:,p); meta.truth(p)];
    vals = vals(isfinite(vals));
    if isempty(vals), vals = meta.truth(p) + [-1;1]; end
    lo = prctile(vals, 5); hi = prctile(vals, 95);
    pad = 0.12*max(1e-9, hi-lo);
    yl = [lo-pad, hi+pad];
    if p==1, yl = [0 1]; end                 % π_1 in [0,1]
    if any(p == [3,5])                       % sigmas > 0
        yl(1) = max(1e-6, yl(1));
    end
    ylim(yl);

    if p==1, legend('Location','best'); end
end
sgtitle('PL Mixture (K=2) with Outliers: Average Estimates vs \epsilon');

end

% ======================================================================
% =============================== HELPERS ===============================
% ======================================================================

function Y = sim_pl_mix(n, piK, mu, sigma)
% Simulate Poisson–Lognormal mixture: pick component, draw lambda~LogN, then Y~Pois(lambda)
K = numel(piK);
u = rand(n,1); edges = [0; cumsum(piK(:))];
z = zeros(n,1);
for k = 1:K
    z = z + (u>edges(k) & u<=edges(k+1))*k;
end
Y = zeros(n,1);
for k = 1:K
    idx = (z==k);
    if any(idx)
        lam = lognrnd(mu(k), sigma(k), sum(idx), 1);
        Y(idx) = poissrnd(lam);
    end
end
end

function fit = align_pl_by_mean(fit)
% Order components so that component 1 has the *larger* mean: E[lambda]=exp(mu+0.5*sigma^2)
mu = fit.theta(1,:); sg = fit.theta(2,:);
m  = exp(mu + 0.5*sg.^2);
[~,ord] = sort(m,'descend');
fit.theta = fit.theta(:,ord);
fit.pi    = fit.pi(ord);
end

% ---------------- EM for PL (known K) ----------------
function fit = em_pl_knownK(Y, K, nStarts, MaxIter, Tol, init, warm, GH)
Y = Y(:); n = numel(Y);
best = struct('ll',-inf,'pi',[],'theta',[],'iters',0,'converged',false);

% Build start set: warm (if provided) + (nStarts-1) fresh
starts = cell(nStarts,1); ns = 0;
if nargin>=8 && ~isempty(warm)
    ns = ns+1; starts{ns} = struct('pi',warm.pi(:), 'theta',warm.theta);
end
for t = (ns+1):nStarts
    [pi0, th0] = init_pl(Y, K, init);
    starts{t} = struct('pi',pi0, 'theta',th0);
end

for s0 = 1:nStarts
    piK = starts{s0}.pi;  th = starts{s0}.theta;
    [ll, piK, th, iters, conv] = em_core_pl(Y, piK, th, MaxIter, Tol, GH);
    if ll > best.ll
        best = struct('ll',ll,'pi',piK,'theta',th,'iters',iters,'converged',conv);
    end
end
fit = struct('pi',best.pi,'theta',best.theta,'loglik',best.ll, ...
             'iters',best.iters,'converged',best.converged);
end

function [ll, piK, theta, iters, converged] = em_core_pl(Y, piK, theta, MaxIter, Tol, GH)
% EM using GH-approximated PL pmf; M-step for (mu,log sigma) via fminunc
n = numel(Y); K = numel(piK); converged=false; prev_ll=-inf;

log_y_fact = gammaln(Y+1);  % for Poisson pmf stabilization
opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                    'MaxIterations',200,'MaxFunctionEvaluations',2e4);

for it = 1:MaxIter
    % ---- E-step: compute Pk(y_i | theta_k) via GH ----
    Pk = zeros(n,K);
    for k = 1:K
        mu = theta(1,k); sg = max(1e-8, theta(2,k));
        Pk(:,k) = pl_pmf_gh_Y(Y, log_y_fact, mu, sg, GH);
    end
    Pk = max(Pk, 1e-300);
    mix = max(Pk*piK(:), 1e-300);
    tau = bsxfun(@rdivide, bsxfun(@times, Pk, piK(:).'), mix);
    ll  = sum(log(mix));

    % ---- M-step: weights ----
    Nk  = sum(tau,1)'; piK = max(eps, Nk/n);

    % ---- M-step: (mu, sigma) via weighted GH-likelihood ----
    th = zeros(2,K);
    for k = 1:K
        wk = tau(:,k);
        mu0 = theta(1,k); sg0 = max(1e-6, theta(2,k));
        v0  = [mu0; log(sg0)];
        obj = @(v) neg_wPL_ll(Y, log_y_fact, wk, v, GH); % numeric grad
        vopt = fminunc(obj, v0, opts);
        th(:,k) = [vopt(1); max(1e-8, exp(vopt(2)))];
    end
    theta = th;

    if it>1 && abs(ll - prev_ll) < Tol*max(1,abs(prev_ll)), converged=true; break; end
    prev_ll = ll;
end
iters = it;
end

function f = neg_wPL_ll(Y, log_y_fact, w, v, GH)
% Weighted negative log-likelihood for PL (GH-approximated)
mu = v(1); sg = max(1e-8, exp(v(2)));
p  = pl_pmf_gh_Y(Y, log_y_fact, mu, sg, GH);
f  = -sum( w .* log( max(p,1e-300) ) );
end

function [pi0, theta0] = init_pl(Y, K, how)
switch lower(how)
    case 'random'
        a = rand(K,1)+0.1; pi0 = a/sum(a);
    otherwise
        pi0 = ones(K,1)/K;    % 'quantile' weights
end
% crude PL starts from quantiles on log(Y+0.5)
qy = prctile(Y+0.5, linspace(5,95,K));
mu0 = log(max(0.1, qy(:)));
sg0 = 0.5*ones(K,1);
theta0 = [mu0.'; sg0.'];
end

% ---------------- HMIX / VNED for PL (known K) ----------------
function fit = dic_pl_knownK(Y, K, method, nStarts, MaxIter, Tol, init, warm, GH)
[YY, g] = empirical(Y); g = max(g,1e-12);
best = struct('obj',Inf,'pi',[],'theta',[]);

% start set: warm + fresh
starts = cell(nStarts,1); ns=0;
if nargin>=9 && ~isempty(warm)
    ns = ns+1; starts{ns} = struct('pi',warm.pi(:), 'theta',warm.theta);
end
for t = (ns+1):nStarts
    [pi0, th0] = init_pl(YY, K, init);
    starts{t} = struct('pi',pi0, 'theta',th0);
end

opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                    'MaxIterations',200,'MaxFunctionEvaluations',2e4);

for s0 = 1:nStarts
    piK = starts{s0}.pi; theta = starts{s0}.theta;
    prev = [piK; theta(:)];
    iters = 0;

    while iters < MaxIter
        iters = iters + 1;
        % E on support
        Pk  = comp_pmf_matrix_pl(YY, theta, GH);     % |support| x K
        num = bsxfun(@times, Pk, piK(:)');
        den = max(sum(num,2), 1e-15);
        R   = bsxfun(@rdivide, num, den);
        R_safe = max(R, 1e-12); g_safe = max(g, 1e-12);

        theta_new = theta;

        switch lower(method)
            case 'hd'    % Hellinger step
                fvals = zeros(K,1);
                for k = 1:K
                    rk = R(:,k);
                    mu0 = theta(1,k); sg0 = max(1e-6, theta(2,k));
                    Hk = @(u) -sum( sqrt( max(1e-12,rk).*g_safe .* ...
                          pl_pmf_gh(YY, u(1), max(1e-8,exp(u(2))), GH) ) );
                    [uk,fk] = fminunc(Hk, [mu0; log(sg0)], opts);
                    theta_new(:,k) = [uk(1); max(1e-8,exp(uk(2)))];
                    fvals(k) = fk;
                end
                w = fvals.^2; if ~all(isfinite(w)), w = ones(K,1); end
                pi_new = max(eps, w/sum(w));

            otherwise   % VNED step
                vned_k = zeros(K,1);
                for k = 1:K
                    rk = R(:,k); rk_safe = R_safe(:,k);
                    mu0 = theta(1,k); sg0 = max(1e-6, theta(2,k));
                    NEDk = @(u) sum( exp( - piK(k).*pl_pmf_gh(YY, u(1), max(1e-8,exp(u(2))), GH) ...
                                        ./ (g_safe.*rk_safe) ) .* g_safe .* rk );
                    uk  = fminunc(NEDk, [mu0; log(sg0)], opts);
                    theta_new(:,k) = [uk(1); max(1e-8,exp(uk(2)))];
                    pmf = pl_pmf_gh(YY, theta_new(1,k), theta_new(2,k), GH);
                    vned_k(k) = sum( exp( - piK(k).*pmf./(g_safe.*rk_safe) ) .* (piK(k).*pmf) );
                    if ~isfinite(vned_k(k)), vned_k(k) = 1e-12; end
                end
                pi_new = max(eps, vned_k/sum(vned_k));
        end

        if max(abs([pi_new; theta_new(:)] - prev)) <= Tol
            piK = pi_new; theta = theta_new; break;
        end
        piK = pi_new; theta = theta_new; prev = [piK; theta(:)];
    end

    % divergence value for selection among starts
    f  = max(mix_pmf_pl(YY, piK, theta, GH), 1e-15);
    if strcmpi(method,'hd')
        D = 2*(1 - sum(sqrt(g.*f)));
    else
        D = sum( exp(-f./g).*g );
    end
    if D < best.obj
        best.obj = D; best.pi = piK; best.theta = theta;
    end
end

fit = struct('pi',best.pi,'theta',best.theta);
end

% ---------- PL pmf via Gauss–Hermite quadrature ----------
function p = pl_pmf_gh_Y(Y, log_y_fact, mu, sigma, GH)
% Vectorized over Y (Nx1). Uses:
%   ∫ Poiss(y|λ) LogN(λ|μ,σ) dλ
% = E[ Poiss(y|exp(μ+σZ)) ],  Z~N(0,1)
% ≈ (1/√π) Σ w_i Poiss(y | exp(μ + σ√2 x_i))
N = numel(Y);
p = zeros(N,1);
scale = GH.scale; s2 = GH.s2; x = GH.x; w = GH.w;

for i = 1:numel(x)
    lam = exp(mu + sigma*s2*x(i));
    % log Poisson pmf: y*log(lam) - lam - log(y!)
    lp = Y.*log(lam) - lam - log_y_fact;
    p  = p + w(i)*exp(lp);
end
p = scale * p;
p = max(p, 1e-300);
end

function p = pl_pmf_gh(YY, mu, sigma, GH)
% Same as above but for a vector of support counts YY (Mx1), without reusing log_y_fact.
log_y_fact = gammaln(YY+1);
p = pl_pmf_gh_Y(YY, log_y_fact, mu, sigma, GH);
end

function P = comp_pmf_matrix_pl(YY, theta, GH)
K = size(theta,2);
P = zeros(numel(YY), K);
for k = 1:K
    mu = theta(1,k); sg = max(1e-8, theta(2,k));
    P(:,k) = pl_pmf_gh(YY, mu, sg, GH);
end
P = max(P, 1e-15);
end

function f = mix_pmf_pl(YY, piK, theta, GH)
P = comp_pmf_matrix_pl(YY, theta, GH);
f = max(P * piK(:), 1e-15);
end

% ---------- GH rule (Golub–Welsch) ----------
function GH = gh_rule(n)
% Nodes/weights for ∫_{-∞}^{∞} e^{-x^2} f(x) dx ≈ Σ w_i f(x_i)
% Map for N(0,1): ∫ f(z) φ(z) dz = (1/√π) Σ w_i f(√2 x_i)
i  = (1:n-1)';
a  = zeros(n,1);
b  = sqrt(i/2);
J  = diag(a) + diag(b,1) + diag(b,-1);
[V,D] = eig(J);
[x,idx] = sort(diag(D));
V = V(:,idx);
w = sqrt(pi) * (V(1,:)').^2;

GH.x = x;
GH.w = w;
GH.scale = 1/sqrt(pi);   % factor for φ(z) mapping
GH.s2    = sqrt(2);      % z = √2 * x
end

% ---------- empirical distribution on counts ----------
function [YY, g] = empirical(Y)
[YY,~,ic] = unique(Y);
cnt = accumarray(ic,1);
g   = cnt / numel(Y);
end
