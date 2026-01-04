function [avg, est, meta] = demo_pg_outlier_sim()
% DEMO_PG_OUTLIER_SIM
% Monte Carlo for PG (Poisson–Gamma) mixture with *known K=2* and
% ε-fraction outliers set to a fixed large count (default 50).
% Methods: EM (MLE), HMIX (Hellinger divergence), VNEDMIX (VNED).
% Returns:
%   avg.(EM|HD|VNED)  -> [E x P] averages across reps (NaN-robust)
%   est.(EM|HD|VNED)  -> [R x E x P] all estimates (per rep & epsilon)
%   meta              -> struct with cont_grid, truth, etc.
%
% Parameter order P = [pi1, a1, b1, a2, b2]  (means = a./b)

% ------------------ Simulation settings ------------------
n            = 5000;                            % sample size per dataset
R            = 5000;                            % Monte Carlo reps
cont_grid    = [0 0.05 0.10 0.15 0.20 0.30];    % ε levels (outlier pct)
outlier_val  = 50;                              % value used as outlier
nStarts0     = 6;                               % multistarts at ε=0
nStartsWarm  = 1;                               % warm-starts for next ε
MaxIter      = 150;                             % inner EM/optimizer iters
Tol          = 1e-6;                            % convergence tolerance
init         = 'quantile';                      % 'quantile' | 'random'
rng(1);                                         % base RNG

% Truth (component 1 has larger mean for stable alignment)
pi_true = [0.3; 0.7];
a_true  = [10; 1];
b_true  = [ 1; 2];
K       = 2;

% ------------------ Storage ------------------
E = numel(cont_grid);
P = 5;  % [pi1, a1, b1, a2, b2]
est_EM   = NaN(R,E,P);
est_HD   = NaN(R,E,P);
est_VNED = NaN(R,E,P);

% ------------------ Parallel pool ------------------
if license('test','Distrib_Computing_Toolbox')
    p = gcp('nocreate'); if isempty(p), parpool('local'); end %#ok<PPOOL>
end

% ------------------ Monte Carlo (parallel over reps) ------------------
parfor r = 1:R
    rng(1000 + r); % independent stream
    local_EM   = NaN(E,P);
    local_HD   = NaN(E,P);
    local_VNED = NaN(E,P);

    % warm starts carried across epsilon levels for this replication
    prev_em   = [];
    prev_hd   = [];
    prev_vned = [];

    cg = cont_grid; % avoid broadcast overhead warnings

    for e = 1:E
        cont = cg(e);

        % --- Generate clean PG mixture, then inject outliers ---
        Y = sim_pg_mix(n, pi_true, a_true, b_true);
        m = round(cont*n);
        if m>0
            idx = randperm(n, m);
            Y(idx) = outlier_val;
        end

        % --- EM (known K) ---
        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_em = em_pg_knownK(Y, K, nStarts, MaxIter, Tol, init, prev_em);
        if fit_em.converged
            fit_em = align_pg_by_mean(fit_em);
            local_EM(e,:) = [fit_em.pi(1), fit_em.theta(1,1), fit_em.theta(2,1), ...
                                         fit_em.theta(1,2), fit_em.theta(2,2)];
            prev_em = fit_em; % warm start forward
        else
            prev_em = [];     % drop warm start if it failed
        end

        % --- HMIX (HD) ---
        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_hd = dic_pg_knownK(Y, K, 'hd', nStarts, MaxIter, Tol, init, prev_hd);
        fit_hd = align_pg_by_mean(fit_hd);
        local_HD(e,:) = [fit_hd.pi(1), fit_hd.theta(1,1), fit_hd.theta(2,1), ...
                                      fit_hd.theta(1,2), fit_hd.theta(2,2)];
        prev_hd = fit_hd;

        % --- VNEDMIX (VNED) ---
        nStarts = (e==1)*nStarts0 + (e>1)*nStartsWarm;
        fit_vn = dic_pg_knownK(Y, K, 'vned', nStarts, MaxIter, Tol, init, prev_vned);
        fit_vn = align_pg_by_mean(fit_vn);
        local_VNED(e,:) = [fit_vn.pi(1), fit_vn.theta(1,1), fit_vn.theta(2,1), ...
                                        fit_vn.theta(1,2), fit_vn.theta(2,2)];
        prev_vned = fit_vn;
    end

    est_EM(r,:,:)   = local_EM;
    est_HD(r,:,:)   = local_HD;
    est_VNED(r,:,:) = local_VNED;
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
meta.param_names = {'\pi_1','\alpha_1','\beta_1','\alpha_2','\beta_2'};
meta.truth       = [pi_true(1), a_true(1), b_true(1), a_true(2), b_true(2)];
meta.n = n; meta.R = R; meta.outlier_val = outlier_val;

% ------------------ Plot (robust y-limits + truth as dotted) ------------------
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
    if p~=1, yl(1) = max(1e-6, yl(1)); end   % positive shape/rate
    ylim(yl);

    if p==1, legend('Location','best'); end
end
sgtitle('PG Mixture (K=2) with Outliers: Average Estimates vs \epsilon');

% Optionally save:
% save('pg_outlier_results.mat','avg','est','meta');

end

% ======================================================================
% =============================== HELPERS ===============================
% ======================================================================

function Y = sim_pg_mix(n, piK, a, b)
% Simulate PG mixture via its NB marginal: r=a, p=b/(1+b)
K = numel(piK);
u = rand(n,1); edges = [0; cumsum(piK(:))];
z = zeros(n,1);
for k = 1:K, z = z + (u>edges(k) & u<=edges(k+1))*k; end
Y = zeros(n,1);
for k = 1:K
    idx = (z==k);
    if any(idx)
        p = b(k)/(1+b(k));
        Y(idx) = nbinrnd(a(k), p, sum(idx), 1);
    end
end
end

function fit = align_pg_by_mean(fit)
% Order components so that component 1 has the *larger* mean (a/b)
mu = fit.theta(1,:)./fit.theta(2,:);
[~,ord] = sort(mu,'descend');
fit.theta = fit.theta(:,ord);
fit.pi    = fit.pi(ord);
end

% ---------------- EM for PG (known K) ----------------
function fit = em_pg_knownK(Y, K, nStarts, MaxIter, Tol, init, warm)
Y = Y(:); n = numel(Y);
best = struct('ll',-inf,'pi',[],'theta',[],'iters',0,'converged',false);

% Build start set: warm (if provided) + (nStarts-1) fresh
starts = cell(nStarts,1); ns = 0;
if nargin>=7 && ~isempty(warm)
    ns = ns+1; starts{ns} = struct('pi',warm.pi(:), 'theta',warm.theta);
end
for t = (ns+1):nStarts
    [pi0, th0] = init_pg(Y, K, init);
    starts{t} = struct('pi',pi0, 'theta',th0);
end

for s0 = 1:nStarts
    piK = starts{s0}.pi;  th = starts{s0}.theta;
    [ll, piK, th, iters, conv] = em_core_pg_fast(Y, piK, th, MaxIter, Tol);
    if ll > best.ll
        best = struct('ll',ll,'pi',piK,'theta',th,'iters',iters,'converged',conv);
    end
end
fit = struct('pi',best.pi,'theta',best.theta,'loglik',best.ll, ...
             'iters',best.iters,'converged',best.converged);
end

function [ll, piK, theta, iters, converged] = em_core_pg_fast(Y, piK, theta, MaxIter, Tol)
% EM using log-sum-exp and gradient-based weighted NB MLE in the M-step.
n = numel(Y); K = numel(piK); converged=false; prev_ll=-inf;

opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                    'MaxIterations',200,'MaxFunctionEvaluations',2e4, ...
                    'OptimalityTolerance',1e-6,'StepTolerance',1e-8, ...
                    'SpecifyObjectiveGradient',true);

for it = 1:MaxIter
    % ---- E-step (log-stabilized) ----
    L = zeros(n,K);
    for k = 1:K
        a = max(1e-8, theta(1,k)); b = max(1e-8, theta(2,k));
        p = b/(1+b);
        L(:,k) = nb_logpmf(Y,a,p) + log(max(1e-300,piK(k)));
    end
    [mL, ~] = max(L,[],2);
    Ls = L - mL;             % subtract row max
    mix = mL + logsumexp_rows(Ls);
    tau = exp(Ls - logsumexp_rows(Ls));  % responsibilities (rows sum to 1)
    ll  = sum(mix);

    % ---- M-step: weights ----
    Nk  = sum(tau,1)'; piK = max(eps, Nk/n);

    % ---- M-step: PG params via weighted NB MLE (BFGS with gradient) ----
    th = zeros(2,K);
    for k = 1:K
        wk = tau(:,k);
        a0 = max(1e-3, theta(1,k)); b0 = max(1e-3, theta(2,k));
        u0 = [log(a0); log(b0)];
        obj = @(u) neg_wNB_ll_and_grad(Y, wk, u);
        uopt = fminunc(obj, u0, opts);
        th(:,k) = [max(1e-8,exp(uopt(1))); max(1e-8,exp(uopt(2)))];
    end
    theta = th;

    if it>1 && abs(ll - prev_ll) < Tol*max(1,abs(prev_ll)), converged=true; break; end
    prev_ll = ll;
end
iters = it;
end

function [f, g] = neg_wNB_ll_and_grad(Y, w, u)
% Weighted NB negative loglik and gradient wrt u = [log a; log b], p=b/(1+b).
a = max(1e-8, exp(u(1))); 
b = max(1e-8, exp(u(2)));
p = b/(1+b);            % p in (0,1)

% log pmf and its derivatives
lg = gammaln(Y+a) - gammaln(a) - gammaln(Y+1) + a*log(1-p) + Y.*log(p);
f  = -sum(w .* lg);

% gradients (chain rule: dp/db = 1/(1+b)^2)
psi_a = psi(a); psi_ya = psi(Y+a);
dL_da = psi_ya - psi_a + log(1-p);               % ∂ log pmf / ∂a
dL_dp = a/(p-1) + Y./p;                           % ∂ log pmf / ∂p
dp_db = 1/(1+b)^2;
dL_db = dL_dp * dp_db;                            % ∂ log pmf / ∂b
g_a   = -sum(w .* dL_da) * a;                     % chain: da/du1 = a
g_b   = -sum(w .* dL_db) * b;                     % chain: db/du2 = b
g = [g_a; g_b];
end

function s = logsumexp_rows(A)
% row-wise log(sum(exp(A),2))
m = max(A,[],2);
s = m + log(sum(exp(A - m),2));
end

function v = nb_logpmf(y,a,p)
% log NB PMF: Y ~ NB(r=a, p), support {0,1,...}
v = gammaln(y+a) - gammaln(a) - gammaln(y+1) + a*log(1-p) + y.*log(p);
end

function [pi0, theta0] = init_pg(Y, K, how)
switch lower(how)
    case 'random'
        a = rand(K,1)+0.1; pi0 = a/sum(a);
    otherwise
        pi0 = ones(K,1)/K;    % 'quantile' weights
end
% crude PG starts from quantiles/variance
qs = linspace(5,95,K); cuts = prctile(Y, qs);
if any(~isfinite(cuts)), cuts = linspace(min(Y), max(Y), K); end
mu_guess = max(1e-2, cuts(:));
vY = max(1, var(Y));
b_guess = max(0.1, mu_guess./max(mu_guess.^2, vY)); % stabilize
a_guess = max(0.2, mu_guess .* b_guess);
theta0  = [a_guess.'; b_guess.'];
end

% ---------------- HMIX / VNED for PG (known K) ----------------
function fit = dic_pg_knownK(Y, K, method, nStarts, MaxIter, Tol, init, warm)
[YY, g] = empirical(Y); g = max(g,1e-12);
best = struct('obj',Inf,'pi',[],'theta',[]);

% start set: warm + fresh
starts = cell(nStarts,1); ns=0;
if nargin>=8 && ~isempty(warm)
    ns = ns+1; starts{ns} = struct('pi',warm.pi(:), 'theta',warm.theta);
end
for t = (ns+1):nStarts
    [pi0, th0] = init_pg(YY, K, init);
    starts{t} = struct('pi',pi0, 'theta',th0);
end

for s0 = 1:nStarts
    f1 = fit_one_dic_pg(YY, g, K, method, MaxIter, Tol, starts{s0});
    f  = max(mix_pmf_pg(YY, f1.pi, f1.theta), 1e-15);
    if strcmpi(method,'hd')
        D = 2*(1 - sum(sqrt(g.*f)));
    else
        D = sum( exp(-f./g).*g );
    end
    if D < best.obj
        best.obj = D; best.pi = f1.pi; best.theta = f1.theta;
    end
end
fit = struct('pi',best.pi,'theta',best.theta);
end

function fit = fit_one_dic_pg(YY, g, K, method, MaxIter, Tol, start)
piK = start.pi(:); theta = start.theta;
prev = [piK; theta(:)]; iters = 0;

opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                    'MaxIterations',200,'MaxFunctionEvaluations',2e4);

while iters < MaxIter
    iters = iters + 1;
    % E on support
    Pk  = comp_pmf_matrix_pg(YY, theta);     % |support| x K
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
                a0 = max(1e-8, theta(1,k)); b0 = max(1e-8, theta(2,k));
                Hk = @(u) -sum( sqrt( max(1e-12,rk).*g_safe .* ...
                          nbinpdf(YY, max(1e-12,exp(u(1))), ...
                                       max(1e-12,exp(u(2)))/(1+max(1e-12,exp(u(2))))) ) );
                [uk,fk] = fminunc(Hk, [log(a0); log(b0)], opts);
                theta_new(:,k) = [max(1e-8,exp(uk(1))); max(1e-8,exp(uk(2)))];
                fvals(k) = fk;
            end
            w = fvals.^2; if ~all(isfinite(w)), w = ones(K,1); end
            pi_new = max(eps, w/sum(w));

        otherwise   % VNED step
            vned_k = zeros(K,1);
            for k = 1:K
                rk = R(:,k); rk_safe = R_safe(:,k);
                a0 = max(1e-8, theta(1,k)); b0 = max(1e-8, theta(2,k));
                NEDk = @(u) sum( exp( - piK(k).*nbinpdf(YY, max(1e-12,exp(u(1))), ...
                              max(1e-12,exp(u(2)))/(1+max(1e-12,exp(u(2))))) ./ (g_safe.*rk_safe) ) .* g_safe .* rk );
                uk  = fminunc(NEDk, [log(a0); log(b0)], opts);
                theta_new(:,k) = [max(1e-8,exp(uk(1))); max(1e-8,exp(uk(2)))];
                pmf = nbinpdf(YY, theta_new(1,k), theta_new(2,k)/(1+theta_new(2,k)));
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
fit = struct('pi',piK,'theta',theta,'iters',iters);
end

function P = comp_pmf_matrix_pg(YY, theta)
K = size(theta,2);
P = zeros(numel(YY), K);
for k = 1:K
    a = theta(1,k); b = theta(2,k);
    P(:,k) = nbinpdf(YY, a, b/(1+b));
end
P = max(P, 1e-15);
end

function f = mix_pmf_pg(YY, piK, theta)
P = comp_pmf_matrix_pg(YY, theta);
f = max(P * piK(:), 1e-15);
end

function [YY, g] = empirical(Y)
tbl = tabulate(Y);
YY = tbl(:,1);
g  = tbl(:,3)/100;
end
