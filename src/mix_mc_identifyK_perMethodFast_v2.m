function out = mix_mc_identifyK_perMethodFast_v2(B, n, model, trueParams, varargin)
% MIX_MC_IDENTIFYK_PERMETHODFAST_V2 (Self-contained; PL only; R2020a)
% Monte Carlo across B datasets. For each method (EM/HD/VNED):
%   (1) Select K^ on a subsample (fast settings; K=1..Kmax).
%   (2) Refit once on full data at that method's K^ (fit settings).
%
% Outputs:
%   out.khat.(EM|HD|VNED)         -> [B x 1] selected K per dataset
%   out.sel_counts.(EM|HD|VNED){b}-> struct with fields K and count
%   out.sel_frac_true.(EM|HD|VNED)-> fraction of splits selecting true K
%   out.est.(EM|HD|VNED){b}       -> struct {K, pi, theta} from full refit
%   out.avg_overall.(method)      -> NaN-robust averages across all datasets
%   out.avg_when_correct.(method) -> averages where that method's K^==trueK
%   out.prop_correct.(method)     -> % datasets where that method's K^==trueK
%   out.settings                  -> echo of settings
%
% NOTE: This self-contained version supports model='pl' only.

% ---------------- Parse inputs ----------------
ip = inputParser;
ip.addRequired('B', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('n', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('model', @(s)(ischar(s)||isstring(s))&&ismember(lower(char(s)),{'pl'}));
ip.addRequired('trueParams', @(s)isstruct(s)&&isfield(s,'pi')&&isfield(s,'theta'));

ip.addParameter('methods', {'em','hd','vned'}, @(c)iscellstr(c)&&~isempty(c));
ip.addParameter('Kmax',   5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% Selection on subsample (fast)
ip.addParameter('n_sub', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>=100);
ip.addParameter('sel_C', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('sel_nStarts', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('sel_MaxIter', 60, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
ip.addParameter('sel_Tol', 1e-5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('sel_rule','min', @(s)ischar(s)&&ismember(lower(s),{'min','stop'})); % we implement 'min'

% Full-data refits at fixed K (per method)
ip.addParameter('fit_C', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('fit_nStarts', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('fit_MaxIter', 120, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
ip.addParameter('fit_Tol', 1e-5, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% PL Gauss–Hermite nodes
ip.addParameter('gh_n', 12, @(x)isnumeric(x)&&isscalar(x)&&x>=8);

% RNG + parallel + verbosity
ip.addParameter('seed', 1, @(x)isnumeric(x)&&isscalar(x));
ip.addParameter('UseParallel', false, @(b)islogical(b)||ismember(b,[0 1]));
ip.addParameter('Verbose', true, @(b)islogical(b)||ismember(b,[0 1]));
ip.parse(B,n,model,trueParams,varargin{:});
S = ip.Results;

% Normalize
model = lower(char(S.model));
methods = lower(S.methods(:)');

% Shape checks (PL only)
Ktrue = numel(S.trueParams.pi(:));
assert(size(S.trueParams.theta,1)==2, 'PL: theta must be 2xK ([mu;sigma]).');
assert(size(S.trueParams.theta,2)==Ktrue, 'PL: size(theta,2)==numel(pi)==K.');

% ---------------- Prealloc sliced outputs ----------------
khat_EM = NaN(B,1); khat_HD = NaN(B,1); khat_VN = NaN(B,1);
est_EM  = cell(B,1); est_HD  = cell(B,1); est_VN  = cell(B,1);
sel_counts_EM = cell(B,1); sel_counts_HD = cell(B,1); sel_counts_VN = cell(B,1);
sel_frac_EM = NaN(B,1); sel_frac_HD = NaN(B,1); sel_frac_VN = NaN(B,1);

% ---------------- Pool ----------------
UsePar = S.UseParallel;
if UsePar
    try
        if license('test','Distrib_Computing_Toolbox')
            p = gcp('nocreate'); if isempty(p), parpool('local'); end
        else
            warning('Parallel toolbox unavailable; running serial.'); UsePar=false;
        end
    catch
        UsePar=false;
    end
end

% ---------------- MC loop ----------------
if UsePar
    parfor b = 1:B
        rng(S.seed + b);
        [khat_EM(b), est_EM{b}, sel_counts_EM{b}, sel_frac_EM(b), ...
         khat_HD(b), est_HD{b}, sel_counts_HD{b}, sel_frac_HD(b), ...
         khat_VN(b), est_VN{b}, sel_counts_VN{b}, sel_frac_VN(b)] = ...
            one_dataset_PL(b, S, methods, Ktrue);
    end
else
    for b = 1:B
        rng(S.seed + b);
        [khat_EM(b), est_EM{b}, sel_counts_EM{b}, sel_frac_EM(b), ...
         khat_HD(b), est_HD{b}, sel_counts_HD{b}, sel_frac_HD(b), ...
         khat_VN(b), est_VN{b}, sel_counts_VN{b}, sel_frac_VN(b)] = ...
            one_dataset_PL(b, S, methods, Ktrue);
    end
end

% ---------------- Summaries ----------------
prop_correct = struct();
prop_correct.EM   = 100*mean(khat_EM==Ktrue,'omitnan');
prop_correct.HD   = 100*mean(khat_HD==Ktrue,'omitnan');
prop_correct.VNED = 100*mean(khat_VN==Ktrue,'omitnan');

avg_overall = struct();
avg_overall.EM   = avg_from_cells(est_EM,  'pl');
avg_overall.HD   = avg_from_cells(est_HD,  'pl');
avg_overall.VNED = avg_from_cells(est_VN,  'pl');

avg_when_correct = struct();
avg_when_correct.EM   = avg_from_cells(est_EM(khat_EM==Ktrue),   'pl');
avg_when_correct.HD   = avg_from_cells(est_HD(khat_HD==Ktrue),   'pl');
avg_when_correct.VNED = avg_from_cells(est_VN(khat_VN==Ktrue),   'pl');

% ---------------- Pack ----------------
out = struct();
out.khat = struct('EM',khat_EM,'HD',khat_HD,'VNED',khat_VN);
out.sel_counts = struct('EM',sel_counts_EM,'HD',sel_counts_HD,'VNED',sel_counts_VN);
out.sel_frac_true = struct('EM',sel_frac_EM,'HD',sel_frac_HD,'VNED',sel_frac_VN);
out.est = struct('EM',est_EM,'HD',est_HD,'VNED',est_VN);

out.avg_overall       = avg_overall;
out.avg_when_correct  = avg_when_correct;
out.prop_correct      = prop_correct;

out.settings = S; out.settings.trueK = Ktrue;
end

% ======================================================================
% ======================= One dataset (PL only) ========================
% ======================================================================
function [K_em, estEM, C_em, F_em, ...
          K_hd, estHD, C_hd, F_hd, ...
          K_vn, estVN, C_vn, F_vn] = one_dataset_PL(b, S, methods, Ktrue)

% simulate full data
Yfull = generate_dataset_PL(S.n, S.trueParams);

% selection splits
n_sub = min(S.n, S.n_sub);
Kmax  = S.Kmax;
C     = S.sel_C;

% containers
Ksel_EM  = NaN(C,1);  Ksel_HD  = NaN(C,1);  Ksel_VN  = NaN(C,1);

for c = 1:C
    if n_sub < S.n
        idx  = randperm(S.n, n_sub);
        Ysel = Yfull(idx);
    else
        Ysel = Yfull;
    end

    GH = gh_rule(S.gh_n);

    % --- EM selection by BIC (with multi-starts) ---
    if any(strcmp(methods,'em'))
        bestScore = Inf; bestK = NaN;
        for K = 1:Kmax
            ll_best = -Inf;
            for s = 1:max(1,S.sel_nStarts)
                try
                    start_init = tern(s==1,'kmeans','random'); % stronger starts
                    fit = em_pl_knownK(Ysel, K, 1, S.sel_MaxIter, S.sel_Tol, start_init, [], GH);
                    if isfield(fit,'loglik') && isfinite(fit.loglik) && fit.loglik > ll_best
                        ll_best = fit.loglik;
                    end
                catch
                end
            end
            if isfinite(ll_best)
                d = 3*K - 1;  % (K-1) weights + 2K PL params
                bic = -2*ll_best + d*log(numel(Ysel));
                if bic < bestScore
                    bestScore = bic; bestK = K;
                end
            end
        end
        Ksel_EM(c) = bestK;
    end

    % --- HMIX selection by min divergence (multi-starts) ---
    if any(strcmp(methods,'hd'))
        bestD = Inf; bestK = NaN;
        for K = 1:Kmax
            Dbest = Inf;
            for s = 1:max(1,S.sel_nStarts)
                try
                    start_init = tern(s==1,'kmeans','random');
                    fit_hd = dic_pl_knownK(Ysel, K, 'hd', 1, S.sel_MaxIter, S.sel_Tol, start_init, [], GH);
                    if isfield(fit_hd,'obj') && isfinite(fit_hd.obj) && fit_hd.obj < Dbest
                        Dbest = fit_hd.obj;
                    end
                catch
                end
            end
            if isfinite(Dbest) && Dbest < bestD
                bestD = Dbest; bestK = K;
            end
        end
        Ksel_HD(c) = bestK;
    end

    % --- VNED selection by min divergence (multi-starts) ---
    if any(strcmp(methods,'vned'))
        bestD = Inf; bestK = NaN;
        for K = 1:Kmax
            Dbest = Inf;
            for s = 1:max(1,S.sel_nStarts)
                try
                    start_init = tern(s==1,'kmeans','random');
                    fit_vn = dic_pl_knownK(Ysel, K, 'vned', 1, S.sel_MaxIter, S.sel_Tol, start_init, [], GH);
                    if isfield(fit_vn,'obj') && isfinite(fit_vn.obj) && fit_vn.obj < Dbest
                        Dbest = fit_vn.obj;
                    end
                catch
                end
            end
            if isfinite(Dbest) && Dbest < bestD
                bestD = Dbest; bestK = K;
            end
        end
        Ksel_VN(c) = bestK;
    end
end

% summarize selection
[K_em, C_em, F_em]   = summarize_selection(Ksel_EM, Ktrue);
[K_hd, C_hd, F_hd]   = summarize_selection(Ksel_HD, Ktrue);
[K_vn, C_vn, F_vn]   = summarize_selection(Ksel_VN, Ktrue);

% ---------------- Refit on FULL data at K^ ----------------
estEM=[]; estHD=[]; estVN=[];
GHfull = gh_rule(S.gh_n);

if ~isnan(K_em) && any(strcmp(methods,'em'))
    try
        fit = em_pl_knownK(Yfull, K_em, S.fit_nStarts, S.fit_MaxIter, S.fit_Tol, 'kmeans', [], GHfull);
        fit = align_pl_by_mean(fit);
        estEM = struct('K',K_em, 'pi',fit.pi(:), 'theta',fit.theta);
    catch
    end
end
if ~isnan(K_hd) && any(strcmp(methods,'hd'))
    try
        fit_hd = dic_pl_knownK(Yfull, K_hd, 'hd', S.fit_nStarts, S.fit_MaxIter, S.fit_Tol, 'kmeans', [], GHfull);
        fit_hd = align_pl_by_mean(fit_hd);
        estHD  = struct('K',K_hd, 'pi',fit_hd.pi(:), 'theta',fit_hd.theta);
    catch
    end
end
if ~isnan(K_vn) && any(strcmp(methods,'vned'))
    try
        fit_vn = dic_pl_knownK(Yfull, K_vn, 'vned', S.fit_nStarts, S.fit_MaxIter, S.fit_Tol, 'kmeans', [], GHfull);
        fit_vn = align_pl_by_mean(fit_vn);
        estVN  = struct('K',K_vn, 'pi',fit_vn.pi(:), 'theta',fit_vn.theta);
    catch
    end
end
end

% ---------------- selection summary helper ----------------
function [Khat, Kcnt, fracTrue] = summarize_selection(Ks, Ktrue)
Ks = Ks(~isnan(Ks));
if isempty(Ks)
    Khat = NaN; Kcnt = struct('K',[],'count',[]); fracTrue = NaN; return;
end
uK = unique(Ks(:)');
cnt = zeros(size(uK));
for i=1:numel(uK), cnt(i) = sum(Ks==uK(i)); end
[~,ix] = max(cnt); Khat = uK(ix);
Kcnt = struct('K',uK(:),'count',cnt(:));
fracTrue = 0; j = find(uK==Ktrue,1); if ~isempty(j), fracTrue = cnt(j)/numel(Ks); end
end

% ======================================================================
% =============================== HELPERS ===============================
% ======================================================================
function Y = generate_dataset_PL(n, truth)
piK = truth.pi(:); K = numel(piK);
u = rand(n,1); edges = [0; cumsum(piK(:))];
z = zeros(n,1);
for k = 1:K, z = z + (u>edges(k) & u<=edges(k+1))*k; end
Y = zeros(n,1);
mu = truth.theta(1,:); sg = truth.theta(2,:);
for k = 1:K
    idx = (z==k);
    if any(idx)
        lam = lognrnd(mu(k), sg(k), sum(idx), 1);
        Y(idx) = poissrnd(lam);
    end
end
Y = double(Y(:));
end

function fit = align_pl_by_mean(fit)
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

% Build starts: warm (if provided) + fresh
starts = cell(nStarts,1); ns = 0;
if nargin>=7 && ~isempty(warm)
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
n = numel(Y); K = numel(piK); converged=false; prev_ll=-inf;

log_y_fact = gammaln(Y+1);
opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                    'MaxIterations',200,'MaxFunctionEvaluations',2e4);

for it = 1:MaxIter
    % E-step
    Pk = zeros(n,K);
    for k = 1:K
        mu = theta(1,k); sg = max(1e-8, theta(2,k));
        Pk(:,k) = pl_pmf_gh_Y(Y, log_y_fact, mu, sg, GH);
    end
    Pk = max(Pk, 1e-300);
    mix = max(Pk*piK(:), 1e-300);
    tau = bsxfun(@rdivide, bsxfun(@times, Pk, piK(:).'), mix);
    ll  = sum(log(mix));

    % M-step: weights
    Nk  = sum(tau,1)'; piK = max(eps, Nk/n); piK = piK/sum(piK);

    % M-step: (mu, sigma) via weighted GH-likelihood
    th = zeros(2,K);
    for k = 1:K
        wk = tau(:,k);
        mu0 = theta(1,k); sg0 = max(1e-6, theta(2,k));
        v0  = [mu0; log(sg0)];
        obj = @(v) neg_wPL_ll(Y, log_y_fact, wk, v, GH);
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
mu = v(1); sg = max(1e-8, exp(v(2)));
p  = pl_pmf_gh_Y(Y, log_y_fact, mu, sg, GH);
f  = -sum( w .* log( max(p,1e-300) ) );
end

% ----------- Robust initializations (kmeans / quantile / random) -------
function [pi0, theta0] = init_pl(Y, K, how)
X = log(Y + 0.5);               % stabilize zeros; approximate log(lambda)
switch lower(how)
    case 'kmeans'
        if exist('kmeans','file') == 2
            try
                % 1D k-means on X
                [lbl, ctr] = kmeans(X, K, 'Replicates', 3, 'MaxIter', 200, 'Display','off');
                pi0 = accumarray(lbl,1,[K,1]) / numel(X);
                mu0 = zeros(1,K); sg0 = zeros(1,K);
                for k=1:K
                    xk = X(lbl==k);
                    if isempty(xk), xk = X; end
                    mu0(k) = mean(xk);
                    sg0(k) = max(0.20, std(xk)); % lower bound for sigma
                end
                theta0 = [mu0; sg0];
                % order by larger PL mean first
                m = exp(mu0 + 0.5*sg0.^2);
                [~,ord] = sort(m,'descend');
                pi0 = pi0(ord); theta0 = theta0(:,ord);
                return
            catch
                % fall through to quantile
            end
        end
        % fall back to quantile if kmeans unavailable/fails
        how = 'quantile';
    case 'quantile'
        qx  = prctile(X, linspace(10,90,K));
        mu0 = qx(:).';
        sg0 = max(0.30, std(X))*ones(1,K);  % slightly larger default spread
        pi0 = ones(K,1)/K;
        theta0 = [mu0; sg0];
        % order by mean(lambda)
        m = exp(mu0 + 0.5*sg0.^2);
        [~,ord] = sort(m,'descend');
        pi0 = pi0(ord); theta0 = theta0(:,ord);
        return
    otherwise % 'random'
        a = rand(K,1)+0.1; pi0 = a/sum(a);
        mu0 = mean(X) + std(X)*randn(1,K);
        sg0 = max(0.25, 0.6*std(X))*ones(1,K);
        theta0 = [mu0; sg0];
        m = exp(mu0 + 0.5*sg0.^2);
        [~,ord] = sort(m,'descend');
        pi0 = pi0(ord); theta0 = theta0(:,ord);
        return
end
end

% ---------------- HMIX / VNED for PL (known K) ----------------
function fit = dic_pl_knownK(Y, K, method, nStarts, MaxIter, Tol, init, warm, GH)
[YY, g] = empirical(Y); g = max(g,1e-12);
best = struct('obj',Inf,'pi',[],'theta',[]);

% starts: warm + fresh
starts = cell(nStarts,1); ns=0;
if nargin>=8 && ~isempty(warm)
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
        R   = bsxfun(@rdivide, num, den);            % responsibilities on support
        R_safe = max(R, 1e-12); g_safe = max(g, 1e-12);

        theta_new = theta;

        switch lower(method)
            case 'hd'    % Hellinger step
                for k = 1:K
                    rk = R(:,k);
                    mu0 = theta(1,k); sg0 = max(1e-6, theta(2,k));
                    Hk = @(u) -sum( sqrt( max(1e-12,rk).*g_safe .* ...
                          pl_pmf_gh(YY, u(1), max(1e-8,exp(u(2))), GH) ) );
                    uk = fminunc(Hk, [mu0; log(sg0)], opts);
                    theta_new(:,k) = [uk(1); max(1e-8,exp(uk(2)))];
                end
            otherwise   % 'vned'
                for k = 1:K
                    rk = R(:,k); rk_safe = R_safe(:,k);
                    mu0 = theta(1,k); sg0 = max(1e-6, theta(2,k));
                    NEDk = @(u) sum( exp( - piK(k).*pl_pmf_gh(YY, u(1), max(1e-8,exp(u(2))), GH) ...
                                        ./ (g_safe.*rk_safe) ) .* g_safe .* rk );
                    uk  = fminunc(NEDk, [mu0; log(sg0)], opts);
                    theta_new(:,k) = [uk(1); max(1e-8,exp(uk(2)))];
                end
        end

        % Robust π update from responsibilities weighted by g
        Nk = sum(bsxfun(@times, R, g), 1)';        % Kx1
        pi_new = max(eps, Nk/sum(Nk));

        if max(abs([pi_new; theta_new(:)] - prev)) <= Tol
            piK = pi_new; theta = theta_new; break;
        end
        piK = pi_new; theta = theta_new; prev = [piK; theta(:)];
    end

    % divergence value for model selection
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

fit = struct('pi',best.pi,'theta',best.theta,'obj',best.obj);
end

% ---------- PL pmf via Gauss–Hermite quadrature ----------
function p = pl_pmf_gh_Y(Y, log_y_fact, mu, sigma, GH)
N = numel(Y);
p = zeros(N,1);
scale = GH.scale; s2 = GH.s2; x = GH.x; w = GH.w;
for i = 1:numel(x)
    lam = exp(mu + sigma*s2*x(i));
    lp = Y.*log(lam) - lam - log_y_fact;
    p  = p + w(i)*exp(lp);
end
p = scale * p;
p = max(p, 1e-300);
end

function p = pl_pmf_gh(YY, mu, sigma, GH)
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

% ---------- GH rule ----------
function GH = gh_rule(n)
i  = (1:n-1)';
a  = zeros(n,1);
b  = sqrt(i/2);
J  = diag(a) + diag(b,1) + diag(b,-1);
[V,D] = eig(J);
[x,idx] = sort(diag(D)); V = V(:,idx);
w = sqrt(pi) * (V(1,:)').^2;
GH.x = x; GH.w = w; GH.scale = 1/sqrt(pi); GH.s2 = sqrt(2);
end

% ---------- empirical ----------
function [YY, g] = empirical(Y)
[YY,~,ic] = unique(Y);
cnt = accumarray(ic,1);
g   = cnt / numel(Y);
end

% ======================================================================
% ============================ Averages ================================
% ======================================================================
function s = safe_first_avg(avg_params_field)
if isempty(avg_params_field)
    s = struct('K',NaN,'pi',[],'theta',[]); return;
end
if iscell(avg_params_field), s = avg_params_field{1}; else, s = avg_params_field; end
if ~isfield(s,'K'), s.K = numel(s.pi); end
if ~isfield(s,'pi'), s.pi = []; end
if ~isfield(s,'theta'), s.theta = []; end
end

function A = avg_from_cells(cells, model) %#ok<INUSD>
if isempty(cells) || all(cellfun(@isempty,cells))
    A = struct('mean_pi',NaN,'mean_theta',NaN); return;
end
Ks = nan(numel(cells),1);
for i=1:numel(cells)
    s = cells{i};
    if ~isempty(s) && isfield(s,'pi'), Ks(i) = numel(s.pi); end
end
Kmax = max(Ks, [], 'omitnan'); if isnan(Kmax), A = struct('mean_pi',NaN,'mean_theta',NaN); return; end
Pis = NaN(Kmax, numel(cells)); Ths = NaN(2*Kmax, numel(cells));
for i=1:numel(cells)
    s = cells{i};
    if isempty(s) || ~isfield(s,'pi') || ~isfield(s,'theta'), continue; end
    k = numel(s.pi); Pis(1:k,i) = s.pi(:);
    th = s.theta; if isvector(th), th = reshape(th(:),2,[]); end
    Ths(1:2*k,i) = th(:);
end
A.mean_pi = mean(Pis, 2, 'omitnan');
tmp = mean(Ths, 2, 'omitnan'); A.mean_theta = reshape(tmp, 2, []);
end

% ---------------- tiny helper ----------------
function val = tern(c,a,b)
if c, val=a; else, val=b; end
end
