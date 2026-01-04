function out = mix_split_select_estimate(Y, varargin)
% MIX_SPLIT_SELECT_ESTIMATE  (parallel-ready, HD legacy weights, safe optim)
% Unified split–select–estimate package for:
%   Methods: 'em'   (EM + BIC/AIC)
%            'hd'   (HMIX + DIC/HIC, H^2 = 2*(1 - sum sqrt(g f)))
%            'vned' (VNEDMIX + DIC/HIC, D = sum g*exp(-f/g))
%   Models : 'poiss' (Poisson), 'pg' (Poisson–Gamma / NB marginal),
%            'pl' (Poisson–Lognormal; fast GH quadrature)
%
% Key choices
% - HMIX/VNEDMIX use empirical kernel g_n on observed support.
% - PL uses fixed Gauss–Hermite quadrature (no integral()).
% - HD weights use *legacy* rule  w = fvals.^2  (restored).
% - All inner optimizations use *safe* wrappers; failures do not stop runs.
% - Hessian (for cov/SD) is computed ONCE per repetition and is optional
%   via 'ComputeHessian' (default false for speed).
%
% Parallelization
%   'UseParallel',true, 'ParallelAxis','R'  -> parfor over repetitions
%   'UseParallel',true, 'ParallelAxis','C'  -> parfor over splits (no nested)
%
% Examples
%   out = mix_split_select_estimate(Y,'Kmax',5,'C',6,'R',20, ...
%         'UseParallel',true,'ParallelAxis','R','ComputeHessian',false);
%   out = mix_split_select_estimate(Y,'models',{'pl'},'methods',{'em'}, ...
%         'R',50,'UseParallel',true,'ParallelAxis','C');

% ---------------- Parse inputs ----------------
ip = inputParser;
ip.addRequired('Y', @(x)isnumeric(x)&&isvector(x)&&all(x>=0)&&all(mod(x,1)==0));
ip.addParameter('models',  {'poiss','pg','pl'}, @(c)iscellstr(c)&&~isempty(c));
ip.addParameter('methods', {'em','hd','vned'},  @(c)iscellstr(c)&&~isempty(c));
ip.addParameter('Kmax',   5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('C',      5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('R',      1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('penalty','bic', @(s)ischar(s)&&ismember(lower(s),{'bic','aic'}));
ip.addParameter('dic_rule','stop', @(s)ischar(s)&&ismember(lower(s),{'stop','min'}));
ip.addParameter('nStarts',4,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('MaxIter',100,@(x)isnumeric(x)&&isscalar(x)&&x>=10);
ip.addParameter('Tol',1e-5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('init','quantile',@(s)ischar(s)&&ismember(lower(s),{'quantile','random'}));
ip.addParameter('seed',1,@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('trueK',[],@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('trueParams', [], @(s)isstruct(s) || isempty(s));
ip.addParameter('gh_n',20,@(x)isnumeric(x)&&isscalar(x)&&x>=5);    % PL GH nodes
ip.addParameter('ComputeHessian',false,@(b)islogical(b)||ismember(b,[0 1]));
% Parallel toggles:
ip.addParameter('UseParallel', false, @(b)islogical(b)||ismember(b,[0 1]));
ip.addParameter('ParallelAxis','R', @(s)ischar(s)&&ismember(lower(s),{'r','c','off'}));
ip.parse(Y,varargin{:});
S = ip.Results;
% 
rng(S.seed);
Y = double(Y(:));
n = numel(Y);

% penalty b(n) for DIC/HIC
switch lower(S.penalty)
    case 'aic', b_n = 1;
    otherwise,  b_n = log(n)/2; % 'bic'
end

% param count functions
pdim_of = @(model) (strcmpi(model,'poiss')*1 + strcmpi(model,'pg')*2 + strcmpi(model,'pl')*2);
vK_of   = @(model) @(K) (pdim_of(model)+1)*K - 1;

% Prepare GH rule once (for PL)
GH = get_gh_rule(S.gh_n);

% ---------------- Optional: spin up a pool ----------------
parAxis = lower(S.ParallelAxis);
if ~S.UseParallel || strcmp(parAxis,'off'), parAxis = 'off'; end
if S.UseParallel && ~strcmp(parAxis,'off')
    try
        if license('test','Distrib_Computing_Toolbox')
            p = gcp('nocreate'); if isempty(p), parpool('local'); end %#ok<PPOOL>
        else
            warning('Parallel toolbox not available; proceeding serially.');
            parAxis = 'off';
        end
    catch
        warning('Could not start parallel pool. Proceeding serially.');
        parAxis = 'off';
    end
end

% ---------------- Runner over (model, method) grid ----------------
models  = lower(S.models(:))';
methods = lower(S.methods(:))';

out = struct();
out.settings = S;
out.results  = struct();

for m = 1:numel(models)
    model = models{m};
    vK    = vK_of(model);
    pdim  = pdim_of(model);

    for t = 1:numel(methods)
        method = methods{t};

        % storage across repetitions
        Khat_mode   = zeros(S.R,1);
        Khat_counts = cell(S.R,1);
        avg_params  = cell(S.R,1);
        curves      = cell(S.R,1);
        cov_params  = cell(S.R,1);
        sd_params   = cell(S.R,1);
        MSE_cell    = cell(S.R,1);

        switch parAxis
            case 'r'
                % --------- PARALLEL across repetitions ---------
                parfor r = 1:S.R
                    rng(S.seed + 1000*r);
                    [Km, Kcnt, avgp, curveM, Ctheta, sdv, MSEr] = ...
                        run_one_rep(Y, n, S, model, method, b_n, vK, pdim, GH);
                    Khat_mode(r)   = Km;
                    Khat_counts{r} = Kcnt;
                    avg_params{r}  = avgp;
                    curves{r}      = curveM;
                    cov_params{r}  = Ctheta;
                    sd_params{r}   = sdv;
                    MSE_cell{r}    = MSEr;
                end
            otherwise
                % --------- SERIAL over repetitions; maybe parallel over splits ---------
                for r = 1:S.R
                    rng(S.seed + 1000*r);
                    [Km, Kcnt, avgp, curveM, Ctheta, sdv, MSEr] = ...
                        run_one_rep(Y, n, S, model, method, b_n, vK, pdim, GH, parAxis);
                    Khat_mode(r)   = Km;
                    Khat_counts{r} = Kcnt;
                    avg_params{r}  = avgp;
                    curves{r}      = curveM;
                    cov_params{r}  = Ctheta;
                    sd_params{r}   = sdv;
                    MSE_cell{r}    = MSEr;
                end
        end

        prop_correct = NaN;
        if ~isempty(S.trueK), prop_correct = mean(Khat_mode==S.trueK); end

        % stash results
        if ~isfield(out.results, model), out.results.(model) = struct(); end
        resS = struct( ...
            'Khat_mode',    Khat_mode, ...
            'prop_correct', prop_correct, ...
            'Khat_counts',  Khat_counts, ...
            'avg_params',   avg_params, ...
            'curves',       curves, ...
            'cov_params',   cov_params, ...
            'sd_params',    sd_params, ...
            'MSE',          MSE_cell, ...
            'settings',     rmfield(S, {'models','methods'}) );
        out.results.(model).(method) = resS;
    end
end
end

% =====================================================================
% =================== one repetition (optionally parallel splits) ======
% =====================================================================
function [K_mode, Kcnt_struct, avgp, curveMat, Ctheta, sdv, MSE_r] = ...
    run_one_rep(Y, n, S, model, method, b_n, vK, pdim, GH, parAxisSplits)

if nargin < 11, parAxisSplits = 'off'; end

Khats = zeros(S.C,1);
ests  = cell(S.C,1);
curveMat = nan(S.C, S.Kmax);

switch lower(parAxisSplits)
    case 'c'
        parfor c = 1:S.C
            rng(S.seed + 100000 + c);  % independent per split
            [Khats(c), curveMat(c,:), ests{c}] = run_one_split(Y, S, model, method, b_n, vK, pdim, GH);
        end
    otherwise
        for c = 1:S.C
            [Khats(c), curveMat(c,:), ests{c}] = run_one_split(Y, S, model, method, b_n, vK, pdim, GH);
        end
end

% modal K across splits & averaged parameters
K_unique = unique(Khats);
counts   = arrayfun(@(k) sum(Khats==k), K_unique);
[~,ix]   = max(counts);   K_mode = K_unique(ix);

keep = find(Khats==K_mode);
Pi = []; Th = [];
for j = 1:numel(keep)
    Pi = [Pi, ests{keep(j)}.pi(:)]; %#ok<AGROW>
    Th = [Th, ests{keep(j)}.theta(:)]; %#ok<AGROW>
end
pi_avg = mean(Pi,2); pi_avg = pi_avg/sum(pi_avg);
theta_avg = mean(Th,2);

Kcnt_struct = struct('K',K_unique,'count',counts);
avgp  = struct('K',K_mode,'pi',pi_avg,'theta',theta_avg);

% ---- Variance/SD via Hessian on full data (optional, once) ----
Ctheta = []; sdv = [];
if S.ComputeHessian
    [~, Jmap, zhat] = map_params_for_info(pi_avg, theta_avg, model);
    switch lower(method)
        case 'em'
            objH = @(z) negloglik_z(Y, z, K_mode, model, GH);
        case 'hd'
            [YYfull,gfull] = empirical(Y);
            objH = @(z) n * hd_objective_z(YYfull, gfull, z, K_mode, model, GH);
        case 'vned'
            [YYfull,gfull] = empirical(Y);
            objH = @(z) n * vned_objective_z(YYfull, gfull, z, K_mode, model, GH);
    end
    Hz = num_hessian(objH, zhat);
    Cz = safe_inv(Hz);
    Ctheta = Jmap * Cz * Jmap.';       % covariance for stacked [pi;theta]
    sdv = sqrt(max(0, diag(Ctheta)));
end

% ---- MSE (if truth provided and modal K equals trueK) ----
MSE_r = [];
if ~isempty(S.trueK) && isstruct(S.trueParams) && isfield(S.trueParams, model) ...
        && K_mode == S.trueK
    truth = S.trueParams.(model);
    [pi_true, th_true] = normalize_truth(truth, model, K_mode);
    pi_est = pi_avg(:);
    th_est = theta_avg(:);
    MSE_r = struct( ...
        'pi',    (pi_est - pi_true).^2, ...
        'theta', (th_est - th_true).^2, ...
        'overall', mean([ (pi_est - pi_true).^2 ; (th_est - th_true).^2 ]) );
end
end

% =====================================================================
% ========================= one split (helper) =========================
% =====================================================================
function [Khat, curveRow, fit2] = run_one_split(Y, S, model, method, b_n, vK, pdim, GH)
n = numel(Y);
idx = randperm(n); n1 = floor(n/2);
D1 = Y(idx(1:n1)); D2 = Y(idx(n1+1:end));

switch lower(method)
    case {'hd','vned'}  % ---- HMIX/VNED selection on D1 ----
        [Khat, ~, DICvals] = chooseK_by_dic(D1, model, method, S.Kmax, b_n, vK, S.dic_rule, S.init, S.nStarts, S.MaxIter, S.Tol, GH);
        curveRow = DICvals(:).';
        % refit on D2 with same method & Khat
        fit2 = fit_mix_dic(D2, Khat, model, method, S.init, S.nStarts, S.MaxIter, S.Tol, GH);
        fit2 = align_components(fit2, model);

    case 'em'           % ---- EM selection on D1 ----
        [Khat, critVals] = chooseK_by_em(D1, model, S.Kmax, pdim, S.penalty, S.nStarts, S.MaxIter, S.Tol, S.init, GH);
        curveRow = critVals(:).';
        % refit on D2 via EM
        fitE = em_mix_fit(D2, Khat, model, S.nStarts, S.MaxIter, S.Tol, S.init, GH);
        fitE = align_components_em(fitE, model);
        fit2 = struct('pi',fitE.pi,'theta',fitE.theta);
end
end

% =====================================================================
% =========================  SELECTION ROUTES  =========================
% =====================================================================
function [Khat, Dvals, DICvals] = chooseK_by_dic(D1, model, method, Kmax, b_n, vK, rule, init, nStarts, MaxIter, Tol, GH)
Kgrid   = 1:Kmax;
Dvals   = nan(Kmax,1);
DICvals = nan(Kmax,1);
[YY, g] = empirical(D1);
g_safe = max(g,1e-12);

for K = 1:Kmax
    best = struct('obj',Inf,'pi',[],'theta',[]);
    for s0 = 1:nStarts %#ok<NASGU>
        f1 = fit_one_dic(YY, g_safe, K, model, method, init, MaxIter, Tol, GH);
        f  = max(mix_pmf(YY, f1.pi, f1.theta, model, GH), 1e-15);
        switch lower(method)
            case 'hd'
                BC = sum( sqrt(g_safe .* f) );
                D  = 2*(1-BC);
            otherwise
                D  = sum( exp(-f./g_safe) .* g_safe );
        end
        if D < best.obj
            best.obj = D; best.pi = f1.pi; best.theta = f1.theta;
        end
    end
    Dvals(K)   = best.obj;
    DICvals(K) = Dvals(K) + (b_n/numel(D1))*vK(K);
end

switch lower(rule)
    case 'stop'
        Khat = Kgrid(end);
        for K = 1:Kmax-1
            etaK = (b_n/numel(D1))*(vK(K+1)-vK(K));
            if Dvals(K) <= Dvals(K+1)+etaK, Khat = K; break; end
        end
    otherwise
        [~,i] = min(DICvals); Khat = Kgrid(i);
end
end

function [Khat, critVals] = chooseK_by_em(D1, model, Kmax, pdim, penalty, nStarts, MaxIter, Tol, init, GH)
ll  = -inf(Kmax,1);
for K = 1:Kmax
    fitK = em_mix_fit(D1, K, model, nStarts, MaxIter, Tol, init, GH);
    ll(K) = fitK.loglik;
end
pK = @(K) pdim*K + (K-1);
switch lower(penalty)
    case 'bic'
        critVals = -2*ll + pK((1:Kmax).')*log(numel(D1));
    otherwise
        critVals =  2*pK((1:Kmax).') - 2*ll;
end
[~,Khat] = min(critVals);
end

% =====================================================================
% ============================  FITTERS  ===============================
% =====================================================================
function fit = fit_mix_dic(Y, K, model, method, init, nStarts, MaxIter, Tol, GH)
[YY, g] = empirical(Y);
g_safe = max(g,1e-12);
best = struct('obj',Inf,'pi',[],'theta',[]);

for s0 = 1:nStarts %#ok<NASGU>
    f1 = fit_one_dic(YY, g_safe, K, model, method, init, MaxIter, Tol, GH);
    f = max(mix_pmf(YY, f1.pi, f1.theta, model, GH), 1e-15);
    switch lower(method)
        case 'hd', D = 2*(1 - sum(sqrt(g_safe.*f)));
        otherwise, D = sum(exp(-f./g_safe).*g_safe);
    end
    if D < best.obj
        best.obj = D; best.pi = f1.pi; best.theta = f1.theta;
    end
end
fit = struct('pi',best.pi,'theta',best.theta);
end

function fit = fit_one_dic(YY, g, K, model, method, init, MaxIter, Tol, GH)
[piK, theta] = init_params(YY, K, model, init);
prev = [piK; theta(:)];
iters = 0;

while iters < MaxIter
    iters = iters + 1;

    % E-step on support
    Pk  = comp_pmf_matrix(YY, theta, model, GH);     % |support| x K
    num = bsxfun(@times, Pk, piK(:)');
    den = max(sum(num,2), 1e-15);
    R   = bsxfun(@rdivide, num, den);
    R_safe = max(R, 1e-12);
    g_safe = max(g, 1e-12);

    theta_new = theta; Kc = K;

    switch lower(method)
        case 'hd'      % ----- HMIX (Hellinger) with legacy weights -----
            fvals = zeros(Kc,1);
            for k = 1:Kc
                rk = R(:,k);
                switch lower(model)
                    case 'poiss'
                        lam0 = max(1e-8, theta(k));
                        Hk = @(t) -sum( sqrt( max(1e-12,rk).*g_safe .* poisspdf(YY, max(1e-12,exp(t))) ) );
                        [tk,fk,ok] = fminunc_safe(Hk, log(lam0));
                        if ~ok, tk = log(lam0); fk = Hk(tk); end
                        theta_new(k) = max(1e-8, exp(tk)); fvals(k) = fk;

                    case 'pg'
                        a0 = max(1e-8, theta(1,k)); b0 = max(1e-8, theta(2,k));
                        Hk = @(tt) -sum( sqrt( max(1e-12,rk).*g_safe .* ...
                                    nbinpdf(YY, max(1e-12,exp(tt(1))), ...
                                         max(1e-12,exp(tt(2)))/(1+max(1e-12,exp(tt(2))))) ) );
                        [ttk,fk,ok] = fminunc_safe(Hk, [log(a0); log(b0)]);
                        if ~ok, ttk = [log(a0); log(b0)]; fk = Hk(ttk); end
                        theta_new(:,k) = [max(1e-8,exp(ttk(1))); max(1e-8,exp(ttk(2)))]; fvals(k) = fk;

                    case 'pl'
                        mu0 = theta(1,k); sig0 = max(1e-6, theta(2,k));
                        Hk = @(uu) -sum( sqrt( max(1e-12,rk).*g_safe .* ...
                                         pl_pmf_gh(YY, uu(1), max(1e-8,exp(uu(2))), GH) ) );
                        [uopt,fk,ok] = fminunc_safe(Hk, [mu0; log(sig0)]);
                        if ~ok, uopt = [mu0; log(sig0)]; fk = Hk(uopt); end
                        theta_new(:,k) = [uopt(1); max(1e-8,exp(uopt(2)))]; fvals(k) = fk;
                end
            end
            % LEGACY weight rule (restored) with safe fallback
            w = fvals.^2;
            if ~all(isfinite(w)) || sum(w)<=0, w = ones(Kc,1); end
            pi_new = max(eps, w/sum(w));

        otherwise      % ----- VNEDMIX (VNED) -----
            vned_k = zeros(Kc,1);
            for k = 1:Kc
                rk = R(:,k); rk_safe = R_safe(:,k);
                switch lower(model)
                    case 'poiss'
                        lam0 = max(1e-8, theta(k));
                        NEDk = @(t) sum( exp( - piK(k)*poisspdf(YY, max(1e-12,exp(t)))./(g_safe.*rk_safe) ) .* g_safe .* rk );
                        [tk,~,ok] = fminunc_safe(NEDk, log(lam0));
                        if ~ok, tk = log(lam0); end
                        theta_new(k) = max(1e-8, exp(tk));
                        pmf  = poisspdf(YY, theta_new(k));
                        v = sum( exp( -piK(k)*pmf./(g_safe.*rk_safe) ) .* (piK(k)*pmf) );
                        vned_k(k) = max(1e-12, v);

                    case 'pg'
                        a0 = max(1e-8, theta(1,k)); b0 = max(1e-8, theta(2,k));
                        NEDk = @(tt) sum( exp( - piK(k).*nbinpdf(YY, max(1e-12,exp(tt(1))), ...
                                   max(1e-12,exp(tt(2)))/(1+max(1e-12,exp(tt(2))))) ./ (g_safe.*rk_safe) ) .* g_safe .* rk );
                        [ttk,~,ok]  = fminunc_safe(NEDk, [log(a0); log(b0)]);
                        if ~ok, ttk = [log(a0); log(b0)]; end
                        theta_new(:,k) = [max(1e-8,exp(ttk(1))); max(1e-8,exp(ttk(2)))];
                        pmf = nbinpdf(YY, theta_new(1,k), theta_new(2,k)/(1+theta_new(2,k)));
                        v = sum( exp( - piK(k).*pmf./(g_safe.*rk_safe) ) .* (piK(k).*pmf) );
                        vned_k(k) = max(1e-12, v);

                    case 'pl'
                        mu0 = theta(1,k); sig0 = max(1e-6, theta(2,k));
                        NEDk = @(uu) sum( exp( - piK(k).*pl_pmf_gh(YY, uu(1), max(1e-8,exp(uu(2))), GH)./(g_safe.*rk_safe) ) .* g_safe .* rk );
                        [uopt,~,ok] = fminunc_safe(NEDk, [mu0; log(sig0)]);
                        if ~ok, uopt = [mu0; log(sig0)]; end
                        theta_new(:,k) = [uopt(1); max(1e-8,exp(uopt(2)))];
                        pmf = pl_pmf_gh(YY, theta_new(1,k), theta_new(2,k), GH);
                        v = sum( exp( - piK(k).*pmf./(g_safe.*rk_safe) ) .* (piK(k).*pmf) );
                        vned_k(k) = max(1e-12, v);
                end
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

function fit = em_mix_fit(Y, K, model, nStarts, MaxIter, Tol, init, GH)
% EM MLE fitter (best over nStarts) with numerical M-steps for PG/PL
Y = double(Y(:));

best = struct('ll',-inf,'pi',[],'theta',[],'iters',0,'converged',false);
for s0 = 1:nStarts %#ok<NASGU>
    [piK, theta] = init_params(Y, K, model, init);
    [ll, piK, theta, iters, conv] = em_core(Y, piK, theta, model, MaxIter, Tol, GH);
    if ll > best.ll
        best.ll = ll; best.pi = piK; best.theta = theta; best.iters = iters; best.converged = conv;
    end
end
fit = struct('pi',best.pi,'theta',best.theta,'loglik',best.ll,'iters',best.iters,'converged',best.converged);
end

function [ll, piK, theta, iters, converged] = em_core(Y, piK, theta, model, MaxIter, Tol, GH)
n = numel(Y); K = numel(piK); converged = false; prev_ll = -inf;
for it = 1:MaxIter
    % E-step
    Pk = comp_pmf_matrix_Y(Y, theta, model, GH);  % n x K
    mix = max(Pk*piK(:), 1e-300);
    tau = bsxfun(@rdivide, bsxfun(@times, Pk, piK(:).'), mix);
    ll  = sum(log(mix));

    % M-step: weights
    Nk  = sum(tau,1)'; piK = max(eps, Nk/n);

    % M-step: parameters
    switch lower(model)
        case 'poiss'
            theta = max(1e-8, (tau.'*Y) ./ max(eps,Nk));

        case 'pg'   % weighted NB MLE via safe fminunc
            th = zeros(2,K);
            for k = 1:K
                a0 = max(1e-3, theta(1,k)); b0 = max(1e-3, theta(2,k));
                wk = tau(:,k);
                obj = @(u) -sum( wk .* log( nbinpdf(Y, max(1e-8,exp(u(1))), ...
                                  max(1e-8,exp(u(2)))/(1+max(1e-8,exp(u(2))))) + 1e-300 ) );
                [uopt,~,ok] = fminunc_safe(obj, [log(a0); log(b0)]);
                if ~ok, uopt = [log(a0); log(b0)]; end
                th(:,k) = [max(1e-8,exp(uopt(1))); max(1e-8,exp(uopt(2)))];
            end
            theta = th;

        case 'pl'   % weighted PL MLE via GH pmf (safe)
            th = zeros(2,K);
            for k = 1:K
                mu0 = theta(1,k); sig0 = max(1e-6, theta(2,k));
                wk  = tau(:,k);
                obj = @(v) -sum( wk .* log( pl_pmf_gh(Y, v(1), max(1e-8,exp(v(2))), GH) + 1e-300 ) );
                [vopt,~,ok] = fminunc_safe(obj, [mu0; log(sig0)]);
                if ~ok, vopt = [mu0; log(sig0)]; end
                th(:,k) = [vopt(1); max(1e-8,exp(vopt(2)))];
            end
            theta = th;
    end

    if it>1 && abs(ll - prev_ll) < Tol*max(1,abs(prev_ll)), converged = true; break; end
    prev_ll = ll;
end
iters = it;
end

% =====================================================================
% ============================  UTILITIES  =============================
% =====================================================================
function [YY, g] = empirical(Y)
tbl = tabulate(Y);
YY = tbl(:,1);
g  = tbl(:,3)/100;
end

function P = comp_pmf_matrix(YY, theta, model, GH)
switch lower(model)
    case 'poiss'
        K = numel(theta);
        P = poisspdf(repmat(YY,1,K), repmat(theta(:).', numel(YY), 1));
    case 'pg'
        K = size(theta,2);
        P = zeros(numel(YY), K);
        for k = 1:K
            a = theta(1,k); b = theta(2,k);
            P(:,k) = nbinpdf(YY, a, b/(1+b));
        end
    case 'pl'
        K = size(theta,2);
        P = zeros(numel(YY), K);
        for k = 1:K
            P(:,k) = pl_pmf_gh(YY, theta(1,k), max(1e-8,theta(2,k)), GH);
        end
end
P = max(P, 1e-15);
end

function f = mix_pmf(YY, piK, theta, model, GH)
P = comp_pmf_matrix(YY, theta, model, GH);
f = max(P * piK(:), 1e-15);
end

function [pi0, theta0] = init_params(Y, K, model, how)
switch lower(how)
    case 'random', a = rand(K,1)+0.1; pi0 = a/sum(a);
    otherwise,     pi0 = ones(K,1)/K; % 'quantile'
end
switch lower(model)
    case 'poiss'
        qs = linspace(5,95,K); lam0 = max(1e-3, prctile(Y, qs).');
        if any(~isfinite(lam0)), mY = max(1e-3, mean(Y)); lam0 = mY*linspace(0.5,1.5,K).'; end
        theta0 = lam0(:);
    case 'pg'
        qs = linspace(5,95,K); cuts = prctile(Y, qs);
        if any(~isfinite(cuts)), cuts = linspace(min(Y), max(Y), K); end
        mu_guess = max(1e-2, cuts(:));
        vY = max(1, var(Y));
        b_guess = max(0.1, mu_guess./max(mu_guess.^2, vY));
        a_guess = max(0.2, mu_guess .* b_guess);
        theta0 = [a_guess.'; b_guess.'];
    case 'pl'
        qs = linspace(5,95,K); qy = prctile(Y+0.5, qs);
        mu0 = log(max(0.1, qy(:))); sig0 = 0.5*ones(K,1);
        theta0 = [mu0.'; sig0.'];
end
end

function fit = align_components(fit, model)
switch lower(model)
    case 'poiss'
        mu = fit.theta(:); [~,ord] = sort(mu,'ascend');
        fit.theta = fit.theta(ord); fit.pi = fit.pi(ord);
    case 'pg'
        a = fit.theta(1,:); b = fit.theta(2,:); mu = a./b;
        [~,ord] = sort(mu,'ascend'); fit.theta = fit.theta(:,ord); fit.pi = fit.pi(ord);
    case 'pl'
        mu = fit.theta(1,:); sig = fit.theta(2,:); m = exp(mu + 0.5*sig.^2);
        [~,ord] = sort(m,'ascend'); fit.theta = fit.theta(:,ord); fit.pi = fit.pi(ord);
end
end

function fit = align_components_em(fit, model)
switch lower(model)
    case 'poiss'
        mu = fit.theta(:); [~,ord] = sort(mu,'ascend');
        fit.theta = fit.theta(ord); fit.pi = fit.pi(ord);
    case 'pg'
        a = fit.theta(1,:); b = fit.theta(2,:); mu = a./b;
        [~,ord] = sort(mu,'ascend'); fit.theta = fit.theta(:,ord); fit.pi = fit.pi(ord);
    case 'pl'
        mu = fit.theta(1,:); sig = fit.theta(2,:); m = exp(mu + 0.5*sig.^2);
        [~,ord] = sort(m,'ascend'); fit.theta = fit.theta(:,ord); fit.pi = fit.pi(ord);
end
end

function P = comp_pmf_matrix_Y(Y, theta, model, GH)
n = numel(Y);
switch lower(model)
    case 'poiss'
        K = numel(theta);
        P = poisspdf(repmat(Y,1,K), repmat(theta(:).', n, 1));
    case 'pg'
        K = size(theta,2); P = zeros(n,K);
        for k = 1:K
            a = theta(1,k); b = theta(2,k);
            P(:,k) = nbinpdf(Y, a, b/(1+b));
        end
    case 'pl'
        K = size(theta,2); P = zeros(n,K);
        for k = 1:K
            P(:,k) = pl_pmf_gh(Y, theta(1,k), max(1e-8,theta(2,k)), GH);
        end
end
P = max(P, 1e-300);
end

% =====================================================================
% =====================  INFORMATION / MSE HELPERS  ===================
% =====================================================================
function [pvec, J, z] = map_params_for_info(pi, theta, model)
pi = pi(:);
K  = numel(pi);

switch lower(model)
    case 'poiss'
        lam = theta(:); lam = max(1e-12, lam(:));
        raw = log(lam);                 Dth = diag(lam);
    case 'pg'
        th = theta; if isvector(th), th = reshape(th(:),2,[]); end
        a = max(1e-12, th(1,:)).';
        b = max(1e-12, th(2,:)).';
        raw = [log(a); log(b)];        Dth = blkdiag(diag(a), diag(b));
    case 'pl'
        th = theta; if isvector(th), th = reshape(th(:),2,[]); end
        mu  = th(1,:).'; sig = max(1e-12, th(2,:).');
        raw = [mu; log(sig)];          Dth = blkdiag(eye(K), diag(sig));
end

alpha = log(max(1e-12, pi(1:K-1) ./ max(1e-12, pi(K))));
z = [alpha; raw];

PiJ = diag(pi) - pi*pi.';    % softmax Jacobian wrt [alpha;0]
Jpi = PiJ(:,1:K-1);

J = blkdiag(Jpi, Dth);
pvec = [pi; theta(:)];
end

function H = num_hessian(fun, z0)
% Central-difference Hessian
f0 = fun(z0);
q  = numel(z0);
H  = zeros(q,q);
eps0 = 1e-4;
for i = 1:q
    ei = zeros(q,1); hi = eps0*(1+abs(z0(i))); ei(i)=hi;
    f_ip = fun(z0+ei); f_im = fun(z0-ei);
    H(i,i) = (f_ip - 2*f0 + f_im) / (hi^2);
    for j = i+1:q
        ej = zeros(q,1); hj = eps0*(1+abs(z0(j))); ej(j)=hj;
        fpp = fun(z0+ei+ej);
        fpm = fun(z0+ei-ej);
        fmp = fun(z0-ei+ej);
        fmm = fun(z0-ei-ej);
        H(i,j) = (fpp - fpm - fmp + fmm) / (4*hi*hj);
        H(j,i) = H(i,j);
    end
end
H = (H+H.')/2 + 1e-9*eye(q);
end

function C = safe_inv(H)
[U,S,~] = svd((H+H.')/2);
s = diag(S);
thr = max(1e-12, 1e-12*max(s));
sInv = 1./max(s,thr);
C = U*diag(sInv)*U.';
end

function val = negloglik_z(Y, z, K, model, GH)
[pi, theta] = unpack_z(z, K, model);
P = comp_pmf_matrix_Y(Y, theta, model, GH);
mix = max(P*pi(:), 1e-300);
val = -sum(log(mix));
end

function val = hd_objective_z(YY, g, z, K, model, GH)
[pi, theta] = unpack_z(z, K, model);
f = max(mix_pmf(YY, pi, theta, model, GH), 1e-15);
BC = sum(sqrt(max(1e-12,g).*f));
val = 2*(1-BC);   % H^2
end

function val = vned_objective_z(YY, g, z, K, model, GH)
[pi, theta] = unpack_z(z, K, model);
f = max(mix_pmf(YY, pi, theta, model, GH), 1e-15);
g_safe = max(g,1e-12);
val = sum( exp(-f./g_safe) .* g_safe );
end

function [pi, theta] = unpack_z(z, K, model)
alpha = z(1:K-1);
raw   = z(K:end);
a = [alpha; 0];
e = exp(a - max(a));
pi = e / sum(e);
switch lower(model)
    case 'poiss'
        theta = max(1e-12, exp(raw(:))).';
    case 'pg'
        a = max(1e-12, exp(raw(1:K)));
        b = max(1e-12, exp(raw(K+1:2*K)));
        theta = [a.'; b.'];
    case 'pl'
        mu  = raw(1:K).';
        sig = max(1e-12, exp(raw(K+1:2*K))).';
        theta = [mu; sig];
end
end

function [pi_true, th_true] = normalize_truth(truth, model, K)
pi_true = truth.pi(:);
if numel(pi_true) ~= K, error('trueParams.%s.pi length mismatch.', model); end
switch lower(model)
    case 'poiss'
        th = truth.theta(:);
        if numel(th) ~= K, error('trueParams.%s.theta length mismatch.', model); end
        th_true = th(:);
    otherwise
        th = truth.theta;
        if isvector(th) && numel(th)==2*K, th = reshape(th,2,K); end
        if ~ismatrix(th) || size(th,1)~=2 || size(th,2)~=K
            error('trueParams.%s.theta must be 2xK or length 2K.', model);
        end
        th_true = th(:);
end
end

% =====================================================================
% ===================  GAUSS–HERMITE & PL PMF (FAST)  =================
% =====================================================================
function GH = get_gh_rule(n)
% Golub–Welsch GH (physicists'), weight exp(-x^2)
% ∫ e^{-x^2} f(x) dx ≈ sum w_i f(x_i)
% PL: p(y) = ∫ φ(z) Pois(y|e^{μ+σ z}) dz, with z=√2 x:
% p(y) ≈ (1/√π) * sum w_i * Pois(y | e^{μ + σ√2 x_i})
beta  = sqrt((1:(n-1))/2);
J     = diag(zeros(n,1)) + diag(beta,1) + diag(beta,-1);
[V,D] = eig(J);
x = diag(D);
w = (V(1,:).^2)'*sqrt(pi);
GH.x = x;
GH.w = w;
GH.scale = 1/sqrt(pi);
end

function p = pl_pmf_gh(YY, mu, sigma, GH)
% Vectorized PL pmf on integer support YY using GH nodes
% p(y) ≈ (1/√π) ∑ w_i Pois(y | exp(mu + sigma*√2*x_i))
x = GH.x; w = GH.w; sc = GH.scale;
lam = exp(mu + sigma*sqrt(2)*x(:)');     % 1 x nGH
P = poisspdf(YY(:), lam);                % |YY| x nGH
p = sc * (P * w);                        % |YY| x 1
p = max(p, 1e-15);
end

% =====================================================================
% ========================== SAFE OPT WRAPPER ==========================
% =====================================================================
function [x,fval,ok] = fminunc_safe(fun, x0)
% Defensive wrapper around fminunc that:
%  - checks fun(x0); if bad, tries small jitters; if still bad, returns x0
%  - catches optimizer failures and returns the best-so-far
ok = true; fval = NaN; x = x0(:);  % ensure column
opts = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
    'MaxIterations',300,'MaxFunctionEvaluations',3e4);

% ensure finite start
[f0,good0] = eval_fun(fun, x);
if ~good0
    jitters = [0, 0.1, -0.1, 0.3, -0.3];
    for j = 1:numel(jitters)
        x_try = x + jitters(j)*ones(size(x));
        [fj, goodj] = eval_fun(fun, x_try);
        if goodj, x = x_try; f0 = fj; break; end
    end
    if ~isfinite(f0)
        ok=false; fval = Inf; return;
    end
end

try
    [x1,f1,flag,~] = fminunc(fun, x, opts); %#ok<ASGLU>
    if ~isfinite(f1), x = x; fval = f0; ok=false; else, x = x1; fval = f1; end
catch
    ok = false; fval = f0; % fall back to start
end
end

function [f,good] = eval_fun(fun, x)
try
    f = fun(x);
    good = isfinite(f);
catch
    f = Inf; good = false;
end
end
