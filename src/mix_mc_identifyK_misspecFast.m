function out = mix_mc_identifyK_misspecFast(B, n, fitModel, genParams, varargin)
% MIX_MC_IDENTIFYK_MISSPECFast
% Monte Carlo over B datasets under MODEL MISSPECIFICATION:
%   - Generate data from genModel/genParams
%   - Fit mixtures using fitModel with EM / HD / VNED
%
% Per method (EM/HD/VNED):
%   1) Select K̂ on a SUBSAMPLE (size n_sub) with lean settings (sel_C splits).
%   2) Refit ONCE on FULL data at that method’s K̂ with fit_C splits.
% Parallelized over datasets (no nested parfor). Robust to optimizer errors.
%
% Required dependency: mix_split_select_estimate.m (your fast version).
%
% Signature:
%   out = mix_mc_identifyK_misspecFast(B, n, fitModel, genParams, 'genModel','poiss', ...)
%
% Example (generate Poisson, fit PG):
%   gen.poiss.pi = [0.3; 0.7];
%   gen.poiss.theta = [10; 20];                     % Poisson rates
%   out = mix_mc_identifyK_misspecFast( ...
%           500, 5000, 'pg', gen.poiss, ...
%           'genModel','poiss', 'methods',{'em','hd','vned'}, ...
%           'Kmax',5, 'sel_C',1, 'fit_C',5, 'sel_nStarts',1, 'fit_nStarts',2, ...
%           'sel_MaxIter',60, 'fit_MaxIter',120, 'gh_n',12, 'seed',11);

% ---------------- Parse inputs ----------------
ip = inputParser;
ip.addRequired('B', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('n', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('fitModel', @(s)ischar(s)&&ismember(lower(s),{'poiss','pg','pl'}));
ip.addRequired('genParams', @(s)isstruct(s) && isfield(s,'pi') && isfield(s,'theta'));

% What we FIT (default PG for misspec example)
ip.addParameter('methods', {'em','hd','vned'}, @(c)iscellstr(c)&&~isempty(c));
ip.addParameter('Kmax',   5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% What we GENERATE (default Poisson)
ip.addParameter('genModel','poiss', @(s)ischar(s)&&ismember(lower(s),{'poiss','pg','pl'}));

% Selection on subsample (fast)
ip.addParameter('n_sub', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>=100);
ip.addParameter('sel_C', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('sel_nStarts', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('sel_MaxIter', 60, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
ip.addParameter('sel_Tol', 1e-5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('sel_rule','stop', @(s)ischar(s)&&ismember(lower(s),{'stop','min'}));
ip.addParameter('penalty','bic', @(s)ischar(s)&&ismember(lower(s),{'bic','aic'}));
ip.addParameter('init','quantile',@(s)ischar(s)&&ismember(lower(s),{'quantile','random'}));

% Full-data refits at fixed K (per method)
ip.addParameter('fit_C', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('fit_nStarts', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('fit_MaxIter', 120, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
ip.addParameter('fit_Tol', 1e-5, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% PL Gauss–Hermite nodes (smaller for speed)
ip.addParameter('gh_n', 12, @(x)isnumeric(x)&&isscalar(x)&&x>=8);

% RNG + parallel
ip.addParameter('seed', 1, @(x)isnumeric(x)&&isscalar(x));
ip.addParameter('UseParallel', true, @(b)islogical(b)||ismember(b,[0 1]));
ip.parse(B,n,fitModel,genParams,varargin{:});
S = ip.Results;

methods  = lower(S.methods(:))';
doEM     = any(strcmp(methods,'em'));
doHD     = any(strcmp(methods,'hd'));
doVNED   = any(strcmp(methods,'vned'));

fitModel = lower(S.fitModel);
genModel = lower(S.genModel);

Ktrue = numel(S.genParams.pi(:));   % true K from the *generator*

% ---------------- Outputs (parfor-sliced) ----------------
khat_EM = NaN(B,1);   khat_HD = NaN(B,1);   khat_VN = NaN(B,1);
est_EM  = cell(B,1);  est_HD  = cell(B,1);  est_VN  = cell(B,1);
sel_counts_EM = cell(B,1); sel_counts_HD = cell(B,1); sel_counts_VN = cell(B,1);
sel_frac_EM = NaN(B,1); sel_frac_HD = NaN(B,1); sel_frac_VN = NaN(B,1);

% ---------------- Pool ----------------
if S.UseParallel
    try
        if license('test','Distrib_Computing_Toolbox')
            p = gcp('nocreate'); if isempty(p), parpool('local'); end %#ok<PPOOL>
        else
            warning('Parallel toolbox not available; proceeding serially.');
            S.UseParallel = false;
        end
    catch
        S.UseParallel = false;
    end
end

% ---------------- MC over B datasets ----------------
if S.UseParallel
    parfor b = 1:B
        rng(S.seed + b);
        [khat_EM(b), est_EM{b}, sel_counts_EM{b}, sel_frac_EM(b), ...
         khat_HD(b), est_HD{b}, sel_counts_HD{b}, sel_frac_HD(b), ...
         khat_VN(b), est_VN{b}, sel_counts_VN{b}, sel_frac_VN(b)] = ...
            run_one_dataset(b, S, fitModel, genModel, S.genParams, methods, Ktrue);
    end
else
    for b = 1:B
        rng(S.seed + b);
        [khat_EM(b), est_EM{b}, sel_counts_EM{b}, sel_frac_EM(b), ...
         khat_HD(b), est_HD{b}, sel_counts_HD{b}, sel_frac_HD(b), ...
         khat_VN(b), est_VN{b}, sel_counts_VN{b}, sel_frac_VN(b)] = ...
            run_one_dataset(b, S, fitModel, genModel, S.genParams, methods, Ktrue);
    end
end

% ---------------- % correct and averages (per method) ----------------
prop_correct = struct();
avg_overall  = struct();
avg_when_correct = struct();

if doEM
    prop_correct.EM = 100*mean(khat_EM==Ktrue,'omitnan');
    avg_overall.EM  = avg_from_cells(est_EM,  fitModel);
    avg_when_correct.EM = avg_from_cells(est_EM(khat_EM==Ktrue), fitModel);
end
if doHD
    prop_correct.HD = 100*mean(khat_HD==Ktrue,'omitnan');
    avg_overall.HD  = avg_from_cells(est_HD,  fitModel);
    avg_when_correct.HD = avg_from_cells(est_HD(khat_HD==Ktrue), fitModel);
end
if doVNED
    prop_correct.VNED = 100*mean(khat_VN==Ktrue,'omitnan');
    avg_overall.VNED  = avg_from_cells(est_VN, fitModel);
    avg_when_correct.VNED = avg_from_cells(est_VN(khat_VN==Ktrue), fitModel);
end

% ---------------- Pack ----------------
out = struct();
out.khat = struct();
if doEM,   out.khat.EM   = khat_EM; end
if doHD,   out.khat.HD   = khat_HD; end
if doVNED, out.khat.VNED = khat_VN; end

out.sel_counts = struct();
if doEM,   out.sel_counts.EM   = sel_counts_EM; end
if doHD,   out.sel_counts.HD   = sel_counts_HD; end
if doVNED, out.sel_counts.VNED = sel_counts_VN; end

out.sel_frac_true = struct();
if doEM,   out.sel_frac_true.EM   = sel_frac_EM; end
if doHD,   out.sel_frac_true.HD   = sel_frac_HD; end
if doVNED, out.sel_frac_true.VNED = sel_frac_VN; end

out.est = struct();
if doEM,   out.est.EM   = est_EM;   end
if doHD,   out.est.HD   = est_HD;   end
if doVNED, out.est.VNED = est_VN;   end

out.avg_overall       = avg_overall;
out.avg_when_correct  = avg_when_correct;
out.prop_correct      = prop_correct;

out.settings = S;
out.settings.trueK     = Ktrue;          % generator's K
out.settings.fitModel  = fitModel;
out.settings.genModel  = genModel;
out.settings.genParams = S.genParams;
end

% ======================================================================
% ======================= One dataset (per method) =====================
% ======================================================================
function [K_em, estEM, C_em, F_em, ...
          K_hd, estHD, C_hd, F_hd, ...
          K_vn, estVN, C_vn, F_vn] = ...
          run_one_dataset(b, S, fitModel, genModel, genParams, methods, Ktrue)

% 0) simulate full data from the *generator*
Yfull = generate_dataset(S.n, genModel, genParams);

% build subsample for selection (fast)
n_sub = min(S.n, S.n_sub);
if n_sub < S.n
    idx = randperm(S.n, n_sub);
    Ysel = Yfull(idx);
else
    Ysel = Yfull;
end

% Defaults
K_em = NaN; estEM=[]; C_em=[]; F_em=NaN;
K_hd = NaN; estHD=[]; C_hd=[]; F_hd=NaN;
K_vn = NaN; estVN=[]; C_vn=[]; F_vn=NaN;

% ---- EM selection on subsample (cheap) ----
if any(strcmp(methods,'em'))
    [K_em, C_em, F_em] = try_select(Ysel, fitModel, 'em', S, b, Ktrue);
    estEM = try_refit(Yfull, fitModel, 'em', S, b, K_em);
end

% ---- HD selection on subsample (cheap) ----
if any(strcmp(methods,'hd'))
    [K_hd, C_hd, F_hd] = try_select(Ysel, fitModel, 'hd', S, b, Ktrue);
    estHD = try_refit(Yfull, fitModel, 'hd', S, b, K_hd);
end

% ---- VNED selection on subsample (cheap) ----
if any(strcmp(methods,'vned'))
    [K_vn, C_vn, F_vn] = try_select(Ysel, fitModel, 'vned', S, b, Ktrue);
    estVN = try_refit(Yfull, fitModel, 'vned', S, b, K_vn);
end
end

% ---------- selection helper (subsample; C=sel_C; no Hessian) ----------
function [Khat, Kcnt, fracTrue] = try_select(Ysel, fitModel, method, S, b, Ktrue)
Khat = NaN; Kcnt = []; fracTrue = NaN;
try
    out_sel = mix_split_select_estimate( ...
        Ysel, ...
        'models', {fitModel}, 'methods', {method}, ...
        'Kmax', S.Kmax, 'C', S.sel_C, 'R', 1, ...
        'penalty', S.penalty, 'dic_rule', S.sel_rule, ...
        'nStarts', S.sel_nStarts, 'MaxIter', S.sel_MaxIter, 'Tol', S.sel_Tol, ...
        'init', S.init, 'seed', S.seed + 101*b + method_seed_offset(method), ...
        'trueK', Ktrue, ...             % reporting only
        'trueParams', [], ...           % <-- misspec: leave empty to skip MSE
        'gh_n', S.gh_n, ...
        'ComputeHessian', false, ...
        'UseParallel', false);
    Rm = out_sel.results.(fitModel).(method);
    Khat = Rm.Khat_mode(1);
    Kcnt = Rm.Khat_counts{1};
    if ~isempty(Kcnt) && isfield(Kcnt,'K') && isfield(Kcnt,'count')
        pos = find(Kcnt.K==Ktrue,1);
        fracTrue = ~isempty(pos) * (Kcnt.count(pos) / S.sel_C);
    end
catch
    % leave NaNs/empty
end
end

% ---------- refit helper (full data; C=fit_C; no Hessian) ----------
function est = try_refit(Yfull, fitModel, method, S, b, Khat)
est = [];
if isnan(Khat) || Khat<1, return; end
try
    out_fit = mix_split_select_estimate( ...
        Yfull, ...
        'models', {fitModel}, 'methods', {method}, ...
        'Kmax', Khat, 'C', S.fit_C, 'R', 1, ...
        'penalty', S.penalty, 'dic_rule', 'stop', ...
        'nStarts', S.fit_nStarts, 'MaxIter', S.fit_MaxIter, 'Tol', S.fit_Tol, ...
        'init', S.init, 'seed', S.seed + 313*b + method_seed_offset(method), ...
        'trueK', Khat, 'trueParams', [], ...
        'gh_n', S.gh_n, ...
        'ComputeHessian', false, ...
        'UseParallel', false);
    est = safe_first_avg(out_fit.results.(fitModel).(method).avg_params);
catch
    % leave empty
end
end

function off = method_seed_offset(m)
switch lower(m)
    case 'em',   off = 11;
    case 'hd',   off = 17;
    otherwise,   off = 23; % 'vned'
end
end

% ======================================================================
% =============================== Utils ================================
% ======================================================================
function Y = generate_dataset(n, genModel, truth)
piK = truth.pi(:); K = numel(piK);
% draw component labels
u = rand(n,1); edges = [0; cumsum(piK)];
z = zeros(n,1);
for k = 1:K, z = z + (u>edges(k) & u<=edges(k+1))*k; end
Y = zeros(n,1);
switch lower(genModel)
    case 'poiss'
        lam = truth.theta(:);
        for k = 1:K
            idx = (z==k);
            if any(idx), Y(idx) = poissrnd(lam(k), sum(idx), 1); end
        end
    case 'pg'
        a = truth.theta(1,:); b = truth.theta(2,:); % NB: p = b/(1+b)
        for k = 1:K
            idx = (z==k);
            if any(idx)
                p = b(k)/(1+b(k));
                Y(idx) = nbinrnd(a(k), p, sum(idx), 1);
            end
        end
    case 'pl'
        mu = truth.theta(1,:); sig = truth.theta(2,:); % lognormal rate
        for k = 1:K
            idx = (z==k);
            if any(idx)
                lam = lognrnd(mu(k), sig(k), sum(idx), 1);
                Y(idx) = poissrnd(lam);
            end
        end
end
Y = double(Y(:));
end

function s = safe_first_avg(avg_params_field)
if iscell(avg_params_field)
    if isempty(avg_params_field)
        s = struct('K',NaN,'pi',[],'theta',[]);
    else
        s = avg_params_field{1};
    end
elseif isstruct(avg_params_field)
    s = avg_params_field;
else
    s = struct('K',NaN,'pi',[],'theta',[]);
end
end
function A = avg_from_cells(cells, model)
% AVG_FROM_CELLS (robust to varying K across datasets)
% Returns NaN-padded, NaN-robust averages of {pi, theta} across cells.
% model: 'poiss', 'pg', or 'pl'

    if isempty(cells) || all(cellfun(@isempty,cells))
        A = struct('mean_pi',NaN,'mean_theta',NaN);
        return;
    end

    % --- find K per dataset and overall Kmax ---
    K_list = [];
    for i = 1:numel(cells)
        s = cells{i};
        if isempty(s) || ~isfield(s,'pi') || isempty(s.pi), continue; end
        K_list(end+1) = numel(s.pi); %#ok<AGROW>
    end
    if isempty(K_list)
        A = struct('mean_pi',NaN,'mean_theta',NaN);
        return;
    end
    Kmax = max(K_list);

    % --- allocate NaN-padded collectors ---
    Pis = NaN(Kmax, numel(cells));
    switch lower(model)
        case 'poiss'
            Ths = NaN(Kmax, numel(cells));
        otherwise  % 'pg' or 'pl'
            Ths = NaN(2*Kmax, numel(cells));
    end

    % --- fill columns with available entries, pad the rest with NaN ---
    for i = 1:numel(cells)
        s = cells{i};
        if isempty(s) || ~isfield(s,'pi') || ~isfield(s,'theta')
            continue;
        end

        % pi
        pi_i = s.pi(:);
        ki = numel(pi_i);
        ki = min(ki, Kmax);          % safety
        if ki >= 1
            Pis(1:ki, i) = pi_i(1:ki);
        end

        % theta
        switch lower(model)
            case 'poiss'
                th_i = s.theta(:);
                if ~isempty(th_i)
                    th_i = th_i(:);
                    ki_th = min(numel(th_i), Kmax);
                    Ths(1:ki_th, i) = th_i(1:ki_th);
                end
            otherwise
                th_i = s.theta;
                if isvector(th_i), th_i = reshape(th_i(:), 2, []); end
                if ~isempty(th_i) && size(th_i,1)==2
                    ki_th = min(size(th_i,2), Kmax);
                    flat = th_i(:,1:ki_th); flat = flat(:);  % 2*ki_th × 1
                    Ths(1:2*ki_th, i) = flat;
                end
        end
    end

    % --- compute NaN-robust means ---
    mean_pi = mean(Pis, 2, 'omitnan');

    switch lower(model)
        case 'poiss'
            mean_theta = mean(Ths, 2, 'omitnan');
        otherwise
            tmp = mean(Ths, 2, 'omitnan');   % 2*Kmax × 1
            mean_theta = reshape(tmp, 2, []);% 2 × Kmax
    end

    A = struct('mean_pi', mean_pi, 'mean_theta', mean_theta);
end

