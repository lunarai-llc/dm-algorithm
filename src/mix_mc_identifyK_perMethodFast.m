function out = mix_mc_identifyK_perMethodFast(B, n, model, trueParams, varargin)
% MIX_MC_IDENTIFYK_PERMETHODFAST
% Monte Carlo over B datasets.
% For EACH method (EM / HD / VNED), we:
%   1) Select K̂ on a SUBSAMPLE (size n_sub) using that method with cheap settings
%      (sel_C = 1 by default; tiny nStarts/iters; no Hessian).
%   2) Refit ONCE on FULL data at that method's K̂ with fit_C splits (default 5).
% Parallelized over B (no nested parfor). Robust to optimizer errors.
%
% Output fields:
%   out.khat.(EM|HD|VNED)               -> [B x 1] K̂ per method
%   out.sel_counts.(EM|HD|VNED){b}      -> selection split count struct for that method
%   out.sel_frac_true.(EM|HD|VNED)      -> fraction of splits selecting true K (uses sel_C)
%   out.est.(EM|HD|VNED){b}             -> struct {K, pi, theta} from final full-data refit
%   out.avg_overall.(method)            -> average of estimates across all datasets (NaN-robust)
%   out.avg_when_correct.(method)       -> average across datasets where that method's K̂==trueK
%   out.prop_correct.(method)           -> % datasets where that method's K̂==trueK
%   out.settings                        -> echo of settings
%
% Requires: mix_split_select_estimate.m (recent version with ComputeHessian & ParallelAxis)

% ---------------- Parse inputs ----------------
ip = inputParser;
ip.addRequired('B', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('n', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('model', @(s)ischar(s)&&ismember(lower(s),{'poiss','pg','pl'}));
ip.addRequired('trueParams', @(s)isstruct(s) && isfield(s,'pi') && isfield(s,'theta'));

ip.addParameter('methods', {'em','hd','vned'}, @(c)iscellstr(c)&&~isempty(c));
ip.addParameter('Kmax',   5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

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
ip.parse(B,n,model,trueParams,varargin{:});
S = ip.Results;

methods = lower(S.methods(:))';
doEM   = any(strcmp(methods,'em'));
doHD   = any(strcmp(methods,'hd'));
doVNED = any(strcmp(methods,'vned'));

model = lower(S.model);
Ktrue = numel(trueParams.pi(:));

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
            run_one_dataset(b, S, model, trueParams, methods, Ktrue);
    end
else
    for b = 1:B
        rng(S.seed + b);
        [khat_EM(b), est_EM{b}, sel_counts_EM{b}, sel_frac_EM(b), ...
         khat_HD(b), est_HD{b}, sel_counts_HD{b}, sel_frac_HD(b), ...
         khat_VN(b), est_VN{b}, sel_counts_VN{b}, sel_frac_VN(b)] = ...
            run_one_dataset(b, S, model, trueParams, methods, Ktrue);
    end
end

% ---------------- % correct and averages (per method) ----------------
prop_correct = struct();
avg_overall  = struct();
avg_when_correct = struct();

if doEM
    prop_correct.EM = 100*mean(khat_EM==Ktrue,'omitnan');
    avg_overall.EM  = avg_from_cells(est_EM,  model);
    avg_when_correct.EM = avg_from_cells(est_EM(khat_EM==Ktrue), model);
end
if doHD
    prop_correct.HD = 100*mean(khat_HD==Ktrue,'omitnan');
    avg_overall.HD  = avg_from_cells(est_HD,  model);
    avg_when_correct.HD = avg_from_cells(est_HD(khat_HD==Ktrue), model);
end
if doVNED
    prop_correct.VNED = 100*mean(khat_VN==Ktrue,'omitnan');
    avg_overall.VNED  = avg_from_cells(est_VN, model);
    avg_when_correct.VNED = avg_from_cells(est_VN(khat_VN==Ktrue), model);
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
out.settings.trueK = Ktrue;
out.settings.trueParams = trueParams;
end

% ======================================================================
% ======================= One dataset (per method) =====================
% ======================================================================
function [K_em, estEM, C_em, F_em, ...
          K_hd, estHD, C_hd, F_hd, ...
          K_vn, estVN, C_vn, F_vn] = ...
          run_one_dataset(b, S, model, truth, methods, Ktrue)

% 0) simulate full data
Yfull = generate_dataset(S.n, model, truth);

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
    [K_em, C_em, F_em] = try_select(Ysel, model, 'em', S, b, Ktrue);
    estEM = try_refit(Yfull, model, 'em', S, b, K_em, truth);
end

% ---- HD selection on subsample (cheap) ----
if any(strcmp(methods,'hd'))
    [K_hd, C_hd, F_hd] = try_select(Ysel, model, 'hd', S, b, Ktrue);
    estHD = try_refit(Yfull, model, 'hd', S, b, K_hd, truth);
end

% ---- VNED selection on subsample (cheap) ----
if any(strcmp(methods,'vned'))
    [K_vn, C_vn, F_vn] = try_select(Ysel, model, 'vned', S, b, Ktrue);
    estVN = try_refit(Yfull, model, 'vned', S, b, K_vn, truth);
end
end

% ---------- selection helper (subsample; C=sel_C; no Hessian) ----------
function [Khat, Kcnt, fracTrue] = try_select(Ysel, model, method, S, b, Ktrue)
Khat = NaN; Kcnt = []; fracTrue = NaN;
try
    out_sel = mix_split_select_estimate( ...
        Ysel, ...
        'models', {model}, 'methods', {method}, ...
        'Kmax', S.Kmax, 'C', S.sel_C, 'R', 1, ...
        'penalty', S.penalty, 'dic_rule', S.sel_rule, ...
        'nStarts', S.sel_nStarts, 'MaxIter', S.sel_MaxIter, 'Tol', S.sel_Tol, ...
        'init', S.init, 'seed', S.seed + 101*b + method_seed_offset(method), ...
        'trueK', Ktrue, 'trueParams', struct(model, S.trueParams), ...
        'gh_n', S.gh_n, 'ComputeHessian', false, ...
        'UseParallel', false, 'ParallelAxis','off');

    Rm = out_sel.results.(model).(method);
    Khat = Rm.Khat_mode(1);
    Kcnt = Rm.Khat_counts{1};
    if ~isempty(Kcnt) && isfield(Kcnt,'K') && isfield(Kcnt,'count')
        pos = find(Kcnt.K==Ktrue,1);
        if ~isempty(pos)
            fracTrue = Kcnt.count(pos) / S.sel_C;
        else
            fracTrue = 0;
        end
    end
catch
    % leave NaNs/empty
end
end

% ---------- refit helper (full data; C=fit_C; no Hessian) ----------
function est = try_refit(Yfull, model, method, S, b, Khat, truth)
est = [];
if isnan(Khat) || Khat<1, return; end
try
    out_fit = mix_split_select_estimate( ...
        Yfull, ...
        'models', {model}, 'methods', {method}, ...
        'Kmax', Khat, 'C', S.fit_C, 'R', 1, ...
        'penalty', S.penalty, 'dic_rule', 'stop', ...
        'nStarts', S.fit_nStarts, 'MaxIter', S.fit_MaxIter, 'Tol', S.fit_Tol, ...
        'init', S.init, 'seed', S.seed + 313*b + method_seed_offset(method), ...
        'trueK', Khat, 'trueParams', struct(model, truth), ...
        'gh_n', S.gh_n, 'ComputeHessian', false, ...
        'UseParallel', false, 'ParallelAxis','off');
    est = safe_first_avg(out_fit.results.(model).(method).avg_params);
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
function Y = generate_dataset(n, model, truth)
piK = truth.pi(:); K = numel(piK);
u = rand(n,1); edges = [0; cumsum(piK)];
z = zeros(n,1);
for k = 1:K, z = z + (u>edges(k) & u<=edges(k+1))*k; end
Y = zeros(n,1);
switch model
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
        mu = truth.theta(1,:); sig = truth.theta(2,:);
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
if isempty(cells) || all(cellfun(@isempty,cells))
    A = struct('mean_pi',NaN,'mean_theta',NaN); return;
end
% infer K from first non-empty
k = NaN; 
for i=1:numel(cells)
    if ~isempty(cells{i}) && isfield(cells{i},'pi')
        k = numel(cells{i}.pi); break;
    end
end
if isnan(k), A = struct('mean_pi',NaN,'mean_theta',NaN); return; end

Pis = NaN(k, numel(cells));
switch model
    case 'poiss', Ths = NaN(k, numel(cells));
    otherwise,    Ths = NaN(2*k, numel(cells));
end
for i=1:numel(cells)
    s = cells{i};
    if isempty(s) || ~isfield(s,'pi') || ~isfield(s,'theta'), continue; end
    Pis(:,i) = s.pi(:);
    if strcmpi(model,'poiss')
        Ths(:,i) = s.theta(:);
    else
        th = s.theta; if isvector(th), th = reshape(th(:),2,[]); end
        Ths(1:2*k,i) = th(:);
    end
end
mean_pi = mean(Pis,2,'omitnan');
if strcmpi(model,'poiss')
    mean_theta = mean(Ths,2,'omitnan');
else
    tmp = mean(Ths,2,'omitnan'); tmp = reshape(tmp,2,[]);
    mean_theta = tmp;
end
A = struct('mean_pi',mean_pi,'mean_theta',mean_theta);
end
