function out = mix_mc_identifyK(B, n, model, trueParams, varargin)
% MIX_MC_IDENTIFYK
% Monte Carlo with R=1 and B datasets, parallelized over B (datasets).
% For each dataset:
%   1) Simulate from Poisson / PG / PL truth.
%   2) Run split-select-estimate (R=1, no inner parallel).
%   3) Record modal-K, per-dataset estimates for EM / HD / VNED.
%
% Outputs:
%   out.prop_correct         -> % of datasets with correctly identified K (per method)
%   out.split_correct_pct    -> [% over splits C that chose true K] per dataset (per method)
%   out.khat                 -> Khat per dataset (per method)
%   out.est                  -> per-dataset estimates (struct {K, pi, theta}) (per method)
%   out.avg_when_correct     -> average estimates over datasets with Khat==Ktrue (per method)
%   out.settings             -> echo of settings used
%
% SPEED CHOICES (default):
%   - dic_rule='stop' (early stop in K)
%   - modest nStarts and MaxIter (tweak via varargin)
%   - GH nodes for PL kept modest (gh_n=25~30)
%
% Example:
% out = mix_mc_identifyK(200, 2000, 'pl', ...
%        struct('pi',[0.3;0.7],'theta',[3 1; 0.5 0.5]), ...
%        'Kmax',5,'C',6,'nStarts',4,'MaxIter',120,'gh_n',25);

% ---------------- Parse inputs ----------------
ip = inputParser;
ip.addRequired('B', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('n', @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addRequired('model', @(s)ischar(s)&&ismember(lower(s),{'poiss','pg','pl'}));
ip.addRequired('trueParams', @(s)isstruct(s)&&isfield(s,'pi')&&isfield(s,'theta'));

ip.addParameter('methods', {'em','hd','vned'}, @(c)iscellstr(c)&&~isempty(c));
ip.addParameter('Kmax',   5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('C',      5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('penalty','bic', @(s)ischar(s)&&ismember(lower(s),{'bic','aic'}));
ip.addParameter('dic_rule','stop', @(s)ischar(s)&&ismember(lower(s),{'stop','min'}));
ip.addParameter('nStarts',4,@(x)isnumeric(x)&&isscalar(x)&&x>=1);        % lean default
ip.addParameter('MaxIter',100,@(x)isnumeric(x)&&isscalar(x)&&x>=10);     % lean default
ip.addParameter('Tol',1e-5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('init','quantile',@(s)ischar(s)&&ismember(lower(s),{'quantile','random'}));
ip.addParameter('seed',1,@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('gh_n',20,@(x)isnumeric(x)&&isscalar(x)&&x>=5);          % lean default for PL
ip.parse(B,n,model,trueParams,varargin{:});
S = ip.Results;

methods = lower(S.methods(:))';
doEM   = any(strcmp(methods,'em'));
doHD   = any(strcmp(methods,'hd'));
doVNED = any(strcmp(methods,'vned'));

% truth and wrappers
pi_true = S.trueParams.pi(:);
Ktrue   = numel(pi_true);
truth_wrap = struct(); truth_wrap.(lower(S.model)) = S.trueParams;

% ---------------- Preallocate parfor-sliced outputs ----------------
khat_EM   = NaN(B,1);
khat_HD   = NaN(B,1);
khat_VNED = NaN(B,1);

% split correctness (% of splits selecting true K) per dataset
splitcorr_EM   = NaN(B,1);
splitcorr_HD   = NaN(B,1);
splitcorr_VNED = NaN(B,1);

% Per-dataset estimates as cells of structs {K, pi, theta}
est_EM   = cell(B,1);
est_HD   = cell(B,1);
est_VNED = cell(B,1);

% ---------------- Spin up a pool (parallel over B) ----------------
if license('test','Distrib_Computing_Toolbox')
    p = gcp('nocreate'); if isempty(p), parpool('local'); end %#ok<PPOOL>
end

% ---------------- Monte Carlo over B datasets (PARFOR) -------------
parfor b = 1:B
    rng(S.seed + b);

    % 1) simulate one dataset
    Y = generate_dataset(S.n, lower(S.model), S.trueParams);

    % 2) run split-select-estimate with R=1 (disable inner parallel!)
    outb = mix_split_select_estimate( ...
        Y, ...
        'models',    {lower(S.model)}, ...
        'methods',   methods, ...
        'Kmax',      S.Kmax, ...
        'C',         S.C, ...
        'R',         1, ...                 % fixed
        'penalty',   S.penalty, ...
        'dic_rule',  S.dic_rule, ...
        'nStarts',   S.nStarts, ...
        'MaxIter',   S.MaxIter, ...
        'Tol',       S.Tol, ...
        'init',      S.init, ...
        'seed',      S.seed + b, ...
        'trueK',     Ktrue, ...
        'trueParams',truth_wrap, ...
        'gh_n',      S.gh_n, ...
        'UseParallel',false, ...        
        'ComputeHessian', false);

    Rm = outb.results.(lower(S.model));

    % 3) Extract Khat, estimates, and split-wise correctness
    if doEM && isfield(Rm,'em')
        kh = scalar_or_first(Rm.em.Khat_mode);
        khat_EM(b) = kh;
        est_EM{b}  = first_avg_struct(Rm.em.avg_params);
        splitcorr_EM(b) = split_correct_pct(Rm.em.Khat_counts, Ktrue, S.C);
    end
    if doHD && isfield(Rm,'hd')
        kh = scalar_or_first(Rm.hd.Khat_mode);
        khat_HD(b) = kh;
        est_HD{b}  = first_avg_struct(Rm.hd.avg_params);
        splitcorr_HD(b) = split_correct_pct(Rm.hd.Khat_counts, Ktrue, S.C);
    end
    if doVNED && isfield(Rm,'vned')
        kh = scalar_or_first(Rm.vned.Khat_mode);
        khat_VNED(b) = kh;
        est_VNED{b}  = first_avg_struct(Rm.vned.avg_params);
        splitcorr_VNED(b) = split_correct_pct(Rm.vned.Khat_counts, Ktrue, S.C);
    end
end

% ---------------- % correct K and averages (when correct) -----------
prop_correct = struct();
avg_when_correct = struct();

if doEM
   prop_correct.EM = 100*mean(khat_EM == Ktrue,'omitnan');
   avg_when_correct.EM = average_when_correct(est_EM, khat_EM, Ktrue, lower(S.model));
end
if doHD
   prop_correct.HD = 100*mean(khat_HD == Ktrue,'omitnan');
   avg_when_correct.HD = average_when_correct(est_HD, khat_HD, Ktrue, lower(S.model));
end
if doVNED
   prop_correct.VNED = 100*mean(khat_VNED == Ktrue,'omitnan');
   avg_when_correct.VNED = average_when_correct(est_VNED, khat_VNED, Ktrue, lower(S.model));
end

% ---------------- Package outputs -----------------------------------
out = struct();
out.prop_correct = prop_correct;

out.split_correct_pct = struct();
if doEM,   out.split_correct_pct.EM   = splitcorr_EM;   end
if doHD,   out.split_correct_pct.HD   = splitcorr_HD;   end
if doVNED, out.split_correct_pct.VNED = splitcorr_VNED; end

out.khat = struct();
if doEM,   out.khat.EM   = khat_EM;   end
if doHD,   out.khat.HD   = khat_HD;   end
if doVNED, out.khat.VNED = khat_VNED; end

out.est = struct();
if doEM,   out.est.EM   = est_EM;   end
if doHD,   out.est.HD   = est_HD;   end
if doVNED, out.est.VNED = est_VNED; end

out.avg_when_correct = avg_when_correct;

out.settings = S;
out.settings.trueK = Ktrue;
out.settings.trueParams = S.trueParams;
end

% ========================= Helpers =========================

function Y = generate_dataset(n, model, truth)
piK = truth.pi(:); K = numel(piK);
u = rand(n,1); edges = [0; cumsum(piK)];
z = zeros(n,1);
for k = 1:K, z = z + (u>edges(k) & u<=edges(k+1))*k; end
Y = zeros(n,1);
switch lower(model)
    case 'poiss'
        lam = truth.theta(:);
        for k = 1:K
            idx = (z==k);
            if any(idx), Y(idx) = poissrnd(lam(k), sum(idx), 1); end
        end
    case 'pg'
        a = truth.theta(1,:); b = truth.theta(2,:);
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

function s = first_avg_struct(avg_params_field)
% Handle cell/struct for avg_params when R=1
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

function v = scalar_or_first(x)
if isempty(x), v = NaN; return; end
if isscalar(x), v = x; else, v = x(1); end
end

function pct = split_correct_pct(Khat_counts_field, Ktrue, C)
% Return fraction (0..100) of splits that chose the true K
% Khat_counts is a 1x1 cell (R=1) of struct('K',K_unique,'count',counts)
if isempty(Khat_counts_field), pct = NaN; return; end
kc = Khat_counts_field;
if iscell(kc), kc = kc{1}; end
if ~isstruct(kc) || ~isfield(kc,'K') || ~isfield(kc,'count')
    pct = NaN; return;
end
ix = find(kc.K==Ktrue, 1);
if isempty(ix), num_ok = 0; else, num_ok = kc.count(ix); end
if nargin<3 || isempty(C) || ~isscalar(C), C = sum(kc.count); end
pct = 100 * num_ok / max(1,C);
end

function A = average_when_correct(ests, khat_vec, Ktrue, model)
ix = find(khat_vec == Ktrue);
if isempty(ix)
    A = struct('mean_pi', NaN(Ktrue,1), 'mean_theta', NaN*zeros(theta_shape(Ktrue,model)));
    return;
end
Pis_ok = NaN(Ktrue, numel(ix));
switch lower(model)
    case 'poiss', Ths_ok = NaN(Ktrue, numel(ix));
    otherwise,    Ths_ok = NaN(2*Ktrue, numel(ix));
end
tcount = 0;
for t = 1:numel(ix)
    s = ests{ix(t)};
    if isempty(s) || ~isfield(s,'pi') || ~isfield(s,'theta') || ~isfield(s,'K') || s.K~=Ktrue
        continue;
    end
    tcount = tcount + 1;
    Pis_ok(:,tcount) = s.pi(:);
    if strcmpi(model,'poiss')
        Ths_ok(:,tcount) = s.theta(:);
    else
        th = s.theta; if isvector(th), th = reshape(th(:),2,[]); end
        Ths_ok(1:2*Ktrue, tcount) = th(:);
    end
end
if tcount==0
    A = struct('mean_pi', NaN(Ktrue,1), 'mean_theta', NaN*zeros(theta_shape(Ktrue,model)));
    return;
end
Pis_ok = Pis_ok(:,1:tcount);
Ths_ok = Ths_ok(:,1:tcount);

mean_pi = mean(Pis_ok,2,'omitnan');
if strcmpi(model,'poiss')
    mean_theta = mean(Ths_ok,2,'omitnan');
else
    tmp = mean(Ths_ok,2,'omitnan'); tmp = reshape(tmp,2,[]);
    mean_theta = tmp;
end
A = struct('mean_pi',mean_pi,'mean_theta',mean_theta);
end

function T = theta_shape(K, model)
if strcmpi(model,'poiss'), T = [K,1]; else, T = [2,K]; end
end
