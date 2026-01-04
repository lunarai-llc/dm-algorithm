function [HIC, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot(Y, K, varargin)
% HMIX_HIC_PLOT  HMIX for Poisson mixtures + plot (no title) + HIC + parameter outputs
%   [HIC, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot(Y, K, ...
%       'ImageSize',[H W], 'AutoImage',true, 'LabelLevels',levels, ...
%       'MaxIter',800, 'Tol',1e-8, 'Init','quantile')
%
% Inputs:
%   Y          : vector of nonnegative integers (counts)
%   K          : number of mixture components (integer >=1)
%
% Name-Value (optional):
%   'ImageSize'  : [H W] — if provided, produce segmented image plot
%   'AutoImage'  : true (default). If ImageSize not given, guess H×W from numel(Y)
%   'LabelLevels': vector of length K with intensities for labels 1..K
%   'MaxIter'    : max HMIX iterations (default 800)
%   'Tol'        : convergence tolerance on [pi; lambda] (default 1e-8)
%   'Init'       : 'quantile' (default) | 'random'
%
% Outputs:
%   HIC : struct with fields:
%         .H2  = 2*(1 - sum sqrt(g .* f_mix))
%         .AIC = H2 + (1/n)        * v(K)      (v(K)=2K-1)
%         .BIC = H2 + (log(n)/(2n))* v(K)
%   hFig   : handle to the generated figure
%   params : struct with fields: .pi, .lambda, .support, .g_n, .iters
%   pi_hat     : Kx1 mixing weights (HMIX estimates)
%   lambda_hat : Kx1 Poisson rates  (HMIX estimates)
%
% Notes:
%   - Requires Optimization Toolbox (fminunc) and Statistics Toolbox (poisspdf, tabulate).
%   - R2019a-compatible. M-step optimizes in log-lambda. No titles/captions. Optimizer display disabled.

% -------- Parse & validate inputs --------
ip = inputParser;
ip.addRequired('Y', @(x) isnumeric(x) && isvector(x) && all(x>=0) && all(mod(x,1)==0));
ip.addRequired('K', @(x) isnumeric(x) && isscalar(x) && x>=1 && mod(x,1)==0);
ip.addParameter('ImageSize', [], @(x) isnumeric(x) && numel(x)==2 && all(x>0));
ip.addParameter('AutoImage', true, @(x) islogical(x) && isscalar(x));
ip.addParameter('LabelLevels', [], @(x) isnumeric(x) && isvector(x));
ip.addParameter('MaxIter', 800, @(x) isnumeric(x) && isscalar(x) && x>=10);
ip.addParameter('Tol', 1e-8, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('Init', 'quantile', @(s) ischar(s) && ismember(lower(s),{'quantile','random'}));
ip.parse(Y, K, varargin{:});
S = ip.Results;

Y = double(Y(:));
n = numel(Y);

if ~isempty(S.LabelLevels) && numel(S.LabelLevels) ~= K
    error('LabelLevels must have length K (%d).', K);
end

% -------- Empirical pmf g_n on observed support --------
a   = tabulate(Y);
YY  = a(:,1);           % support (unique sorted)
g_n = a(:,3)/100;       % empirical pmf (>0 on support)

% -------- Initialize (pi, lambda) --------
switch lower(S.Init)
    case 'quantile'
        pi  = ones(K,1)/K;
        qs  = linspace(5,95,K);
        lam = max(1e-3, prctile(Y, qs).');
        if any(~isfinite(lam))
            mY  = max(1e-3, mean(Y));
            lam = mY * linspace(0.5, 1.5, K).';
        end
    case 'random'
        a0  = rand(K,1) + 0.1;
        pi  = a0 / sum(a0);
        mY  = max(1e-3, mean(Y));
        sY  = max(1e-3, std(Y));
        lam = max(1e-3, mY + sY*randn(K,1));
end

% -------- HMIX iterations --------
% Display is hard-coded OFF (silent) for fminunc
opts_fmin = optimoptions(@fminunc, 'Display','off', ...
    'Algorithm','quasi-newton','MaxIterations',500,'MaxFunctionEvaluations',5e4);

prev  = [pi; lam];
iters = 0;
while iters < S.MaxIter
    iters = iters + 1;

    % E-step: responsibilities on support
    Pk   = poisspdf(repmat(YY,1,K), repmat(lam(:).',numel(YY),1));  % |YY|×K
    rnum = bsxfun(@times, Pk, pi(:).');                             % |YY|×K
    rden = max(sum(rnum,2), 1e-15);
    R    = bsxfun(@rdivide, rnum, rden);                            % |YY|×K

    % M-step: per-component minimization in log-lambda
    lam_new = lam;
    f_k     = zeros(K,1);
    for k = 1:K
        rk = R(:,k);
        Hk = @(t) -sum( sqrt( max(1e-12, rk).*g_n .* poisspdf(YY, max(1e-12,exp(t))) ) );
        tk0 = log(max(1e-8, lam(k)));
        [tk_hat, fk] = fminunc(Hk, tk0, opts_fmin);
        lam_new(k) = max(1e-8, exp(tk_hat));
        f_k(k)     = fk;   % negative is fine; we use fk^2
    end

    % Weight update: pi_k ? f_k^2
    w = f_k.^2;
    if ~all(isfinite(w)), w = ones(K,1); end
    pi_new = max(eps, w / sum(w));

    theta_new = [pi_new; lam_new];
    if max(abs(theta_new - prev)) <= S.Tol
        pi = pi_new; lam = lam_new; break;
    end
    pi = pi_new; lam = lam_new; prev = theta_new;
end

% -------- Final mixture pmf and HIC --------
Pk_final = poisspdf(repmat(YY,1,K), repmat(lam(:).',numel(YY),1));  % |YY|×K
f_mix    = max(Pk_final * pi(:), 1e-15);                            % |YY|×1

BC = sum( sqrt(g_n .* f_mix) );
H2 = 2 * (1 - BC);

vK      = 2*K - 1;                  % Poisson mixture params
HIC.H2  = H2;
HIC.AIC = H2 + (1/n) * vK;
HIC.BIC = H2 + (log(n)/2 / n) * vK;

% -------- Decide image size (auto by default) --------
imgSize = S.ImageSize;
if isempty(imgSize) && S.AutoImage
    imgSize = guess_image_size(n);      % returns a 1×2 vector
end

% -------- Plot (no title) --------
if ~isempty(imgSize)
    % MAP component on support (final posteriors)
    rnum = bsxfun(@times, Pk_final, pi(:).');        % |YY|×K
    rden = max(sum(rnum,2), 1e-15);
    Rfin = bsxfun(@rdivide, rnum, rden);
    [~, lab_YY] = max(Rfin, [], 2);                  % |YY|×1

    % Propagate labels to all entries of Y
    Y_label = zeros(n,1);
    for i = 1:numel(YY)
        idx = find(Y == YY(i));
        if ~isempty(idx), Y_label(idx) = lab_YY(i); end
    end

    % Labels ? intensities
    if ~isempty(S.LabelLevels)
        levels = S.LabelLevels(:);                    % exact mapping
    else
        levels = round(linspace(20,235,K)).';
    end
    new_Y = levels(Y_label);                          % n×1

    img = uint8(reshape(new_Y, imgSize(1), imgSize(2)));
    hFig = figure;
    if exist('imshow','file')
        imshow(img);
    else
        imagesc(img); axis image off; colormap(gray(256));
    end
else
    % PMF plot: empirical vs mixture and components (no title)
    hFig = figure;
    bar(YY, g_n, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none'); hold on;
    plot(YY, f_mix, 'LineWidth', 2);
    for k = 1:K
        plot(YY, pi(k)*Pk_final(:,k), '--', 'LineWidth', 1);
    end
    grid on; xlim([min(YY) max(YY)]);
    xlabel('Value'); ylabel('Probability / pmf');
    legend([{'Empirical pmf g_n','Mixture f'} arrayfun(@(k)sprintf('Component %d',k),1:K,'UniformOutput',false)], ...
           'Location','best');
end

% -------- Outputs --------
pi_hat     = pi(:);
lambda_hat = lam(:);
params.pi      = pi_hat;
params.lambda  = lambda_hat;
params.support = YY;
params.g_n     = g_n;
params.iters   = iters;
end

% ---------- helpers ----------
function sz = guess_image_size(n)
% Find factor pair (h,w) closest to square; return [h w]
s = floor(sqrt(n));
for d = s:-1:1
    if mod(n,d)==0
        h = d; w = n/d;
        if h > w, t = h; h = w; w = t; end
        sz = [h, w];
        return
    end
end
% Fallback (n is prime): 1×n
sz = [1, n];
end
