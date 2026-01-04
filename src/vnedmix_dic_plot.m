function [DIC, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot(Y, K, varargin)
% VNEDMIX_DIC_PLOT  VNEDMIX for Poisson mixtures + plot (no title) + DIC + parameter outputs
%   [DIC, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot(Y, K, ...
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
%   'MaxIter'    : max VNEDMIX iterations (default 800)
%   'Tol'        : convergence tolerance on [pi; lambda] (default 1e-8)
%   'Init'       : 'quantile' (default) | 'random'
%
% Outputs:
%   DIC : struct with fields (Disparity IC using VNED disparity):
%         .D   = VNED disparity = sum_i g_i * exp( - f_mix(YY_i) / g_i )
%         .AIC = D + (1/n)        * v(K)      (v(K)=2K-1 for Poisson mix)
%         .BIC = D + (log(n)/(2n))* v(K)
%   hFig   : handle to the generated figure
%   params : struct with fields: .pi, .lambda, .support, .g_n, .iters
%   pi_hat     : Kx1 mixing weights (VNEDMIX estimates)
%   lambda_hat : Kx1 Poisson rates  (VNEDMIX estimates)
%
% Notes:
%   - VNEDMIX E-step uses standard responsibilities on the observed support.
%   - M-step (for each k) minimizes in log-lambda:
%         NED_k(t) = sum_i exp( - pi_k * Pois(YY_i; exp(t)) / (g_i * r_{ik}) ) * g_i * r_{ik}
%     and updates weights via:
%         vned_k = sum_i exp( - pi_k * Pois(YY_i; ?_k) / (g_i * r_{ik}) ) * pi_k * Pois(YY_i; ?_k)
%         pi_k   = vned_k / sum_j vned_j
%   - Requires Optimization Toolbox (fminunc) and Statistics Toolbox (poisspdf, tabulate).
%   - MATLAB R2019a–compatible. Optimizer display disabled. No titles/captions.

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

% -------- VNEDMIX iterations --------
opts_fmin = optimoptions(@fminunc, 'Display','off', ...
    'Algorithm','quasi-newton','MaxIterations',500,'MaxFunctionEvaluations',5e4);

prev  = [pi; lam];
iters = 0;
while iters < S.MaxIter
    iters = iters + 1;

    % E-step: responsibilities on support
    % Pk(i,k) = Pois(YY_i; lam_k)
    Pk   = poisspdf(repmat(YY,1,K), repmat(lam(:).',numel(YY),1));   % |YY|×K
    rnum = bsxfun(@times, Pk, pi(:).');                              % |YY|×K
    rden = max(sum(rnum,2), 1e-15);
    R    = bsxfun(@rdivide, rnum, rden);                             % |YY|×K

    % M-step: per-component VNED objective in log-lambda
    lam_new  = lam;
    vned_val = zeros(K,1);
    for k = 1:K
        rk = R(:,k);

        % NED_k(t): sum_i exp( - pi_k * Pois(YY_i; exp(t)) / (g_i * r_{ik}) ) * g_i * r_{ik}
        NEDk = @(t) sum( exp( - pi(k) .* poisspdf(YY, max(1e-12,exp(t))) ...
                               ./ max(1e-15, g_n .* max(1e-12,rk)) ) ...
                         .* g_n .* max(1e-12, rk) );
        tk0 = log(max(1e-8, lam(k)));
        [tk_hat, ~] = fminunc(NEDk, tk0, opts_fmin);
        lam_new(k) = max(1e-8, exp(tk_hat));

        % vned_k for weight update (evaluate at new ?_k)
        Pk_k  = poisspdf(YY, lam_new(k));
        expo  = exp( - pi(k) .* Pk_k ./ max(1e-15, g_n .* max(1e-12,rk)) );
        vned_val(k) = sum( expo .* pi(k) .* Pk_k );
    end

    % Weight update: pi_k ? vned_k
    v = max(vned_val, 0);
    if sum(v)==0 || any(~isfinite(v)), v = ones(K,1); end
    pi_new = v / sum(v);

    theta_new = [pi_new; lam_new];
    if max(abs(theta_new - prev)) <= S.Tol
        pi = pi_new; lam = lam_new; break;
    end
    pi = pi_new; lam = lam_new; prev = theta_new;
end

% -------- Final mixture pmf and VNED DIC --------
Pk_final = poisspdf(repmat(YY,1,K), repmat(lam(:).',numel(YY),1));  % |YY|×K
f_mix    = max(Pk_final * pi(:), 1e-15);                            % |YY|×1

% VNED disparity: D = sum_i g_i * exp( - f_i / g_i )
ratio = f_mix ./ max(1e-15, g_n);
D_vned = sum( g_n .* exp( - ratio ) );

vK      = 2*K - 1;                  % Poisson mixture params
DIC.D   = D_vned;
DIC.AIC = D_vned + (1/n) * vK;
DIC.BIC = D_vned + (log(n)/2 / n) * vK;

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
    legend([{'Empirical pmf g_n','Mixture f'} ...
            arrayfun(@(k)sprintf('Component %d',k),1:K,'UniformOutput',false)], ...
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
