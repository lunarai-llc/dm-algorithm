% This is Lema Image Analysis with Poisson Mixture


%%
load image_a.mat

Y_c = rgb2gray(cdata);
Y = double(Y_c);

Y = Y(:);


%% plot original image Y

imshow(uint8(reshape(Y, 778, 780)))


%% choose # of K: EM algorithm

out_em = selectK_DIC(Y, ...
    'do_dic', false, ...              % skip Hellinger/VNED DIC
    'do_em',  true,  ...              % run EM + BIC/AIC
    'Kmax', 6, ...
    'em_nStarts', 10, ...             % multistart EM
    'em_maxIter', 1000, ...
    'em_tol', 1e-6, ...
    'em_init','quantile', ...
    'seed', 1);


% return K=3   ==> use EM 3 component below           

%% choose # of K: HMIX algorithm

out_hd = selectK_DIC(Y,'method','hd','penalty','bic','Kmax',6,...
                  'nStarts',10,'seed',1,'dic_solver','auto',...
                  'em_nStarts',10,'do_dic',true,'do_em',true);
              
% return K=3   ==> use HMIX 3 compnent below           

%% choose # of K: VNEDMIX algorithm
out_vned = selectK_DIC(Y, ...
    'method','vned', ...           % use VNED disparity
    'penalty','bic', ...           % b(n)=log(n)/2 (BIC-style)
    'Kmax',5, ...
    'nStarts',10, ...              % multistart for robustness
    'dic_maxIter',800, 'dic_tol',1e-8, ...
    'do_dic',true,  ...            % run VNEDMIX+DIC
    'do_em',false, ...             % (optional) skip EM/BIC here
    'seed',1);

% return K=3   ==> use VNEDMIX 3 compnent below    

%% EM 2 component
[EMIC2, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
    Y, 2, 'LabelLevels',[1 150], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC2.BIC

%% EM 3 component
[EMIC3, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
     Y, 3, 'LabelLevels',[1 100 200], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC3.BIC

%% EM 4 component
[EMIC4, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
     Y, 4, 'LabelLevels',[1 100 150 225], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC4.BIC

%% EM 5 component
[EMIC5, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
    Y, 5, 'LabelLevels', [1 50 150 175 225], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC5.BIC

%% HMIX 2 component
[HIC2, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 2, 'LabelLevels',[1 150], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);

HIC2

%% HMIX 3 component
[HIC3, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 3, 'LabelLevels',[1 100 200], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);

HIC3

%% HMIX 4 component
[HIC4, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 4, 'LabelLevels',[1 100 150 225], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);

HIC4
%% HMIX 5 component

[HIC5, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 5, 'LabelLevels',[1 50 150 175 225], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);

HIC5

%% VNEDMIX 2 component
[VNEDIC2, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 2, 'LabelLevels',[1 150], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);
VNEDIC2

%% VNEDMIX 3 component
[VNEDIC3, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 3, 'LabelLevels',[1 100 200], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);

VNEDIC3
%% VNEDMIX 4 component
[VNEDIC4, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 4, 'LabelLevels',[1 100 150 225], 'AutoImage', true, ...
    'MaxIter',800, 'Tol',1e-8, 'Init','quantile');

VNEDIC4

%% VNEDMIX 5 component
[VNEDIC5, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 5, 'LabelLevels',[1 50 150 175 225], 'AutoImage', true, ...
    'MaxIter',800, 'Tol',1e-8, 'Init','quantile');

VNEDIC5

%% Build BIC table across algorithms and K = 1..5
% Assumes Y exists. Uses/creates: EMICk, HICk, VNEDICk for k=1..5
Kvec = 1:5;

% Helpers to fetch or compute results (and close figures)
get_EM_BIC = @(k) fetch_or_run( ...
    sprintf('EMIC%d',k), @() em_bic_plot(Y, k, 'MaxIter',800,'Tol',1e-8,'Init','quantile'), 'BIC');

get_HMIX_BIC = @(k) fetch_or_run( ...
    sprintf('HIC%d',k), @() hmix_hic_plot(Y, k, 'MaxIter',800,'Tol',1e-8,'Init','quantile'), 'BIC');

get_VNED_BIC = @(k) fetch_or_run( ...
    sprintf('VNEDIC%d',k), @() vnedmix_dic_plot(Y, k, 'MaxIter',800,'Tol',1e-8,'Init','quantile'), 'BIC');

% Collect BICs (NaN if something goes wrong)
bic_em   = arrayfun(@(k) get_EM_BIC(k),   Kvec);
bic_hmix = arrayfun(@(k) get_HMIX_BIC(k), Kvec);
bic_vned = arrayfun(@(k) get_VNED_BIC(k), Kvec);

% Build table
RowNames = {'# of components','BIC','HIC','VNEDIC'};
VarNames = arrayfun(@(k) sprintf('K%d',k), Kvec, 'UniformOutput', false);

data = [
    Kvec;          % first row: 1..5
    bic_em;        % EM BIC
    bic_hmix;      % HMIX HIC with BIC-style penalty
    bic_vned       % VNEDMIX DIC with BIC-style penalty
];

T = array2table(data, 'RowNames', RowNames, 'VariableNames', VarNames)



%% Parameter table across methods (EM / HDMIX / VNEDMIX) and K = 3..5
% Requirements on path: em_bic_plot.m, hmix_hic_plot.m, vnedmix_dic_plot.m
% Input: Y (vector of nonnegative integers)

Kset    = [2 3 4 5];
methods = {'EM','HDMIX','VNEDMIX'};
Kmax    = max(Kset);

% Build parameter column names: pi1 lam1 pi2 lam2 ... pi_{Kmax-1} lam_{Kmax-1} lam_{Kmax}
paramNames = {};
for k = 1:Kmax
    if k < Kmax, paramNames{end+1} = sprintf('pi%d',k); end %#ok<SAGROW>
    paramNames{end+1} = sprintf('lambda%d',k); %#ok<SAGROW>
end

% Preallocate containers
nRows = numel(Kset) * numel(methods);
MethodCol = cell(nRows,1);
KCol      = zeros(nRows,1);
DataMat   = nan(nRows, numel(paramNames));

row = 0;
for kk = 1:numel(Kset)
    K = Kset(kk);

    for mm = 1:numel(methods)
        method = methods{mm};
        row = row + 1;

        % --- run estimator (no plotting kept) ---
        switch method
            case 'EM'
                [~, hFig, ~, pi_hat, lambda_hat] = em_bic_plot( ...
                    Y, K, 'MaxIter',800, 'Tol',1e-8, 'Init','quantile');
            case 'HDMIX'
                [~, hFig, ~, pi_hat, lambda_hat] = hmix_hic_plot( ...
                    Y, K, 'MaxIter',800, 'Tol',1e-8, 'Init','quantile');
            case 'VNEDMIX'
                [~, hFig, ~, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
                    Y, K, 'MaxIter',800, 'Tol',1e-8, 'Init','quantile');
        end
        if ishghandle(hFig), close(hFig); end

        % --- pack into one row: [pi1 lam1 pi2 lam2 ... pi_{K-1} lam_{K} , pad to Kmax] ---
        DataMat(row,:) = packParams(pi_hat, lambda_hat, K, Kmax);

        MethodCol{row} = method;
        KCol(row)      = K;
    end
end

% Build the table (first cols: K and Method, then parameter columns)
T_params = array2table(DataMat, 'VariableNames', paramNames);
T_params = addvars(T_params, KCol, MethodCol, 'Before', 1, ...
                   'NewVariableNames', {'K','Method'});

% Nice viewing order
T_params = sortrows(T_params, {'K','Method'});

% Show the table
disp(T_params);

% Optional: write to CSV or LaTeX
% writetable(T_params, 'param_estimates.csv');
% (For LaTeX you can format yourself or use any table writer you prefer.)

%% ------- helper: interleave params and pad to Kmax -------
function rowvec = packParams(pi_hat, lambda_hat, K, Kmax)
% Returns 1×(2*Kmax-1) row: [pi1 lam1 pi2 lam2 ... pi_{Kmax-1} lam_{Kmax}]
    out = [];
    for j = 1:Kmax
        if j < K
            out = [out, pi_hat(j)];           %#ok<AGROW>
        elseif j >= K && j < Kmax
            out = [out, NaN];                 %#ok<AGROW> % pi_j not present (or j>=K)
        end
        if j <= K
            out = [out, lambda_hat(j)];       %#ok<AGROW>
        else
            out = [out, NaN];                 %#ok<AGROW>
        end
    end
    % Remove the extra NaN inserted when j==Kmax and j>=K (we designed names as 2*Kmax-1)
    % Our loop already matches the 2*Kmax-1 layout, so nothing else to trim.
    rowvec = out;
end


%% --- Helper function: fetch or compute (and close plots) ---
function val = fetch_or_run(varname, runner, fieldname)
    % If varname exists (struct with .(fieldname)), use it
    if evalin('base', sprintf('exist(''%s'',''var'')', varname))
        S = evalin('base', varname);
        if isstruct(S) && isfield(S, fieldname)
            val = S.(fieldname);
            return;
        end
    end
    % Otherwise, compute; close figure; stash to base; return field
    try
        [S, hFig] = runner();               % #ok<NASGU> S is first output struct
        % Close the figure if it exists
        if exist('hFig','var') && ishghandle(hFig), close(hFig); end
        % Save the struct under varname in base workspace
        assignin('base', varname, S);
        if isfield(S, fieldname)
            val = S.(fieldname);
        else
            val = NaN;
        end
    catch
        val = NaN;
    end
end


