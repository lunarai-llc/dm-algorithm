function T = summarize_pg_results_data_splitting(out)
% SUMMARIZE_PG_RESULTS
% Build a table for EM / HMIX / VNEDMIX:
% - % correctly identifying true K
% - averages and SDs for [pi1, a1, b1, a2, b2] when K is correct,
%   after 3*IQR trimming and with components aligned so comp 1 has larger mean (a/b).

methods = {'EM','HD','VNED'};
Ktrue = numel(out.settings.trueParams.pi(:));

rows = {};
for m = 1:numel(methods)
    tag = methods{m};
    if ~isfield(out.khat, tag) || ~isfield(out.est, tag), continue; end

    kh = out.khat.(tag);
    estCells = out.est.(tag);

    % overall % correct
    pct_correct = 100*mean(kh == Ktrue, 'omitnan');

    % collect estimates for runs with correct K
    idx = find(kh == Ktrue);
    P = [];   % columns: [pi1, a1, b1, a2, b2] after alignment
    for i = 1:numel(idx)
        s = estCells{idx(i)};
        if isempty(s) || ~isfield(s,'pi') || ~isfield(s,'theta') || ~isfield(s,'K'), continue; end
        if s.K ~= Ktrue, continue; end
        p = extract_pg_row(s);
        if any(~isfinite(p)) || any(p <= 0) || p(1) >= 1   % pi1 in (0,1)
            continue;
        end
        P(end+1,:) = p; %#ok<AGROW>
    end

    n_correct_valid = size(P,1);

    % 3*IQR trimming across all 5 columns
    if ~isempty(P)
        q1 = prctile(P,25,1);
        q3 = prctile(P,75,1);
        iqr_ = q3 - q1;
        lower = q1 - 3*iqr_;
        upper = q3 + 3*iqr_;
        keep = all(bsxfun(@ge,P,lower) & bsxfun(@le,P,upper), 2);
        Ptrim = P(keep,:);
    else
        Ptrim = P;
    end

    n_used     = size(Ptrim,1);
    n_trimmed  = n_correct_valid - n_used;

    mu = nan(1,5); sd = nan(1,5);
    if ~isempty(Ptrim)
        mu = mean(Ptrim,1,'omitnan');
        sd = std(Ptrim,0,1,'omitnan');
    end

    rows(end+1,:) = { tag, pct_correct, ...
        mu(1), sd(1), mu(2), sd(2), mu(3), sd(3), mu(4), sd(4), mu(5), sd(5), ...
        n_correct_valid, n_used, n_trimmed }; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', { ...
    'Method','PctCorrect', ...
    'pi1_mean','pi1_sd','a1_mean','a1_sd','b1_mean','b1_sd','a2_mean','a2_sd','b2_mean','b2_sd', ...
    'N_correct_valid','N_used_afterTrim','N_trimmed'});

disp(T);
end

% ---- helpers ----
function row = extract_pg_row(s)
% Align so component 1 has the larger mean a/b and return [pi1 a1 b1 a2 b2].
pi = s.pi(:); if ~isempty(pi), pi = pi / sum(pi); end
th = s.theta;
if isvector(th), th = reshape(th(:),2,[]); end
a = th(1,:); b = th(2,:);
mu = a ./ b;
[~, ord] = sort(mu, 'descend');   % comp-1 = larger mean
pi = pi(ord);
a  = a(ord);
b  = b(ord);
row = [pi(1), a(1), b(1), a(2), b(2)];
end
