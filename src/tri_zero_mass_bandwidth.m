

function [c_hat, info] = tri_zero_mass_bandwidth(Y, a, renorm)
% TRI_ZERO_MASS_BANDWIDTH  Bandwidth (shape) selection for discrete triangle kernel
% by zero-mass matching with fixed span a.
%
%   [c_hat, info] = tri_zero_mass_bandwidth(Y, a, renorm)
%
% Inputs:
%   Y      : column or row vector of nonnegative integer counts
%   a      : fixed span (half-width) of triangle kernel support (positive integer)
%   renorm : (optional) true/false, boundary renormalization (default: true)
%
% Output:
%   c_hat  : selected shape/bandwidth parameter c > 0
%   info   : struct with details (p0, feasible, bracket, g0_at_c, residual, etc.)
%
% Method:
%   Choose c so that g(0; c) = mean(Y==0),
%   where g(0; c) = (1/n) * sum_i K_{Y_i, c}(0) with triangle kernel of span a.
%
% Notes:
%   - If there are no observations with Y <= a, then g(0; c) == 0 for all c,
%     so exact matching is only possible if mean(Y==0) == 0. We return c=1.
%   - If a sign-change bracket cannot be found, we minimize |g(0; c) - p0|
%     over c in [1e-3, 20] as a fallback.

    if nargin < 3 || isempty(renorm), renorm = true; end
    Y = Y(:);
    if any(Y < 0) || any(abs(Y - round(Y)) > 0)
        error('Y must contain nonnegative integers.');
    end
    if ~(isscalar(a) && a == round(a) && a >= 1)
        error('a must be a positive integer.');
    end

    n   = numel(Y);
    p0  = mean(Y == 0);

    % counts n_y for y = 0..a (only these can contribute to mass at zero)
    maxy   = max(max(Y), a);
    counts = accumarray(Y+1, 1, [maxy+1, 1]);  % 1-based indexing
    ny     = counts(1:a+1);                    % n_y for y=0..a

    info = struct();
    info.p0 = p0;

    if sum(ny) == 0
        % No Y <= a => g0(c) = 0 for all c; exact match only if p0==0
        c_hat = 1;
        info.feasible = (p0 == 0);
        info.note = 'No observations with Y<=a; g(0;c)=0 for all c.';
        info.g0_at_c = 0;
        info.residual = abs(0 - p0);
        return
    end

    % g0(c): estimated mass at zero as a function of c
    yvec = (0:a).';             % column vector 0..a
    k    = (1:a).';             % distances 1..a (for power sums)
    nn   = n;                   % capture for nested function

    function g0 = g0_eval(c)
        A  = (a+1)^c;           % (a+1)^c
        kc = k.^c;              % k^c, k=1..a
        Sa = sum(kc);           % sum_{k=1}^a k^c
        if renorm
            % Boundary-renormalized normalizer:
            % P_y(a,c) = (a+y+1)A - [sum_{k=1}^a k^c + sum_{k=1}^y k^c]
            Ty  = [0; cumsum(kc)];                 % Ty(y) = sum_{k=1}^y k^c, Ty(0)=0
            Py  = (a + yvec + 1) .* A - (Sa + Ty); % denom per y = 0..a
            num = A - yvec.^c;                     % numerator per y (|0-y|=y)
            term = ny .* (num ./ Py);
            g0 = sum(term) / nn;
        else
            % Global normalizer (slightly faster, less boundary-accurate)
            P  = (2*a + 1) * A - 2 * Sa;          % P(a,c)
            num = ny .* (A - yvec.^c);
            g0 = sum(num) / (nn * P);
        end
    end

    f = @(c) g0_eval(c) - p0;

    % Try to find a sign-change bracket by scanning a coarse grid
    c_grid = [0.05 0.1 0.2 0.5 1 1.5 2 3 4 6 8 10 15 20];
    f_vals = arrayfun(f, c_grid);
    idx    = find(f_vals(1:end-1) .* f_vals(2:end) <= 0, 1, 'first');

    if ~isempty(idx) && all(isfinite([f_vals(idx), f_vals(idx+1)]))
        c_lo = c_grid(idx);  c_hi = c_grid(idx+1);
        c_hat = fzero(f, [c_lo, c_hi]);
        info.bracket = [c_lo, c_hi];
        info.feasible = true;
    else
        % No bracket found: minimize absolute mismatch over a wide interval
        obj   = @(c) abs(f(c));
        lb = 1e-3; ub = 20;
        c_hat = fminbnd(obj, lb, ub);
        info.bracket = [NaN, NaN];
        info.feasible = false;
        info.note = 'No sign-change bracket; returned minimizer of |g(0;c)-p0|.';
    end

    info.g0_at_c = g0_eval(c_hat);
    info.residual = abs(info.g0_at_c - p0);
end
