%% example_riesz.m
%  Solves an equation involving the Riesz fractional integral operator:
%   (I + I_l + I_r) u = f, where I_l and I_r are left/right FIOs of order 1/2.
%  The Riesz operator is (I_l + I_r) / (2 cos(pi*mu/2)).

N = 10:2:200;
err = [];
for j = 1:length(N)
    n = N(j);
    mu = 1/2;
    b = 3.9;
    phi = @(t) tanh(pi/2* sinh(b*t));       % DE forward map

    % Left and right FIO matrices
    Il = frac_coeffs(n, mu, b, -1);
    Ir = frac_coeffs(n, mu, b, 1);

    % Scale to Riesz form
    Il = Il./(2*cos(pi*mu/2));
    Ir = Ir./(2*cos(pi*mu/2));

    % Exact solution: (1+x)^{1/2} scaled appropriately
    x = phi(chebpts(n));
    xx = chebpts(n);
    y = pi * sinh(b * xx);
    ln_ff = zeros(size(y));
    idx_pos = (y >= 0);
    ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
    idx_neg = ~idx_pos;
    ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));
    yy = 2*gamma(1/2)*cos(pi/4)*exp(1/2 * ln_ff);
    solc = chebtech2.vals2coeffs(yy);

    % Right-hand side f (computed stably) plus the exact solution term
    f = compute_fx_stable(xx, b);
    fc = chebtech2.vals2coeffs(f);
    fc = fc + solc;

    I = eye(n) + Il + Ir;
    sol = I \ fc;

    err = [err norm(solc - sol, 'inf')];
end

% ---------- Plot ----------
figure('Units', 'inches', 'Position', [1 1 8 4]);
semilogy(N, err, 'LineWidth', 2);
ax = gca; ax.FontSize = 12;
xlabel('$N$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('error', 'FontSize', 16, 'Interpreter', 'latex');

% ---------- Stable computation of the forcing function f ----------
function f_val = compute_fx_stable(y, omega)
    % Evaluate f(y) using numerically stable formulas that avoid overflow
    % in sinh / cosh for large arguments.
    v = (pi/2) * sinh(omega .* y);
    f_val = zeros(size(v));

    idx_pos = (v >= 0);
    if any(idx_pos)
        vp = v(idx_pos);
        E = exp(-vp);  E2 = E.^2;
        term1 = pi ./ (1 + E2);
        term2 = (2 .* E) ./ sqrt(1 + E2);
        term3 = 2 .* asinh(E) ./ (1 + E2);
        f_val(idx_pos) = term1 + term2 + term3;
    end

    idx_neg = (v < 0);
    if any(idx_neg)
        vn = v(idx_neg);
        E = exp(vn);  E2 = E.^2;
        term1 = (pi .* E2) ./ (1 + E2);
        term2 = 2 ./ sqrt(1 + E2);
        asinh_safe = -vn + log(1 + sqrt(1 + E2));   % stable asinh for negative large v
        term3 = (2 .* E2 .* asinh_safe) ./ (1 + E2);
        f_val(idx_neg) = term1 + term2 + term3;
    end
end
