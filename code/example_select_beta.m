%% example_select_beta.m
%  Demonstrates how to select the optimal parameter beta (omega) for the
%  double-exponential transform. The optimal beta depends on the singularity
%  strength mu. The dashed lines mark the theoretically optimal values
%   beta_opt(mu) = max(3.15, asinh(...)).

% ---------- Parameter sweep ----------
MU = [1 1/10 1/100 1e-3 1e-4 1e-5 1e-6];
beta_opt = [];
BB = max(3.15, asinh((log(2) - log(eps)/MU(end) + log1p(-0.5 * eps^(1/MU(end)))) / pi));
ERR = [];

for j = 1:length(MU)
    mu = MU(j);
    f = @(x) (1+x).^mu;                % test function with endpoint singularity
    B = 2:0.05:BB+8;                   % range of beta values to try
    errdm = [];
    dom = [-1,1];
    xcheck = linspace(dom(1), dom(2), 100);
    fx = f(xcheck);
    n = 700;

    for i = 1:length(B)
        b = B(i);
        % DE-transformed grid points
        xx = chebpts(n);
        y = pi * sinh(b * xx);
        % Numerically stable computation of (1+x)^mu on the DE grid
        ln_ff = zeros(size(y));
        idx_pos = (y >= 0);
        ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
        idx_neg = ~idx_pos;
        ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));
        yy = exp(mu * ln_ff);
        fc = chebtech2.vals2coeffs(yy);

        % Map physical points back through the DE transform and evaluate
        tinv = @(y) asinh(2/pi .* atanh(y))./b;
        xx = tinv(xcheck);
        xx(xx > 1) = 1;
        xx(xx < -1) = -1;

        fcx = clenshaw_chebyshev(fc, 2/(dom(2)-dom(1)).*(xx - (dom(1)+dom(2))/2));
        errdm = [errdm norm(fx - fcx, 'Inf')];
    end
    ERR = [ERR; errdm];
    % Theoretical optimal beta for this mu
    beta_opt = [beta_opt max(3.154, asinh((log(2) - log(eps)/mu + log1p(-0.5 * eps^(1/mu))) / pi))];
end

% ---------- Plotting ----------
col = [
    0.10, 0.35, 0.65;  % Navy blue
    0.85, 0.20, 0.25;  % Brick red
    0.15, 0.60, 0.35;  % Forest green
    0.95, 0.55, 0.15;  % Orange
    0.50, 0.30, 0.65;  % Purple
    0.20, 0.75, 0.85;  % Teal
    0.90, 0.40, 0.60;  % Pink
];

MU_str = {'1', '10^{-1}', '10^{-2}', '10^{-3}', '10^{-4}', '10^{-5}', '10^{-6}'};

figure('Units', 'inches', 'Position', [1 1 8 4]);
for i = 1:7
    leg_str = sprintf('$\\mu = %s$', MU_str{i});
    semilogy(B, ERR(i,:), '-', 'Color', col(i,:), 'LineWidth', 2, 'DisplayName', leg_str);
    if i == 1, hold on; end
    semilogy([beta_opt(i) beta_opt(i)], [1e-16 1], '--', 'Color', col(i,:), 'LineWidth', 2, 'HandleVisibility', 'off');
end

legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
ylim([1e-17, 1]); xlim([2 BB+6]);
ax = gca; ax.FontSize = 12; box on;
xlabel('$\omega$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('error', 'FontSize', 16, 'Interpreter', 'latex');
yticks([1e-15 1e-10 1e-5 1]);

% ---------- Clenshaw's algorithm ----------
function y = clenshaw_chebyshev(c, x)
    N = length(c) - 1;
    y = zeros(size(x));
    for k = 1:length(x)
        b0 = 0; b1 = 0;
        for n = N:-1:1
            b2 = b1; b1 = b0;
            b0 = c(n+1) + 2 * x(k) * b1 - b2;
        end
        y(k) = c(1) + x(k) * b0 - b1;
    end
end
