%% example_eigenvalue.m
%  Computes eigenvalues of a fractional integral operator using the DE method.
%  Demonstrates Cauchy error convergence of the computed eigenvalues as N increases.

six_colors = [
    1 0 0;       % Red
    0 0.5 0;     % Dark green
    0 0 1;       % Blue
    1 0.5 0;     % Orange
    0.75 0 0.75; % Purple
    0 0.75 0.75  % Cyan
];

mu = 1.23456789;
mu2 = 0.123456789;
b = max(3.15, asinh((mu2*log(2) - log(eps))/(mu2*pi)));

n = 500;

op = frac_coeffs(n, mu, b, -1);
op2 = frac_coeffs(n, mu - mu2, b, -1);
B = ones(1, n);

xx = chebpts(n);
y = pi * sinh(b * xx);
ln_ff = zeros(size(y));
idx_pos = (y >= 0);
ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
idx_neg = ~idx_pos;
ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));
yy = exp((mu-1) * ln_ff);
g = chebtech2.vals2coeffs(yy);

C = (gamma(mu - mu2)/gamma(mu)) * 2^(1 - mu + mu2) * g * B * op2;
OPeig = C - op;

% Compute eigenvalue convergence
ce = [];
l = 20:1:n;
for n_k = l
    [v, d] = eigs(OPeig(1:n_k, 1:n_k), 11);
    d = 1 ./ diag(d);
    d = sort(d);
    ce = [ce d(1:11)];
end

[m1, m2] = size(ce);
ced = zeros(m1, m2-1);
for k = 1:11
    ced(k, :) = diff(ce(k, :)) ./ norm(ce(k, end));
end
ced = abs(ced);
ced(ced == 0) = 1e-16;

ll = length(l) - 1;

figure('Units', 'inches', 'Position', [1 1 8 4]);
semilogy(l(1:ll), ced(1, 1:ll), 'Color', six_colors(1,:), 'LineWidth', 2); hold on;
semilogy(l(1:ll), ced(2, 1:ll), 'Color', six_colors(2,:), 'LineWidth', 2);
semilogy(l(1:ll), ced(3, 1:ll), 'Color', six_colors(3,:), 'LineWidth', 2);
semilogy(l(1:ll), ced(4, 1:ll), 'Color', six_colors(3,:), 'LineWidth', 2);
semilogy(l(1:ll), ced(5, 1:ll), 'Color', six_colors(5,:), 'LineWidth', 2);
semilogy(l(1:ll), ced(6, 1:ll), 'Color', six_colors(4,:), 'LineWidth', 2);
semilogy(l(1:ll), ced(8, 1:ll), 'Color', six_colors(5,:), 'LineWidth', 2);
semilogy(l(1:ll), ced(10, 1:ll), 'Color', six_colors(6,:), 'LineWidth', 2);
xlim([10, n]); ylim([1e-16, 10]);
legend({'\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5', '\lambda_6'}, ...
    'NumColumns', 2, 'FontSize', 12);
xlabel('$N$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Cauchy error', 'FontSize', 12, 'Interpreter', 'latex');
