%% example_left_fie_eq.m
%  Solves the fractional integral equation:  u(x) + I^mu[u](x) = f(x)
%  for various mu.  The exact solution is (1+x)^mu.
%  Demonstrates spectral convergence as N increases.

% ---------- Problem parameters ----------
Mu = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1]; Mu = flip(Mu);
Gn = [920 660 350 170 120 100 80]; Gn = flip(Gn);   % NN for kernel resolution
N = 10:2:300;
err = zeros(length(Mu), length(N));

for ii = 1:length(Mu)
    mu = Mu(ii);
    b = max(3.15, asinh((mu*log(2) - log(eps))/(mu*pi)));   % optimal DE parameter
    IN = frac_coeffs(N(end), mu, b, -1, Gn(ii));

    for jj = 1:length(N)
        n = N(jj);

        % Right-hand side: f = (1+x)^mu + Gamma(1+mu)/Gamma(1+2mu) * (1+x)^{2mu}
        xx = chebpts(n);
        y = pi * sinh(b * xx);
        ln_ff = zeros(size(y));
        idx_pos = (y >= 0);
        ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
        idx_neg = ~idx_pos;
        ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));

        % First term: (1+x)^mu
        yy = exp(mu * ln_ff);
        fc = chebtech2.vals2coeffs(yy);

        % Second term: Gamma(1+mu)/Gamma(1+2mu) * (1+x)^{2mu}
        yy2 = gamma(1+mu)/gamma(1+2*mu) * exp(2*mu * ln_ff);
        fc2 = chebtech2.vals2coeffs(yy2);

        f = fc + fc2;
        In = IN(1:n, 1:n);

        % Solve and measure error against exact solution fc
        u = (eye(n) + In) \ f;
        err(ii, jj) = norm(u - fc, 'inf');
    end
end

% ---------- Plot ----------
figure('Units', 'inches', 'Position', [1 1 8 4]);
semilogy(N, err, 'LineWidth', 2);
legend({'$\mu =~~1$', '$\mu = 10^{-1}$', '$\mu = 10^{-2}$', '$\mu = 10^{-3}$', ...
    '$\mu = 10^{-4}$', '$\mu = 10^{-5}$', '$\mu = 10^{-6}$'}, ...
    'Location', 'northeast', 'NumColumns', 1, 'FontSize', 12, 'Interpreter', 'latex');
xlabel('$N$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('error', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca; ax.FontSize = 12;
ylim([1e-16, 1]);
