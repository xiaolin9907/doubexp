%% example_approx.m
%  Demonstrates spectral convergence of the double-exponential (DE) transform
%  for approximating functions with endpoint singularities.
%  Compares naive DE approximation vs. the improved formula (3.7) from the paper.

mu = 1/3;
f = @(x) (1-x).^(1/3).*(1+x).^(1/2);
b = max(3.15, asinh((mu*log(2) - log(eps))/(mu*pi)));

errd = [];  % naive DE error
errdm = []; % improved formula (3.7) error
dom = [-1,1];

xcheck = linspace(dom(1), dom(2), 20000);
fx = f(xcheck);
tinv = @(y) asinh(2/pi .* atanh(y))./b;
xx = tinv(xcheck);
xx(xx > 1) = 1;
xx(xx < -1) = -1;

n = 1:1:500;

for k = n
    fc = doubexp(f, k, b, dom);
    fcx = clenshaw_chebyshev(fc, 2/(dom(2)-dom(1)).*(xx - (dom(1)+dom(2))/2));
    errd = [errd norm(fx - fcx, 'Inf')];
end

for k = n
    x = chebpts(k);
    y = pi * sinh(b * x);
    ln_ff = zeros(size(y));
    idx_pos = (y >= 0);
    ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
    idx_neg = ~idx_pos;
    ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));
    yy1 = exp(mu * ln_ff);

    x = chebpts(k);
    y = pi * sinh(b * x);
    ln_ff = zeros(size(y));
    idx_pos = (y >= 0);
    ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
    idx_neg = ~idx_pos;
    ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));
    yy2 = exp(1/2 * ln_ff);
    fc = chebtech2.vals2coeffs(flip(yy1).*yy2);

    fcx = clenshaw_chebyshev(fc, 2/(dom(2)-dom(1)).*(xx - (dom(1)+dom(2))/2));
    errdm = [errdm norm(fx - fcx, 'Inf')];
end

figure('Units', 'inches', 'Position', [1 1 8 4]);
semilogy(n, errd, 'LineWidth', 2); hold on;
semilogy(n, errdm, 'LineWidth', 2);
ax = gca; ax.FontSize = 12;
legend({'naive', '(3.7)'}, 'NumColumns', 1, 'FontSize', 12, 'Interpreter', 'latex', 'Location', 'best');
xlabel('$N$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('error', 'FontSize', 16, 'Interpreter', 'latex');

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
