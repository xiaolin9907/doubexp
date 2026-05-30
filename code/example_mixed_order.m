%% example_mixed_order.m
%  Solves a FIE with multiple fractional operators of different orders:
%   u + a(x) I^{sqrt(2)}[u] + b(x) I^{e/3}[u] + c(x) I^{pi/4}[u] = 1.
%  Demonstrates the DE method's ability to handle mixed irrational orders
%  with variable-coefficient multipliers.

b = 3.65;
N = 5:2:300;
mu = sqrt(2);

% Precompute operators at a larger N (extra padding for banded structure)
Ill = frac_coeffs(N(end)+10, mu, b, -1);
Irr = frac_coeffs(N(end)+10, mu, b, 1);
I11 = frac_coeffs(N(end)+10, exp(1)/3, b, 1);
I22 = frac_coeffs(N(end)+10, pi/4, b, -1);

u = zeros(N(end), length(N));
for r = 1:length(N)
    n = N(r);

    % Map Chebyshev points through the DE transform to physical space
    x = chebpts(n);
    phi = @(t) tanh(pi/2* sinh(b*t));
    x = phi(x);

    % Variable coefficient a(x) — a smooth function of x
    ff = @(x) 2*exp(pi/2*exp(b*x)) ./ (exp(pi/2*exp(b*x)) + exp(pi/2*exp(-b*x)));
    xx = chebpts(n);
    yy = ff(xx).^(2/3);
    p1c = chebtech2.vals2coeffs(yy);

    % Variable coefficient c(x) = (1-x)^{sqrt(3)} — endpoint singularity
    p2 = @(x) (1-x).^(sqrt(3));
    p2c = doubexp(p2, n, b, [-1, 1]);

    % Multiplication operators
    M1 = ultraS.multmat(n, p1c, 0);
    M2 = ultraS.multmat(n, p2c, 0);

    % Slice operators to current N
    I1 = I11(1:n, 1:n);
    I2 = I22(1:n, 1:n);
    Il = Ill(1:n, 1:n);
    Ir = Irr(1:n, 1:n);

    % Riesz combination: (I_l + I_r) / (2 cos(pi*mu/2))
    Il = Il./(2*cos(pi*mu/2));
    Ir = Ir./(2*cos(pi*mu/2));
    IR = Il + Ir;

    % Assemble: I + M1*IR + I1*M2 + I2
    I_mat = eye(n) + M1*IR + I1*M2 + I2;

    f = @(x) 1;
    fc = doubexp(f, n, b, [-1, 1]);
    utemp = I_mat \ fc;
    u(1:n, r) = utemp;
end

% Cauchy error
ud = diff(u')';
err = zeros(length(N)-1, 1);
for k = 1:length(N)-1
    err(k) = norm(ud(:, k));
end

% ---------- Plot Cauchy error ----------
figure('Units', 'inches', 'Position', [1 1 8 4]);
semilogy(N(1:end-1), err, 'LineWidth', 2);
xlabel('$N$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('Cauchy error', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca; ax.FontSize = 12;

% ---------- Plot final solution ----------
figure('Units', 'inches', 'Position', [1 1 8 4]);
x = chebpts(N(end));
tinv = @(y) asinh(2/pi .* atanh(y))./b;
ys = clenshaw_chebyshev(u(:, end), tinv(x));
plot(x, ys, 'LineWidth', 2);
ax = gca; ax.FontSize = 12;
xlabel('$x$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$u(x)$', 'FontSize', 16, 'Interpreter', 'latex');

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
