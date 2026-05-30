%% example_airy.m
%  Solves a fractional Airy-type equation with the DE spectral method.
%  Equation involves a fractional derivative of order 3/2.
%  The solution is obtained by solving a bordered linear system.

ep = 1e-4;
n = 2000;
mu = 3/2;
b = max(3.15, asinh((mu*log(2) - log(eps))/(mu*pi)));
fprintf('n: %d\n', n);

tic;
IN = frac_coeffs(n, 3/2, b, -1);
toc;

tic;
xx = chebpts(150);
phixx = tanh(pi/2 .* sinh(b .* xx));
y = pi * sinh(b * xx);
ln_ff = zeros(size(y));
idx_pos = (y >= 0);
ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
idx_neg = ~idx_pos;
ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));

yy12 = exp(1/2 * ln_ff);
fc12 = chebtech2.vals2coeffs(yy12);

yyx12 = phixx .* exp(1/2 * ln_ff);
fcx12 = chebtech2.vals2coeffs(yyx12);

yyx32 = phixx .* exp(3/2 * ln_ff) + ep * (1i)^(3/2) / gamma(1/2);
fcg = chebtech2.vals2coeffs(yyx32);

M12 = ultraS.multmat(n, fc12, 0);
Mx12 = ultraS.multmat(n, fcx12, 0);
toc;

op = ep*(1i)^(3/2) * M12 - Mx12 * IN;

B = ones(1, n);

OP = [B*IN, 2; op, [fcg; zeros(n - length(fcg), 1)]];
l = zeros(n+1, 1); l(1) = 1;

tic;
u = OP \ l; a = u(end);
toc;
u = IN * u(1:end-1);

ucheb = chebfun(u, 'coeffs');

% Uncomment to plot solution:
xxx = chebpts(30000); uval = ucheb(xxx);
phixxx = tanh(pi/2 .* sinh(b .* xxx));
plot(phixxx, real(uval), 'LineWidth', 1.5); hold on;
plot(phixxx, imag(uval), 'LineWidth', 1.5);
legend('real', 'imag');
xlabel('x', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('u(x)', 'FontSize', 16, 'Interpreter', 'latex');

% Plot coefficient decay
% semilogy(1:length(u), abs(u), 'LineWidth', 1.5);
% ax = gca; ax.FontSize = 12;
% xlabel('n', 'FontSize', 16, 'Interpreter', 'latex');
% ylabel('coefficients', 'FontSize', 16, 'Interpreter', 'latex');
