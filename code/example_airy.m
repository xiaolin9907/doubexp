%% example_airy.m
%  Solves a fractional Airy-type equation with the DE spectral method.
%  The equation involves a fractional derivative of order 3/2.
%  The solution is obtained by solving a bordered linear system that
%  enforces the boundary condition.

ep = 1e-4;              % small perturbation parameter
n = 2000;
mu = 3/2;
b = max(3.15, asinh((mu*log(2) - log(eps))/(mu*pi)));
fprintf('n: %d\n', n);

% ---------- Construct the FIO matrix ----------
tic;
IN = frac_coeffs(n, 3/2, b, -1);
toc;

% ---------- Precompute multiplication operators on the DE grid ----------
tic;
xx = chebpts(150);
phixx = tanh(pi/2 .* sinh(b .* xx));      % physical coordinate x on the DE grid
y = pi * sinh(b * xx);

% Numerically stable computation of singular factors
ln_ff = zeros(size(y));
idx_pos = (y >= 0);
ln_ff(idx_pos) = log(2) - log1p(exp(-y(idx_pos)));
idx_neg = ~idx_pos;
ln_ff(idx_neg) = log(2) - (-y(idx_neg) + log1p(exp(y(idx_neg))));

% (1+x)^{1/2}
yy12 = exp(1/2 * ln_ff);
fc12 = chebtech2.vals2coeffs(yy12);

% x * (1+x)^{1/2}
yyx12 = phixx .* exp(1/2 * ln_ff);
fcx12 = chebtech2.vals2coeffs(yyx12);

% x * (1+x)^{3/2} + boundary correction term
yyx32 = phixx .* exp(3/2 * ln_ff) + ep * (1i)^(3/2) / gamma(1/2);
fcg = chebtech2.vals2coeffs(yyx32);

% Multiplication operators
M12 = ultraS.multmat(n, fc12, 0);
Mx12 = ultraS.multmat(n, fcx12, 0);
toc;

% ---------- Assemble the operator ----------
op = ep*(1i)^(3/2) * M12 - Mx12 * IN;

% ---------- Bordered system for the boundary condition ----------
B = ones(1, n);                          % boundary condition row

OP = [B*IN, 2;                           % [  B*IN  |  2  ]
      op, [fcg; zeros(n - length(fcg), 1)]];  % [  op    | fcg ]
l = zeros(n+1, 1); l(1) = 1;             % RHS

tic;
u = OP \ l; a = u(end);                  % a = Lagrange multiplier
toc;
u = IN * u(1:end-1);                     % solution coefficients

ucheb = chebfun(u, 'coeffs');            % represent as chebfun

% ---------- Plot solution ----------
xxx = chebpts(30000); uval = ucheb(xxx);
phixxx = tanh(pi/2 .* sinh(b .* xxx));
plot(phixxx, real(uval), 'LineWidth', 1.5); hold on;
plot(phixxx, imag(uval), 'LineWidth', 1.5);
legend('real', 'imag');
xlabel('x', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('u(x)', 'FontSize', 16, 'Interpreter', 'latex');

% Uncomment to also plot coefficient decay:
% semilogy(1:length(u), abs(u), 'LineWidth', 1.5);
% ax = gca; ax.FontSize = 12;
% xlabel('n', 'FontSize', 16, 'Interpreter', 'latex');
% ylabel('coefficients', 'FontSize', 16, 'Interpreter', 'latex');
