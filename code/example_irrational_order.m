%% example_irrational_order.m
%  Solves a fractional integral equation with irrational-order fractional
%  operators: (I + I^{e/3} + I^{pi/4}) u = 1.
%  Demonstrates that the DE-based method handles irrational orders naturally.

b = 3.65;
N = 5:2:200;
I11 = frac_coeffs(N(end), exp(1)/3, b, -1);
I22 = frac_coeffs(N(end), pi/4, b, -1);

u = zeros(N(end), length(N));
for r = 1:length(N)
    n = N(r);

    I1 = I11(1:n, 1:n);
    I2 = I22(1:n, 1:n);

    I_mat = eye(n) + I1 + I2;
    f = @(x) 1;
    fc = doubexp(f, n, b, [-1, 1]);
    utemp = I_mat \ fc;
    u(1:n, r) = utemp;
end

ud = diff(u')';
err = zeros(length(N)-1, 1);
for k = 1:length(N)-1
    err(k) = norm(ud(:, k));
end

figure('Units', 'inches', 'Position', [1 1 6 5]);
semilogy(N(1:end-1), err + eps, 'LineWidth', 1.5);
xlabel('N', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('Cauchy error', 'FontSize', 16, 'Interpreter', 'latex');
