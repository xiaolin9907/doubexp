%% example_basis_plot.m
%  Visualizes the transplanted Chebyshev polynomials (TCPs) Q_n(x) under
%  the double-exponential (DE) and algebraic (JFP) variable transforms.
%  TCPs are defined as  Q_n(x) = T_n(psi^{-1}(x))  where psi maps [-1,1]
%  onto itself.

b = 3.15;

% Standard Chebyshev polynomials (as chebfun objects)
T0 = chebpoly(0);
T1 = chebpoly(1);
T2 = chebpoly(2);
T3 = chebpoly(3);
T4 = chebpoly(4);

x = chebpts(1000);

% ========  DE-transformed Chebyshev polynomials  ========
tinv = @(y) asinh(2/pi .* atanh(y))./b;
xx_de = tinv(x);
xx_de(xx_de > 1) = 1;                 % clamp to domain
xx_de(xx_de < -1) = -1;

y0 = T0(xx_de); y1 = T1(xx_de); y2 = T2(xx_de);
y3 = T3(xx_de); y4 = T4(xx_de);
Y_de = [y0 y1 y2 y3 y4];

figure('Units', 'inches', 'Position', [1 1 6*4/3 5*2/3]);
plot(x, y0, 'LineWidth', 2); hold on;
for k = 2:5, plot(x, Y_de(:, k), 'LineWidth', 2); end
xlim([-1.1 1.1]); ylim([-1.1 1.1]);
legend({'Q_0(x)', 'Q_1(x)', 'Q_2(x)', 'Q_3(x)', 'Q_4(x)'}, ...
    'Position', [0.36, 0.24, 0.2, 0.15], 'NumColumns', 2, 'FontSize', 12);
title('DE-transformed Chebyshev polynomials (TCPs)');

% ========  Algebraic (JFP) transformed Chebyshev polynomials  ========
xx_jfp = 2*sqrt((x+1)/2) - 1;         % algebraic map
y0_j = T0(xx_jfp); y1_j = T1(xx_jfp); y2_j = T2(xx_jfp);
y3_j = T3(xx_jfp); y4_j = T4(xx_jfp);
Y_jfp = [y0_j y1_j y2_j y3_j y4_j];

figure('Units', 'inches', 'Position', [1 1 6*4/3 5*2/3]);
plot(x, y0_j, 'LineWidth', 2); hold on;
for k = 2:5, plot(x, Y_jfp(:, k), 'LineWidth', 2); end
xlim([-1.1 1.1]); ylim([-1.1 1.1]);
legend({'Q_0(x)', 'Q_1(x)', 'Q_2(x)', 'Q_3(x)', 'Q_4(x)'}, ...
    'Position', [0.792, 0.29, 0.1, 0.1], 'NumColumns', 1, 'FontSize', 12);
title('Algebraic (JFP) transformed Chebyshev polynomials');
