# Double-Exponential Transform for Fractional Calculus

This repository contains MATLAB code accompanying the paper:

> **Fractional calculus via variable-transform-based spectral approximations**
> by Xiaolin Liu and Kuan Xu

## Overview

We present a unifying framework for constructing spectral approximations to fractional integral operators (FIOs) using **transplanted Chebyshev polynomials (TCPs)**. The key idea is to compose Chebyshev polynomials with a variable transform that maps $[-1,1]$ onto itself:

$$Q_n(x) = T_n(\psi^{-1}(x))$$

When an exponential (or double-exponential) transform is used for $\psi$, the framework yields a versatile spectral approximation applicable to a much broader class of fractional calculus problems than existing methods — including equations with irrational-order fractional operators and mixed derivative types.

## Dependencies

The code requires **[Chebfun](https://www.chebfun.org/)** for MATLAB.

## Repository Structure

```
.
├── src/                          # Core algorithm implementations
│   ├── doubexp.m                 # DE transform approximation of singular functions
│   ├── frac_coeffs.m             # Main FIO matrix construction (low-rank + ultraspherical)
│   ├── frac_DE.m                 # Direct FIO matrix via DE quadrature
│   ├── two_order_diff.m          # Differentiation matrices under DE transform
│   ├── chebvals2chebcoeffs.m     # Chebyshev value-to-coefficient conversion
│   └── IMT_map.m                 # IMT transform (for comparison)
│
├── examples/                     # Demonstration scripts
│   ├── example_approx.m          # DE convergence for singular functions
│   ├── example_select_beta.m     # Optimal parameter β selection
│   ├── example_DE_vs_SE.m        # DE vs single-exponential comparison
│   ├── example_ode_second_order.m    # Integer-order ODE with DE method
│   ├── example_ode_variable.m        # Variable-coefficient ODE
│   ├── example_left_fie.m            # Left-sided FIO accuracy test
│   ├── example_left_fie_eq.m         # Left-sided fractional integral equation
│   ├── example_irrational_order.m    # Irrational-order FIE
│   ├── example_mixed_order.m         # Mixed-order FIE
│   ├── example_riesz.m               # Riesz fractional integral
│   ├── example_eigenvalue.m          # Eigenvalue analysis
│   ├── example_airy.m                # Fractional Airy equation
│   ├── example_pseudospectra.m       # Pseudospectra of fractional operators
│   ├── example_kernel_test.m         # Kernel approximation quality
│   ├── example_matrix_viz.m          # Matrix sparsity visualization
│   ├── example_basis_plot.m          # TCP basis function plots
│   ├── example_weight.m              # Weight function ω(t)
│   ├── example_mlf.m                 # Mittag-Leffler function plot
│   └── example_kernel_plot.m         # Kernel G(y,t) analysis
│
└── paper/                        # Paper (for reference)
    ├── de.pdf                    # Preprint PDF
    └── de.tex                    # LaTeX source
```

## Quick Start

### 1. Function approximation with the DE transform

```matlab
addpath('src');
mu = 1/3;
f = @(x) (1-x).^(1/3) .* (1+x).^(1/2);
b = max(3.15, asinh((mu*log(2) - log(eps)) / (mu*pi)));

% Approximate f(x) using N DE-transformed Chebyshev coefficients
N = 200;
fc = doubexp(f, N, b, [-1, 1]);
```

### 2. Constructing the fractional integral operator matrix

```matlab
mu = 1/2;
b = 3.9;
N = 100;

% Left-sided FIO of order mu
I_left = frac_coeffs(N, mu, b, -1);

% Right-sided FIO of order mu
I_right = frac_coeffs(N, mu, b, 1);
```

### 3. Solving a fractional integral equation

```matlab
% Solve: u + I^{1/2}[u] = f
I = frac_coeffs(N, mu, b, -1);
u = (eye(N) + I) \ f_coeffs;
```

## Key Algorithms

### `frac_coeffs.m` — Main algorithm

Constructs the spectral approximation matrix for the left- or right-sided FIO based on the DE transform. The algorithm:

1. Computes a **low-rank SVD approximation** of the kernel $G(y,t)$
2. Solves a **three-term recurrence ODE** for the moments $\varphi_n^j(y)$ using the ultraspherical spectral method
3. Assembles the final operator matrix

The matrix is **lower-banded** with bandwidth $K+1$ where $K$ is the degree of the TCP approximation of $(1\pm\psi(y))^\mu$.

### `doubexp.m` — Function approximation

Approximates a function $f(x)$ with endpoint singularities using the DE-transformed Chebyshev expansion. The DE transform clusters points exponentially near the endpoints, resolving singularities at spectral convergence rates.

### `two_order_diff.m` — Differentiation matrices

Constructs the first and second-order differentiation matrices under the DE variable transform, needed for solving differential equations.

## License

This code is provided for research purposes. Please cite the accompanying paper if you use it in your work.

## Reference

```bibtex
@article{liu2026spectral,
  title   = {{Spectral approximation to fractional integral operators}},
  author  = {Liu, Xiaolin and Xu, Kuan},
  journal = {Mathematics of Computation},
  year    = {2026},
  doi     = {10.1090/mcom/4185},
}
```
