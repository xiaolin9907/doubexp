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
├── README.md
└── code/
    ├── frac_coeffs.m               # Main FIO matrix construction (low-rank + ultraspherical)
    ├── chebvals2chebcoeffs.m       # Chebyshev value-to-coefficient conversion
    ├── example_approx.m            # DE convergence for singular functions
    ├── example_select_beta.m       # Optimal parameter β selection
    ├── example_left_fie_eq.m       # Left-sided fractional integral equation
    ├── example_irrational_order.m  # Irrational-order FIE
    ├── example_mixed_order.m       # Mixed-order FIE
    ├── example_riesz.m             # Riesz fractional integral
    ├── example_eigenvalue.m        # Eigenvalue analysis
    ├── example_airy.m              # Fractional Airy equation
    └── example_basis_plot.m        # TCP basis function plots
```

## Quick Start

```matlab
addpath('code');

% Construct the fractional integral operator matrix
mu = 1/2;      % order of fractional integral
b = 3.9;       % DE transform parameter
N = 100;       % number of Chebyshev points

% Left-sided FIO
I_left = frac_coeffs(N, mu, b, -1);

% Right-sided FIO
I_right = frac_coeffs(N, mu, b, 1);

% Solve a fractional integral equation: u + I^{1/2}[u] = f
I = frac_coeffs(N, mu, b, -1);
f_coeffs = randn(N, 1);  % right-hand side coefficients
u = (eye(N) + I) \ f_coeffs;
```

## Key Algorithm

### `frac_coeffs.m` — Main algorithm

Constructs the spectral approximation matrix for the left- or right-sided fractional integral operator based on the double-exponential transform. The algorithm:

1. Computes a **low-rank SVD approximation** of the kernel $G(y,t)$
2. Solves a **three-term recurrence ODE** for the moments using the ultraspherical spectral method
3. Assembles the final operator matrix

The matrix is **lower-banded** with bandwidth $K+1$ where $K$ is the degree of the TCP approximation of $(1\pm\psi(y))^\mu$.

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
