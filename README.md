Arbitrary-precision C++ implementation for evaluating hydrogenic radial probability integrals using GMP.
# Hydrogenic Radial Probability Evaluation Using GMP

## Overview

This repository contains a C++ implementation using the GNU Multiple Precision Arithmetic Library (GMP) to evaluate hydrogenic radial probability integrals using arbitrary-precision arithmetic.

The program is intended for numerical investigations of highly excited hydrogenic states, where cancellation in finite-sum expressions can cause severe numerical errors in standard double-precision calculations.

The arbitrary-precision implementation evaluates the analytical expression at increased numerical precision. It is intended to support numerical-stability investigations and reproducibility; it does not constitute an independent proof of the analytical formula.

## Files

- `compute_prob_GMP_v5.cpp` — C++ source code.
- `README_github_compute_prob_GMP_v7.txt` — detailed instructions and usage information.

## Requirements

- A C++14-compatible compiler.
- Boost C++ Libraries, including Boost.Multiprecision.
- The GNU Multiple Precision Arithmetic Library (GMP), including its development headers and libraries.

See `README_github_compute_prob_GMP_v7.txt` for installation guidance, compilation instructions, command-line syntax, input conventions, examples, and numerical limitations.

## Licence

The original code in this repository is made available under the Creative Commons Attribution 4.0 International licence (CC BY 4.0), subject to the terms of that licence. The GMP library is a separate third-party dependency and is governed by its own licence.

Licence: https://creativecommons.org/licenses/by/4.0/

## Citation

If you use this code in research, please cite the repository and identify the version used. If a separate publication or archived release is provided later, include that citation as well.
