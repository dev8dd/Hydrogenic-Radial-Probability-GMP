================================================================
HYDROGENIC PROBABILITY FRAMEWORK - BOOST GMP C++ IMPLEMENTATION
================================================================
Author: Devdatta Mankame
 
Associated Manuscript: "A Unified Analytical Framework for Hydrogenic
Probability Distributions: From Orthogonal Polynomials to 3D Nodal
Topology"


1. OVERVIEW
----------------------------------------------------------------
This program computes hydrogenic radial probabilities over arbitrary
finite or semi-infinite radial intervals using the closed-form
finite-sum expression developed in the associated manuscript.

For an interval specified by the dimensionless radial boundaries
m1 and m2, the program evaluates the radial probability

    P(m1 <= r <= m2)
        = integral from m1 to m2 of |R_nl(r)|^2 r^2 dr,

where the input boundaries m1 and m2 are expressed in units of a0, so that r1 = m1 a0 and r2 = m2 a0. The nuclear charge Z enters the radial wavefunction and the probability calculation, not the input-unit convention.

The program supports:

    - finite intervals [m1, m2],
    - lower-tail probabilities [0, m2],
    - upper-tail probabilities [m1, infinity),
    - the full normalization interval [0, infinity).

The implementation uses Boost.Multiprecision with the GNU Multiple
Precision (GMP) arithmetic backend:

    boost::multiprecision::gmp_float<256>

providing 256 decimal digits of working precision.

The higher-precision arithmetic is intended to reduce numerical errors
associated with large factorial/gamma structures, large intermediate
terms, and cancellation in the finite-sum evaluation of the analytical
probability expression.

The analytical expression remains exact independently of the numerical
precision used for its evaluation. The numerical value returned by
the program, however, is subject to the working precision and the
conditioning of the finite-sum calculation.

No universal principal-quantum-number threshold is assumed for the
256-digit implementation. The precision required for reliable
numerical evaluation is state- and interval-dependent and should be
assessed by numerical convergence when demanding high-n calculations
are performed.


2. PREREQUISITES & DEPENDENCIES
----------------------------------------------------------------
The program requires:

    - A C++14-compliant compiler.
    - Boost C++ Libraries with Boost.Multiprecision.
    - GNU Multiple Precision Arithmetic Library (GMP).
    - GMP development headers and libraries.

The program uses the Boost.Multiprecision GMP backend through:

    <boost/multiprecision/gmp.hpp>

The main standard-library headers used by the source include:

    <iostream>
    <vector>
    <cmath>
    <algorithm>
    <string>
    <iomanip>

Installation examples:

Linux (Ubuntu/Debian):
    sudo apt-get install libgmp-dev libboost-dev

macOS (Homebrew):
    brew install gmp boost

Windows (MSYS2 / MinGW-w64):
    pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-boost

The exact package names may vary between operating systems and
distribution versions.


3. COMPILATION INSTRUCTIONS
----------------------------------------------------------------
The program requires C++14 support and must be linked against the
GMP library.

Linux / macOS:
    g++ -O3 -std=c++14 compute_prob_GMP_v5.cpp -o compute_prob_GMP_v5 -lgmp

Windows (MSYS2 / MinGW):
    g++ -O3 -std=c++14 compute_prob_GMP_v5.cpp -o compute_prob_GMP_v5.exe -lgmp

The -O3 compiler option enables optimization and can improve execution
time for high-precision evaluations.

The Boost.Multiprecision GMP backend is header-based from the Boost side,
while the numerical GMP library itself is linked through:

    -lgmp


4. USAGE & SYNTAX
----------------------------------------------------------------
The executable accepts five positional arguments:

    1. n   = principal quantum number
    2. l   = orbital angular momentum quantum number
    3. Z   = nuclear charge or effective nuclear charge
    4. m1  = lower radial boundary in units of a0
    5. m2  = upper radial boundary in units of a0

Syntax:

    ./compute_prob_GMP_v5 <n> <l> <Z> <m1> <m2>


5. RADIAL BOUNDARY CONVENTION
----------------------------------------------------------------
The program uses the following convention for the upper boundary:

    m2 = -1.0

to represent

    m2 -> infinity.

Therefore:

    m1 = 0, m2 = -1.0

represents the complete radial interval

    [0, infinity).

The program evaluates this case directly as the full normalization
and returns

    P = 1.

For a lower-tail probability,

    0 <= m1 < m2,

the program evaluates the cumulative probability up to m2 and subtracts
the cumulative probability up to m1.

For an upper-tail probability,

    m2 = -1.0,

the program evaluates

    P(r >= m1) = 1 - P(0 <= r <= m1).

For a finite interval,

    0 <= m1 < m2,

the program evaluates the difference between the cumulative
probabilities at the two boundaries.


6. EXAMPLES
----------------------------------------------------------------
Example 1:
Probability between r = 0 and r = 50 a0 for n = 100, l = 1, Z = 1.

    ./compute_prob_GMP_v5 100 1 1.0 0 50


Example 2:
Upper-tail probability for r >= 15000 a0 for n = 100, l = 1,
Z = 1.

    ./compute_prob_GMP_v5 100 1 1.0 15000 -1.0


Example 3:
Probability in the finite interval [1000, 5000] a0 for an effective
charge Z_eff = 1.6875.

    ./compute_prob_GMP_v5 100 1 1.6875 1000 5000


Example 4:
Full normalization check for n = 2, l = 0, Z = 1.

    ./compute_prob_GMP_v5 2 0 1.0 0 -1.0

The expected result is numerically consistent with

    P = 1.


Sample output format:

    P = 0.00005101121113773009332300999818784042679651849466
    for radial interval [0, 50] in units of a0 (m2 = -1 denotes infinity);
    n = 100, l = 1, Z = 1.0


7. PRECISION ARCHITECTURE
----------------------------------------------------------------
The numerical type used throughout the main probability calculation
is defined as:

    using mp_float =
        boost::multiprecision::number<
            boost::multiprecision::gmp_float<256>
        >;

This provides 256 decimal digits of working precision.

The implementation uses logarithmic factorial/gamma quantities when
constructing the finite-sum terms. In particular:

    log_fact()

computes logarithmic factorial values, while

    mp_lgamma_int()

provides the corresponding logarithmic gamma quantity for integer
arguments.

This avoids directly constructing extremely large factorial values
during the intermediate calculations.

The Laguerre polynomial is evaluated through a finite alternating sum.
The individual terms are constructed using logarithmic factorial
quantities before exponentiation.

This approach substantially extends the usable numerical range compared
with ordinary double precision, but it does not eliminate cancellation
in an alternating finite sum. Consequently, numerical reliability must
still be assessed for demanding high-n states or extreme radial
intervals.


8. PRECISION AND NUMERICAL LIMITATIONS
----------------------------------------------------------------
The 256-digit implementation is designed to provide a high-precision
numerical realization of the analytical probability formula.

However, arbitrary-precision arithmetic should not be interpreted as
an automatic guarantee of numerical accuracy for every possible
quantum number or radial interval.

The numerical result can be affected by:

    - cancellation between large alternating terms,
    - the magnitude of the radial argument,
    - the degree of the associated Laguerre polynomial,
    - the selected radial interval,
    - the dynamic range of intermediate quantities.

For demanding calculations, numerical convergence should be checked
by repeating the calculation with increased working precision and
comparing the resulting probability.

For example, the source can be modified from

    gmp_float<256>

to a higher-precision configuration when additional precision is
required, subject to the capabilities of the selected Boost/GMP
configuration.


9. INPUT CONSTRAINTS & VALIDATION
----------------------------------------------------------------
The physical quantum-number constraints are:

    n >= 1

    0 <= l < n

    Z > 0

    m1 >= 0

For the upper boundary:

    m2 > m1

or

    m2 = -1.0

to denote infinity.

The program checks the integer quantum-number conditions and verifies
the orbital constraint l < n before evaluating the probability.

The nuclear/effective charge Z may be an integer or a positive
fractional value.

Examples include:

    Z = 1.0
    Z = 1.6875


10. OUTPUT FORMAT
----------------------------------------------------------------
The probability is printed using fixed-point notation with

    50 decimal places.

The output also reports the numerical radial interval in units of a0.
The upper-boundary sentinel m2 = -1 is displayed as supplied and denotes
infinity; it is not a physical radial coordinate.

Example:

    P = 0.XXXXXXXX...
    for radial interval [m1, m2] in units of a0 (m2 = -1 denotes infinity);
    n = n, l = l, Z = Z

The displayed decimal digits should not automatically be interpreted
as significant digits. The number of reliable digits depends on the
working precision and the numerical conditioning of the calculation.


11. SPECIAL NORMALIZATION CHECK
----------------------------------------------------------------
The command

    ./compute_prob_GMP_v5 n l Z 0 -1.0

represents the complete radial interval

    0 <= r < infinity.

For a correctly normalized hydrogenic radial wavefunction, the
analytical result is

    P = 1.

This provides a useful numerical normalization check for the
implementation.

The program explicitly returns 1 for this special input rather than
evaluating a difference of two cumulative expressions.


12. TROUBLESHOOTING
----------------------------------------------------------------
Problem:
    fatal error: boost/multiprecision/gmp.hpp: No such file or directory

Solution:
    Install the Boost development libraries.

    Ubuntu/Debian:
        sudo apt-get install libboost-dev

    macOS:
        brew install boost

    MSYS2:
        pacman -S mingw-w64-x86_64-boost


Problem:
    fatal error: gmp.h: No such file or directory

Solution:
    Install the GMP development package.

    Ubuntu/Debian:
        sudo apt-get install libgmp-dev

    macOS:
        brew install gmp

    MSYS2:
        pacman -S mingw-w64-x86_64-gmp


Problem:
    undefined reference to `__gmpf_...`

Solution:
    Ensure that the GMP library is linked at the end of the compilation
    command:

        -lgmp


Problem:
    The probability becomes numerically unstable for a high-n state.

Possible causes include:

    - insufficient working precision,
    - cancellation among large alternating finite-sum terms,
    - large radial arguments,
    - demanding radial intervals.

Possible remedies:

    - increase the working precision,
    - repeat the calculation at a higher precision,
    - compare the numerical results at two or more precisions,
    - independently verify the result using numerical quadrature where
      practical.

Do not assume that increasing precision alone guarantees a correct
result. Convergence of the numerical value should be checked.


Problem:
    The result is close to, but not exactly, 1 for a full normalization
    check.

Use:

    ./compute_prob_GMP_v5 n l Z 0 -1.0

The program recognizes this interval as the full radial domain and
returns the normalization value directly.


13. RELATION TO THE STANDARD DOUBLE-PRECISION IMPLEMENTATION
----------------------------------------------------------------
This program provides a higher-precision implementation of the same
analytical radial probability framework used by the corresponding
standard C++ probability program.

The principal numerical difference is the arithmetic precision:

    Standard implementation:
        IEEE-754 double precision

    This implementation:
        Boost.Multiprecision with GMP backend,
        256 decimal digits by default

The analytical finite-sum expression being evaluated is the same
underlying formulation. The GMP implementation is intended for
calculations where standard double precision may be insufficient due
to cancellation, large intermediate quantities, or the numerical
demands of highly excited states.


14. REPRODUCIBILITY
----------------------------------------------------------------
For reproducible calculations, record at minimum:

    - n
    - l
    - Z
    - m1
    - m2
    - numerical precision used
    - compiler and compiler version
    - Boost version
    - GMP version

For demanding high-n calculations, numerical convergence should also
be documented by comparing results obtained at different working
precisions.

When possible, an independent numerical quadrature calculation can
provide an additional check on the analytical finite-sum evaluation.


15. CITATION & LICENSE
----------------------------------------------------------------
If you use this probability framework or its implementation in
research, please cite the associated manuscript:

Devdatta  Mankame,
"A Unified Analytical Framework for Hydrogenic Probability
Distributions: From Orthogonal Polynomials to 3D Nodal Topology,"
2026.

License:

    Creative Commons Attribution 4.0 International (CC BY 4.0).

================================================================
END OF README
================================================================