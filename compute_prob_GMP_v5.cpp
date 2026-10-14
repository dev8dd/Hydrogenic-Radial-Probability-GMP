//================================================================
//HYDROGENIC PROBABILITY FRAMEWORK - C++ IMPLEMENTATION (GMP)
//================================================================
// Author: Devdatta Mankame

//Associated Manuscript: "A Unified Analytical Framework for Hydrogenic Probability Distributions: From Orthogonal Polynomials to 3D Nodal Topology" 
//

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <cctype>
#include <iomanip>
#include <boost/multiprecision/gmp.hpp>

// Define arbitrary precision float (256 decimal digits) mapping directly to GMP backend
using mp_float = boost::multiprecision::number<boost::multiprecision::gmp_float<256>>;

// Compute log(n!) at 256-digit precision with cache
mp_float log_fact(int n) {
    if (n <= 1) return mp_float(0.0);
    static std::vector<mp_float> cache = {mp_float(0.0), mp_float(0.0)};
    if (n < static_cast<int>(cache.size())) return cache[n];
    for (size_t i = cache.size(); i <= static_cast<size_t>(n); ++i) {
        cache.push_back(cache.back() + log(mp_float(i)));
    }
    return cache[n];
}

// Compute log(Gamma(N)) = log((N-1)!) for integer N >= 1
mp_float mp_lgamma_int(int N) {
    if (N <= 1) return mp_float(0.0);
    return log_fact(N - 1);
}

bool is_strict_int(const std::string& s, bool allow_zero = false) {
    if (s.empty()) return false;
    for (char const &c : s) if (!std::isdigit(c)) return false;
    try {
        long long val = std::stoll(s);
        return allow_zero ? (val >= 0) : (val > 0);
    } catch (...) { return false; }
}

bool is_positive_double(const std::string& s) {
    if (s.empty()) return false;
    try {
        double val = std::stod(s);
        return val > 0.0;
    } catch (...) { return false; }
}

mp_float laguerre(int q, int s, mp_float x) {
    if (q < 0) return mp_float(0.0);
    if (q == 0) return mp_float(1.0);
    mp_float sum = 0.0;
    for (int k = 0; k <= q; k++) {
        mp_float log_f = mp_lgamma_int(q + s + 1) 
                       - mp_lgamma_int(q - k + 1) 
                       - mp_lgamma_int(s + k + 1) 
                       - mp_lgamma_int(k + 1);
        mp_float max_val = x > mp_float("1e-18") ? x : mp_float("1e-18");
        mp_float term = exp(log_f + k * log(max_val));
        if (k % 2 == 0) sum += term; else sum -= term;
    }
    return sum;
}

mp_float sum_k_term(int N, mp_float alpha, mp_float m) {
    mp_float sum = 0.0;
    mp_float x = alpha * m;
    for (int k = 0; k <= N; k++) {
        mp_float log_term = -k * log(x) - mp_lgamma_int(N - k + 1);
        sum += exp(log_term);
    }
    return sum;
}

// FIX: Changed 'int Z' to 'mp_float Z' to prevent double-to-int truncation
mp_float evaluate_F(int n, int l, mp_float Z, mp_float m, bool include_residuals) {
    int q = n - l - 1; 
    int s = 2 * l + 1;
    mp_float alpha = (mp_float(2.0) * Z) / mp_float(n);
    mp_float x = alpha * m;
    mp_float sum_i = 0.0;

    for (int i = 0; i <= q - 1; i++) {
        mp_float m_term = (m > mp_float("1e-15")) ? (mp_float(n * i) / (mp_float(2.0) * Z * m)) : mp_float(0.0);
        mp_float bracket = laguerre(q - i, s + i, x) - m_term * laguerre(q - i + 1, s + i - 1, x);  
        sum_i += pow(x, i) * exp(mp_lgamma_int(q - i)) * laguerre(q - i - 1, s + i + 1, x) * bracket;
    }

    mp_float term1 = exp(-x) * pow(m, s + 2) * sum_i;
    mp_float I3_part = exp(-x) * pow(m, s + q + 1) * sum_k_term(s + q + 1, alpha, m);
    if (include_residuals) I3_part -= pow(mp_float(1.0)/alpha, s + q + 1);
    mp_float term_I3 = pow(alpha, q - 1) * mp_float(q + 1) * exp(mp_lgamma_int(s + q + 2)) * I3_part;

    mp_float I4_part = exp(-x) * pow(m, s + q) * sum_k_term(s + q, alpha, m);
    if (include_residuals) I4_part -= pow(mp_float(1.0)/alpha, s + q);
    mp_float term_I4 = pow(alpha, q - 2) * mp_float(q * (q + s)) * exp(mp_lgamma_int(s + q + 1)) * I4_part;

    mp_float prefactor = pow(alpha, s + 2) / (mp_float(2.0 * n) * exp(mp_lgamma_int(q + s + 1)));
    return prefactor * (term1 - term_I3 + term_I4);
}

int main(int argc, char* argv[]) {
    if (argc != 6) { std::cerr << "USAGE: n l Z m1 m2" << std::endl; return 1; }

    if (!is_strict_int(argv[1], false) || !is_strict_int(argv[2], true) || !is_positive_double(argv[3])) {
        std::cerr << "Error: n and l must be integers, and Z must be positive (Z > 0)." << std::endl; 
        return 1;
    }

    int n = std::stoi(argv[1]), l = std::stoi(argv[2]);
    mp_float Z(argv[3]);  // Accepts floating-point Z_eff (e.g., 1.6875 or 1.0)
    mp_float m1(argv[4]); // FIX: Added m1 declaration
    mp_float m2(argv[5]); // FIX: Added m2 declaration

    if (l >= n || m1 < 0 || (m2 != mp_float("-1.0") && m2 <= m1)) {
        std::cerr << "Error: Orbital constraint or interval [m1, m2] invalid." << std::endl; return 1;
    }

    mp_float result = 0.0;
    if (m1 == 0 && m2 == mp_float("-1.0")) result = 1.0;
    else if (m1 == 0) result = evaluate_F(n, l, Z, m2, true);
    else if (m2 == mp_float("-1.0")) result = 1.0 - evaluate_F(n, l, Z, m1, true);
    else result = evaluate_F(n, l, Z, m2, false) - evaluate_F(n, l, Z, m1, false);

// Radial interval inputs m1 and m2 specify the boundaries
// r1 = m1*a0 and r2 = m2*a0.
// Thus, the inputs are in units of a0, not in units
// of the scaled plotting coordinate Z*r/a0.
// The upper-bound sentinel m2 = -1 denotes infinity.

    std::cout << "P = " << std::fixed << std::setprecision(50) << result
              << " for radial interval [" << m1 << ", " << m2
              << "] in units of a0 (m2 = -1 denotes infinity); "
              << "n = " << n << ", l = " << l << ", Z = " << Z
              << std::endl;
    return 0;
}
