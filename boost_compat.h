#ifndef BOOST_COMPAT_H_
#define BOOST_COMPAT_H_

/**
 * Boost-compatible RNG implementations (no Boost dependency).
 * Matches Boost.Random's uniform_01 and gamma_distribution algorithms
 * so that results are identical to the original Boost-based RSEM.
 */

#include <cmath>
#include <limits>

// Boost's uniform_01 formula: (eng() - min) / (max - min + 1) for [0, 1)
// Matches boost/random/uniform_01.hpp new_uniform_01
template<typename Engine>
inline double boost_uniform_01(Engine& eng) {
	typedef typename Engine::result_type base_result;
	const double factor = 1.0 / (static_cast<double>(eng.max() - eng.min()) +
		(std::numeric_limits<base_result>::is_integer ? 1.0 : 0.0));
	double result = static_cast<double>(eng() - eng.min()) * factor;
	// For integer engines, max output gives result < 1, so no retry needed
	return result;
}

// Boost's gamma_distribution (Knuth algorithm) - matches boost/random/gamma_distribution.hpp
template<typename Engine>
inline double boost_gamma(Engine& eng, double alpha, double beta) {
	using std::tan; using std::sqrt; using std::exp; using std::log; using std::pow;

	if (alpha == 1.0) {
		double u = boost_uniform_01(eng);
		if (u >= 1.0) u = 1.0 - 1e-10;
		return (-std::log(1.0 - u)) * beta;
	} else if (alpha > 1.0) {
		const double pi = 3.14159265358979323846;
		for (;;) {
			double y = std::tan(pi * boost_uniform_01(eng));
			double x = sqrt(2.0 * alpha - 1.0) * y + alpha - 1.0;
			if (x <= 0.0) continue;
			if (boost_uniform_01(eng) > (1.0 + y * y) * exp((alpha - 1.0) * log(x / (alpha - 1.0)) - sqrt(2.0 * alpha - 1.0) * y))
				continue;
			return x * beta;
		}
	} else {
		// alpha < 1: use gamma(alpha+1) * u^(1/alpha) transformation
		double p = std::exp(1.0) / (alpha + std::exp(1.0));
		for (;;) {
			double u = boost_uniform_01(eng);
			double y;  // exponential(1)
			{
				double eu = boost_uniform_01(eng);
				if (eu >= 1.0) eu = 1.0 - 1e-10;
				y = -std::log(1.0 - eu);
			}
			double x, q;
			if (u < p) {
				x = exp(-y / alpha);
				q = p * exp(-x);
			} else {
				x = 1.0 + y;
				q = p + (1.0 - p) * pow(x, alpha - 1.0);
			}
			if (u >= q) continue;
			return x * beta;
		}
	}
}

#endif /* BOOST_COMPAT_H_ */
