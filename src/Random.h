#ifndef INC_teh_RandomZufall_h
#define INC_teh_RandomZufall_h

/*!
@file		Random.h
@brief 		Random number generator for uniform, normal, and Poisson distributions
		using a lagged Fibonacci series generator.
@author		Thomas E. Hilinski <https://github.com/tehilinski>
@details{
	The zufall code is wrapped in a class, and the seed size is
	larger (unsigned int) versus int in the original.

	An accessory class RandomSeed can generate a seed in the
	valid range for zufall.

	Original description of zufall:

	This package contains a portable random number generator set
	for: uniform (u in [0,1)), normal (<g> = 0, <g^2> = 1), and
	Poisson distributions. The basic module, the uniform generator,
	uses a lagged Fibonacci series generator:

		t    = u(n-273) + u(n-607)
		u(n) = t - float(int(t))

	where each number generated, u(k), is floating point. Since
	the numbers are floating point, the left end boundary of the
	range contains zero. This package is portable except that
	the test package contains some machine dependent timing data.
	These are cycle times (in seconds) for NEC SX-3, Fujitsu VP2200,
	Cray Y-MP, and Sun-4. Select your favorite and comment out the
	others. There are also vectorization directives for Cray Y-MP
	machines in the form of "pragma _CRI", which should be ignored
	although perhaps complained about by other compilers. Otherwise
	the package is portable and returns the same set of floating
	point numbers up to word precision on any machine.

	The maximum seed size is max(unsigned int) - 144.

	External documentation, "Lagged Fibonacci Random Number Generators
	for the NEC SX-3," is to be published in the International
	Journal of High Speed Computing (1994). Otherwise, ask the
	author:

		W. P. Petersen
		IPS, RZ F-5
		ETHZ
		CH 8092, Zurich
		Switzerland
		e-mail:  wpp@ips.ethz.ch.
}
@copyright{
    Copyright 2020 Thomas E. Hilinski. All rights reserved.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
    and in the accompanying file LICENSE.md.

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
}
*/

#include <vector>
#include <utility>
#include <climits>

namespace teh {
  namespace zufall {

    class RandomSeed
    {
      public:

	/// Constructor
	/// @param entropy	An entropy value; default is std::time(NULL)
	/// @param minSeed	Minimum allowed seed > zero
	/// @param maxSeed	Maximum allowed seed; < RandomSeed::GetUpperLimit()
	RandomSeed (
	    unsigned int const entropy = 0U,
	    unsigned int const minSeed = 1U,
	    unsigned int const maxSeed = RandomSeed::GetUpperLimit() - 1 );

	/// Get the seed
	unsigned int Get () const
	  { return seed; }

	/// Get the maximum seed value + 1 for the zufall algorithms.
	static unsigned int GetUpperLimit ()
	  { return seedUpperLimit; }

	/// Get the valid seed range.
	std::pair<unsigned int, unsigned int> GetRange ()
	  { return std::make_pair<>( seedMin, seedMax ); }

      private:

	static unsigned int const seedUpperLimit;
	unsigned int const seedMin, seedMax, seed;
    };

    class Random
    {
      public:

	/// Constructor
	/// @param seed 	0 < seed < max(unsigned int) - 144; default is 1802
	Random (
	  unsigned int const seed );

	/// Generate n uniform random values u[0], ..., u[n-1].
	///
	/// A portable lagged Fibonacci series uniform random number
	/// generator with "lags" -273 und -607:
	/// W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
	///
	/// @param n	number of values to generated
	/// @param u	vector to hold generated random values
	void Uniform (
	    unsigned int const n,
	    std::vector<double> & u );

	/// Get one uniform random value.
	double Uniform ();

	/// Box-Muller method to generate Gaussian random numbers.
	///
	/// Generate n normally distributed values g[0], ..., g[n-1] such that
	/// mean = g = 0, and variance = g**2 = 1.
	///
	/// @param n	number of values to generated
	/// @param g	vector to hold generated random values
	void Normal (
	    unsigned int const n,
	    std::vector<double> & g );

	/// Get one normally distributed value.
	double Normal ();

	/// Generate Poisson-distributed values.
	///
	/// Generate n random integers q from a Poisson density distribution,
	/// with density p(q, mu) = exp(-mu) (mu**q) / q!
	///
	/// @param n	number of values to generated
	/// @param mu	positive real number == expected value of q
	/// @param q	vector to hold generated random values
	void Poisson (
	    unsigned int const n,
	    double const mu,
	    std::vector<unsigned int> & q );

	/// Get one Poisson-distributed value.
	/// @param mu	positive real number == expected value of q
	unsigned int Poisson (
	    double const mu );

      private:

	void InitializeSeedBuffer (unsigned int const seed);

	// Saves seed buffer in saveBuffer[608] for later restarts
	// using function zufallrs.
	void zufallsv (
	  std::vector<double> & buffer);	//  size = klotz0_1_.size + 1

	// Restores seed buffer in saveBuffer[608] which was previously
	// saved using function zufallsv.
	void zufallrs (
	  std::vector<double> const & buffer);	//  size = klotz0_1_.size + 1

	// Saves seed buffer in saveBuffer[1634] for later restarts
	// using function normalrs.
	void normalsv (
	  std::vector<double> & buffer);	//  size = 1634

	// Restores seed buffer in saveBuffer[1634] which was previously
	// saved using function normalsv.
	void normalrs (
	  std::vector<double> const & buffer);	//  size = 1634

	void normal00 ();

    };

  } // namespace zufall
} // namespace teh

#endif // INC_teh_RandomZufall_h
