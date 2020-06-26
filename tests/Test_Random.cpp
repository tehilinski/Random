//---------------------------------------------------------------------------
// Test of class teh::Random
// Build:
//     g++ -g -I../src -o Test_Random Test_Random.cpp ../src/Random.cpp
//
// Author: Thomas E. Hilinski <https://github.com/tehilinski>
//
// Copyright 2020 Thomas E. Hilinski. All rights reserved.
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//    http://www.apache.org/licenses/LICENSE-2.0
//    and in the accompanying file LICENSE.md.
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//---------------------------------------------------------------------------

#include "Random.h"
#include <ctime>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>

char const * const appName = "Test_Random";
char const * const appTitle = "Test of class teh::Random";

std::size_t const sizeSet = 10000;

char const * const fileNames[] =
    { "Random_Uniform.csv", "Random_Normal.csv", "Random_Poisson.csv" };
enum { IUniform, INormal, IPoisson };

//---------------------------------------------------------------------------
template< class VectorType > void WriteRandomSet (
	std::ofstream & os,
	VectorType const & v)
{
	char const fieldSeparator = ',';
	os << "\"N\",\"Nprev\"" << endl;
	os << v[0] << fieldSeparator << 0.0 << endl;
	for ( std::size_t i = 1; i < v.size(); ++i )
	    os << v[i] << fieldSeparator << v[i - 1] << endl;
}

//---------------------------------------------------------------------------
int main ()
{
	using teh::zufall::RandomSeed;
	using teh::zufall::Random;

	cout << '\n' << appName << ": " << appTitle << endl;

	// seed
	RandomSeed seedGen;
	cout << "Seed: " << seedGen.Get()
	     << " (valid range = "
	     << seedGen.GetRange().first
	     << " - "
	     << seedGen.GetRange().second
	     << ')'
	     << endl;
	Random rng ( seedGen.Get() );

	cout << "Uniform distribution: ";
	std::vector<double> u;
	rng.Uniform (sizeSet, u);
	cout << "size = " << u.size() << ", ";
	std::ofstream osUniform ( fileNames[IUniform] );
	WriteRandomSet (osUniform, u);
	cout << "file = " << fileNames[IUniform] << endl;
	osUniform.close ();
	cout << "One random value = " << rng.Uniform() << endl;

	cout << "Normal distribution: ";
	std::vector<double> g;
	rng.Normal (sizeSet, g);
	cout << "size = " << g.size() << ", ";
	std::ofstream osNormal ( fileNames[INormal] );
	WriteRandomSet (osNormal, g);
	cout << "file = " << fileNames[INormal] << endl;
	osNormal.close ();
	cout << "One random value = " << rng.Normal() << endl;

	/// @todo Poisson distribution
	cout << "Poisson distribution: ";
	std::vector<unsigned int> p;
	double const mu = 2.0;
	rng.Poisson (sizeSet, mu, p);
	cout << "size = " << p.size() << ", ";
	std::ofstream osPoisson ( fileNames[IPoisson] );
	WriteRandomSet (osPoisson, p);
	cout << "file = " << fileNames[IPoisson] << endl;
	osPoisson.close ();
	cout << "One random value = " << rng.Poisson(mu) << endl;

	cout << "\n   all done!" << endl;
	return 0;
}
//---------------------------------------------------------------------------
