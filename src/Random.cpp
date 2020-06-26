//------------------------------------------------------------------------------------------------------------
//	File:	  Random.cpp
//	Class:	  teh::Random
//
//	Description:
//	Random number generator using the code from zufall.h/.c.
//	See header file for more details.
//
//	Author:	Thomas E. Hilinski <https://github.com/tehilinski>
//
//	Copyright 2020 Thomas E. Hilinski. All rights reserved.
//	Licensed under the Apache License, Version 2.0 (the "License");
//	you may not use this file except in compliance with the License.
//	You may obtain a copy of the License at
//	    http://www.apache.org/licenses/LICENSE-2.0
//	    and in the accompanying file LICENSE.md.
//	Unless required by applicable law or agreed to in writing, software
//	distributed under the License is distributed on an "AS IS" BASIS,
//	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//	See the License for the specific language governing permissions and
//	limitations under the License.
}
//------------------------------------------------------------------------------------------------------------

#include "Random.h"
#include <cmath>
using std::exp;
using std::cos;
using std::sin;
using std::sqrt;
using std::log;
#include <stdexcept>
#include <string>
#include <ctime>

namespace teh {
  namespace zufall {

    //---------------------------------------- class RandomSeed ----------------------------------------

    unsigned int const RandomSeed::seedUpperLimit = UINT_MAX - 144U;

    RandomSeed::RandomSeed (
	unsigned int const entropy,
	unsigned int const minSeed,
	unsigned int const maxSeed )
	: seedMin ( minSeed > 0 ? minSeed : 1 ),
	  seedMax ( maxSeed < seedUpperLimit ? maxSeed : seedUpperLimit - 1 ),
	  seed    ( entropy == 0 ?
		    std::time(NULL) % GetRange().second :
		    entropy % GetRange().second )
    {
	if ( seed < seedMin || seed >= seedMax )
	{
	    std::string msg = "Random: seed is out of range.";
	    throw std::runtime_error( msg );
	}
	if ( seedMin >= seedMax )
	{
	    std::string msg = "Either minimum or maximum seed (or both) is invalid.";
	    throw std::runtime_error( msg );
	}
    }

    //---------------------------------------- private to Random ----------------------------------------

    struct klotz0_1_
    {
	static short const size = 607;
	std::vector<double> buff;
	unsigned int ptr;

	klotz0_1_ ()
	  : buff (size, 0.0),
	    ptr (0)
	  {
	  }
    };

    klotz0_1_ klotz0_1;

    // private to Random
    struct klotz1_1_
    {
	static short const size = 1024;
	std::vector<double> xbuff;
	unsigned int first, xptr;

	klotz1_1_ ()
	  : xbuff (size, 0.0),
	    first (0),
	    xptr (0)
	  {
	  }
    };

    static klotz1_1_ klotz1_1;

    //---------------------------------------- class Random ----------------------------------------

    Random::Random (
	unsigned int const seed )
    {
	InitializeSeedBuffer (seed);
    }

    // ----------------------------------------------------------------------------
    //	InitializeSeedBuffer (formerly zufalli)
    //	Generates initial seed buffer by linear congruential
    //	method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
    //	Variable seed should be 0 < seed < ( numeric_limits<unsigned int> - 144 )
    // ----------------------------------------------------------------------------
    void Random::InitializeSeedBuffer (
	unsigned int const seed)
    {
	unsigned int const ij = seed;
	unsigned int i = ij / 177 % 177 + 2;
	unsigned int j = ij % 177 + 2;
	unsigned int k = 9373 / 169 % 178 + 1;
	unsigned int l = 9373 % 169;
	for (unsigned int ii = 0; ii < (unsigned int)klotz0_1.size; ++ii)
	{
	    double s = 0.0;
	    double t = 0.5;
	    for (unsigned int jj = 1; jj <= 24; ++jj)
	    {
		unsigned int m = i * j % 179 * k % 179;
		i = j;
		j = k;
		k = m;
		l = (l * 53 + 1) % 169;
		if (l * m % 64 >= 32)
		    s += t;
		t *= 0.5;
	    }
	    klotz0_1.buff[ii] = s;
	}
    }

    // ----------------------------------------------------------------------------
    //	Uniform (formerly zufall)
    //	portable lagged Fibonacci series uniform random number
    //	generator with "lags" -273 und -607:
    //	W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
    // ----------------------------------------------------------------------------
    double Random::Uniform ()
    {
	std::vector<double> a;
	Uniform( 1, a );
	return a[0];
    }

    void Random::Uniform (
	unsigned int const n,
	std::vector<double> & a)
    {
	if (n == 0)
	    return;

	a.assign (n, 0.0);

	unsigned int left, bptr, aptr0;
	double t;
	unsigned int vl, k273, k607, kptr;
	unsigned int aptr = 0;
	unsigned int nn = n;

	L1:

	if (nn == 0)
	    return;

	/* factor nn = q*607 + r */

	unsigned int q = (nn - 1) / klotz0_1.size;
	left = klotz0_1.size - klotz0_1.ptr;

	if (q <= 1)
	{

	    /* only one or fewer full segments */

	    if (nn < left)
	    {
		kptr = klotz0_1.ptr;
		for (unsigned int i = 0; i < nn; ++i)
		    a[i + aptr] = klotz0_1.buff[kptr + i];
		klotz0_1.ptr += nn;
		return;
	    }
	    else
	    {
		kptr = klotz0_1.ptr;
		#pragma _CRI ivdep
		for (unsigned int i = 0; i < left; ++i)
		    a[i + aptr] = klotz0_1.buff[kptr + i];
		klotz0_1.ptr = 0;
		aptr += left;
		nn -= left;
		/*  buff -> buff case */
		vl = 273;
		k273 = 334;
		k607 = 0;
		for (unsigned int k = 0; k < 3; ++k)
		{
		    #pragma _CRI ivdep
		    for (unsigned int i = 0; i < vl; ++i)
		    {
		       t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
		       klotz0_1.buff[k607+i] = t - (double) ((unsigned int) t);
		    }
		    k607 += vl;
		    k273 += vl;
		    vl = 167;
		    if (k == 0)
		    {
			k273 = 0;
		    }
		}
		goto L1;
	    }
	}
	else
	{

	    /* more than 1 full segment */

	    kptr = klotz0_1.ptr;
	    #pragma _CRI ivdep
	    for (unsigned int i = 0; i < left; ++i)
		a[i + aptr] = klotz0_1.buff[kptr + i];
	    nn -= left;
	    klotz0_1.ptr = 0;
	    aptr += left;

	    /* buff -> a(aptr0) */

	    vl = 273;
	    k273 = 334;
	    k607 = 0;
	    for (unsigned int k = 0; k < 3; ++k)
	    {
		if (k == 0)
		{
		    #pragma _CRI ivdep
		    for (unsigned int i = 0; i < vl; ++i)
		    {
			t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
			a[aptr + i] = t - (double) ((unsigned int) t);
		    }
		    k273 = aptr;
		    k607 += vl;
		    aptr += vl;
		    vl = 167;
		}
		else
		{
		    #pragma _CRI ivdep
		    for (unsigned int i = 0; i < vl; ++i)
		    {
			t = a[k273 + i] + klotz0_1.buff[k607 + i];
			a[aptr + i] = t - (double) ((unsigned int) t);
		    }
		    k607 += vl;
		    k273 += vl;
		    aptr += vl;
		}
	    }
	    nn += -klotz0_1.size;

	    /* a(aptr-607) -> a(aptr) for last of the q-1 segments */

	    aptr0 = aptr - klotz0_1.size;
	    vl = klotz0_1.size;

	    for (unsigned int qq = 0; qq < q-2; ++qq)
	    {
		k273 = aptr0 + 334;
		#pragma _CRI ivdep
		for (unsigned int i = 0; i < vl; ++i)
		{
		    t = a[k273 + i] + a[aptr0 + i];
		    a[aptr + i] = t - (double) ((unsigned int) t);
		}
		nn += -klotz0_1.size;
		aptr += vl;
		aptr0 += vl;
	    }

	    /* a(aptr0) -> buff, last segment before residual */

	    vl = 273;
	    k273 = aptr0 + 334;
	    k607 = aptr0;
	    bptr = 0;
	    for (unsigned int k = 0; k < 3; ++k)
	    {
		if (k == 0)
		{
		    #pragma _CRI ivdep
		    for (unsigned int i = 0; i < vl; ++i)
		    {
			t = a[k273 + i] + a[k607 + i];
			klotz0_1.buff[bptr + i] = t - (double) ((unsigned int) t);
		    }
		    k273 = 0;
		    k607 += vl;
		    bptr += vl;
		    vl = 167;
		}
		else
		{
		    #pragma _CRI ivdep
		    for (unsigned int i = 0; i < vl; ++i)
		    {
			t = klotz0_1.buff[k273 + i] + a[k607 + i];
			klotz0_1.buff[bptr + i] = t - (double) ((unsigned int) t);
		    }
		    k607 += vl;
		    k273 += vl;
		    bptr += vl;
		}
	    }
	    goto L1;
	}
    }

    // ----------------------------------------------------------------------------
    //	Normal (formerly normalen)
    //	Box-Muller method for Gaussian random numbers.
    // ----------------------------------------------------------------------------
    double Random::Normal ()
    {
	std::vector<double> x;
	Uniform( 1, x );
	return x[0];
    }

    void Random::Normal (
	unsigned int const n,
	std::vector<double> & x)
    {
	if (n == 0)
	    return;

	x.assign (n, 0.0);

	if (klotz1_1.first == 0)
	{
	    normal00();
	    klotz1_1.first = 1;
	}
	unsigned int ptr = 0;
	unsigned int nn = n;

      L1:
	unsigned int const left = klotz1_1.size - klotz1_1.xptr;
	unsigned int const kptr = klotz1_1.xptr;
	if (nn < left)
	{
	    #pragma _CRI ivdep
	    for (unsigned int i = 0; i < nn; ++i)
		x[i + ptr] = klotz1_1.xbuff[kptr + i];
	    klotz1_1.xptr += nn;
	    return;
	}
	else
	{
	    #pragma _CRI ivdep
	    for (unsigned int i = 0; i < left; ++i)
		x[i + ptr] = klotz1_1.xbuff[kptr + i];
	    klotz1_1.xptr = 0;
	    ptr += left;
	    nn -= left;
	    normal00();
	    goto L1;
	}
    }

    // ----------------------------------------------------------------------------
    //	Poisson (formerly fische)
    //	Poisson generator for distribution function of p's:
    //	q(mu,p) = exp(-mu) mu**p/p!
    // ----------------------------------------------------------------------------
    unsigned int Random::Poisson (
	double const mu)
    {
	std::vector<unsigned int> p;
	Poisson( 1, mu, p );
	return p[0];
    }

    void Random::Poisson (
	unsigned int const n,
	double const mu,
	std::vector<unsigned int> & p)
    {
	if (n == 0)
	    return;

	p.assign (n, 0);

	double const pmu = exp(-mu);
	unsigned int p0 = 0;
	unsigned int nsegs = (n - 1) / klotz1_1.size;
	unsigned int left = n - (nsegs << 10);
	++nsegs;
	unsigned int nl0 = left;

	std::vector<unsigned int> indx (klotz1_1.size);
	std::vector<double> q (klotz1_1.size);
	std::vector<double> u (klotz1_1.size);

	for (unsigned int k = 0; k < nsegs; ++k)
	{
	    for (unsigned int i = 0; i < left; ++i)
	    {
		indx[i] = i;
		p[p0 + i] = 0;
		q[i] = 1.;
	    }

	    /* Begin iterative loop on segment of p's */
	    do
	    {
		Uniform (left, u);	// Get the needed uniforms
		unsigned int jj = 0;
		for (unsigned int i = 0; i < left; ++i)
		{
		    unsigned int const ii = indx[i];
		    double const q0 = q[ii] * u[i];
		    q[ii] = q0;
		    if (q0 > pmu)
		    {
			indx[jj++] = ii;
			++p[p0 + ii];
		    }
		}
		left = jj;		// any left in this segment?
	    }
	    while (left > 0);

	    p0  += nl0;
	    nl0  = klotz1_1.size;
	    left = klotz1_1.size;
	}
    }

    // ----------------------------------------------------------------------------
    //	zufallsv
    //	saves common blocks klotz0, containing seeds and
    //	pointer to position in seed block. IMPORTANT: svblk must be
    //	dimensioned at least 608 in driver. The entire contents
    //	of klotz0 (pointer in buff, and buff) must be saved.
    // ----------------------------------------------------------------------------
    void Random::zufallsv (
	std::vector<double> & saveBuffer)			//  size = klotz0_1_.size + 1
    {
	saveBuffer.resize (klotz0_1_::size + 1);
	saveBuffer[0] = (double) klotz0_1.ptr;
	#pragma _CRI ivdep
	for (short i = 0; i < klotz0_1.size; ++i)
	    saveBuffer[i + 1] = klotz0_1.buff[i];
    }

    // ----------------------------------------------------------------------------
    //	zufallrs
    //	restores common block klotz0, containing seeds and pointer
    //	to position in seed block. IMPORTANT: saveBuffer must be
    //	dimensioned at least 608 in driver. The entire contents
    //	of klotz0 must be restored.
    // ----------------------------------------------------------------------------
    void Random::zufallrs (
	std::vector<double> const & saveBuffer)			//  size = klotz0_1_.size + 1
    {
	if ( saveBuffer.size() != klotz0_1_::size + 1 )
	{
	    std::string msg;
	    msg = "Random::zufallrs, restore of uninitialized block.";
	    throw std::runtime_error (msg);
	}
	klotz0_1.ptr = (unsigned int) saveBuffer[0];
	#pragma _CRI ivdep
	for (short i = 0; i < klotz0_1.size; ++i)
	    klotz0_1.buff[i] = saveBuffer[i + 1];
    }

    // ----------------------------------------------------------------------------
    //	normalsv
    //	save zufall block klotz0
    //	IMPORTANT: svbox must be dimensioned at
    //	least 1634 in driver. The entire contents of blocks
    //	klotz0 (via zufallsv) and klotz1 must be saved.
    // ----------------------------------------------------------------------------
    void Random::normalsv (
	    std::vector<double> & saveBuffer)				// size = 1634
    {
	if (klotz1_1.first == 0)
	{
	    std::string msg;
	    msg = "Random::normalsv, save of uninitialized block.";
	    throw std::runtime_error (msg);
	}

	saveBuffer.resize (1634);
	zufallsv (saveBuffer);

	saveBuffer[klotz0_1.size + 1] = (double) klotz1_1.first;	// [608]
	saveBuffer[klotz0_1.size + 2] = (double) klotz1_1.xptr;		// [609]
	unsigned int const k = klotz0_1.size + 3;			// 610
	#pragma _CRI ivdep
	for (short i = 0; i < klotz1_1.size; ++i)
	    saveBuffer[i + k] = klotz1_1.xbuff[i];
    }

    // ----------------------------------------------------------------------------
    //	normalrs
    //	restore zufall blocks klotz0 and klotz1
    //	IMPORTANT: saveBuffer must be dimensioned at
    //	least 1634 in driver. The entire contents
    //	of klotz0 and klotz1 must be restored.
    // ----------------------------------------------------------------------------
    void Random::normalrs (
	std::vector<double> const & saveBuffer)				// size = 1634
    {
	zufallrs (saveBuffer);

	klotz1_1.first = (unsigned int) saveBuffer[klotz0_1.size + 1];	// [608]
	if (klotz1_1.first == 0)
	{
	    std::string msg;
	    msg = "Random::normalrs, restore of uninitialized block.";
	    throw std::runtime_error (msg);
	}

	klotz1_1.xptr = (unsigned int) saveBuffer[klotz0_1.size + 2];	// [609]
	unsigned int const k = klotz0_1.size + 3;			// 610
	#pragma _CRI ivdep
	for (short i = 0; i < klotz1_1.size; ++i)
	    klotz1_1.xbuff[i] = saveBuffer[i + k];
    }

    // ----------------------------------------------------------------------------
    //	normal00
    // ----------------------------------------------------------------------------
    void Random::normal00 ()
    {
	/* Builtin functions */
	/* double cos(), sin(), log(), sqrt(); */

	double const twopi = 6.2831853071795862;

	Uniform (klotz1_1.size, klotz1_1.xbuff);

// 	std::vector<double>::iterator i   = klotz1_1.xbuff.begin();
// 	std::vector<double>::iterator ip1 = klotz1_1.xbuff.begin();
// 	while ( ip1 != klotz1_1.xbuff.end() )
// 	{
// 	    double const r1 = twopi * (*i);
// 	    double const t1 = cos(r1);
// 	    double const t2 = sin(r1);
// 	    double const r2 = sqrt(-2.0 * ( log(1.0 - (*ip1) ) ) );
// 	    (*i)   = t1 * r2;
// 	    (*ip1) = t2 * r2;
// 	    ++i;   ++i;
// 	    ++ip1; ++ip1;
// 	}

	#pragma _CRI ivdep
	for (short i = 0; i < klotz1_1.size - 1; i += 2)
	{
	    double const r1 = twopi * klotz1_1.xbuff[i];
	    double const t1 = cos(r1);
	    double const t2 = sin(r1);
	    double const r2 = sqrt(-2.0 * ( log(1.0 - klotz1_1.xbuff[i+1] ) ) );
	    klotz1_1.xbuff[i]   = t1 * r2;
	    klotz1_1.xbuff[i+1] = t2 * r2;
	}

    } // class Random

  } // namespace zufall
} // namespace teh
