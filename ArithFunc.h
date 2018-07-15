#pragma once

#include <cmath>
#include <vector>
#include <functional>
#include <complex>

namespace AF 
{
typedef std::vector<uint> kTuple;
typedef std::complex<double> cpx;
typedef std::function<cpx(uint)> primAF; 

//The summation limit for non-graphical series computation
uint compSeriesLength;

struct ArithFunc
	{
	//Returns the Dirichlet Identity ArithFunc
	ArithFunc();
	ArithFunc( primAF _func );

	//Component-wise negation
	ArithFunc operator-() const;

	//Component-wise sum
	ArithFunc operator+( ArithFunc g ) const;

	//Component-wise difference
	ArithFunc operator-( ArithFunc g ) const;

	//Dirichlet Product
	ArithFunc operator*( ArithFunc g ) const;

	//Dirichlet Inverse
	ArithFunc operator~( ) const;

	//Function Evaluation
	cpx operator()( uint n ) const;

	//Evalutate at nth prime, only works for primes under 1000 (168 of them). Go set up primegen if you want more.
	cpx atNthPrime( uint n ) const;

	//Evaluate L-Series at s. This is the basic summation form so don't expect reasonable data from divergent points.
	cpx LSeries( cpx s ) const;

	//Render function up to length as a sequence of complex, double, or float
	//Important! The indices get borked by 1, func(1) == func.render(30)[0]
	//Don't blame me, blame programmers and mathematicians for using different conventions
	std::vector<cpx> render(uint length) const;//todo: add time as secondary param
	std::vector<float> renderf(uint length) const;

	primAF func;
	}; 

///ArithFunc generators ------------------------------------------------------------------
//Note that these are not Arthfuncs -- they are functions returning ArithFuncs

//Takes a k-tuple, returns the phi_k ArithFunc for that tuple
ArithFunc phi_k( kTuple k );

//Returns an Arithfunc which returns its input raised to the kth power
ArithFunc n_k( int k );

//Returns an ArithFunc which takes n and returns the sum of the positive divisors of n, each raised to the kth power
ArithFunc sigma_k( int k );

//The k-unit function of the provided k-tuple
ArithFunc u_k( kTuple k );

///ArithFuncs ----------------------------------------------------------------------------

//Constantly returns 1
const ArithFunc u( [](uint){ return 1.0; } );

//The Dirichlet inverse of u
const ArithFunc mu( ~u );

//The standard relative primality counting phi
const ArithFunc phi(phi_k( { 0 } ));

//The identity function which returns its input (not the Dirichlet Identity)
const ArithFunc n( [](uint n ){ return n; } );

//Sum of positive divisors
const ArithFunc sigma( sigma_k(1) );

//Number of positive divisors
const ArithFunc d( u*u );

//It's the Mangoldt function
const ArithFunc mangoldt( ArithFunc(  []( uint n ){ return log( n ); }  )*mu );

}
