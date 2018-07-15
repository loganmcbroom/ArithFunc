#include "../include/ArithFunc.h"

using namespace std;
using namespace AF;

uint primes[] =
	{
	2,		3,		5,		7,		11,		13,		17,		19,		23,		29,
	31,		37,		41,		43,		47,		53,		59,		61,		67,		71,
	73,		79,		83,		89,		97,		101,	103,	107,	109,	113,
	127,	131,	137,	139,	149,	151,	157,	163,	167,	173,
	179,	181,	191,	193,	197,	199,	211,	223,	227,	229,
	233,	239,	241,	251,	257,	263,	269,	271,	277,	281,
	283,	293,	307,	311,	313,	317,	331,	337,	347,	349,
	353,	359,	367,	373,	379,	383,	389,	397,	401,	409,
	419,	421,	431,	433,	439,	443,	449,	457,	461,	463,
	467,	479,	487,	491,	499,	503,	509,	521,	523,	541,
	547,	557,	563,	569,	571,	577,	587,	593,	599,	601,
	607,	613,	617,	619,	631,	641,	643,	647,	653,	659,
	661,	673,	677,	683,	691,	701,	709,	719,	727,	733,
	739,	743,	751,	757,	761,	769,	773,	787,	797,	809,
	811,	821,	823,	827,	829,	839,	853,	857,	859,	863,
	877,	881,	883,	887,	907,	911,	919,	929,	937,	941,
	947,	953,	967,	971,	977,	983,	991,	997,	1009
	};
	
uint AF::compSeriesLength = 1024;

uint gcd( uint a, uint b ) 
	{ return b == 0 ? a : gcd( b, a % b ); }

AF::ArithFunc::ArithFunc() : func( []( uint n ){ return n == 1 ? 1 : 0; } ) {}
AF::ArithFunc::ArithFunc( primAF _func ) : func(_func) {}
ArithFunc AF::ArithFunc::operator-() const 
	{ 
	auto f = *this;
	return ArithFunc( [=](uint n)
		{ 
		return -f(n); 
		}); 
	}
ArithFunc AF::ArithFunc::operator+( ArithFunc g ) const
	{ 
	auto f = *this;
	return ArithFunc( [=](uint n)
		{ 
		return f(n) + g(n); 
		}); 
	}
ArithFunc AF::ArithFunc::operator-( ArithFunc g ) const
	{ 
	auto f = *this;
	return ArithFunc( [=](uint n)
		{ 
		return f(n) - g(n); 
		}); 
	}
ArithFunc AF::ArithFunc::operator*( AF::ArithFunc g ) const
	{
	//This allows f and g to be treated symmetrically for clearer code
	auto f = *this;
	return ArithFunc ([=]( uint n )
		{
		cpx sum = 0;

		//We only need to check divisors up to the sqrt of n, one sqrt op is worth the savings on large numbers
		float nf = float( n );
		float sqrt_n = sqrt( nf );

		//Check up to sqrt(n) for divisors
		for( float i = 1; i < sqrt_n; ++i )
			{
			//if i|n, then d=n/i|n. Add both parts of the Dirichlet sum to the total
			float d = nf / i;
			if( d == floor( d ) )
				{
				sum += f( i ) * ( g )( d );
				sum += f( d ) * ( g )( i );
				}
			}
		//If n is a perfect square, we double added the symmetric point in the previous sum, fix that
		if( sqrt_n == floor(sqrt_n) )
			sum += f( sqrt_n ) * g( sqrt_n );

		return sum;
		});
	}
ArithFunc AF::ArithFunc::operator~( ) const
	{
	//Pretty much the same mechanics as the product function here, look there for info
	auto f = *this;
	return ArithFunc ([=]( uint n )
		{
		float nf = float( n );
		if( n == 1 ) return 1.0 / f(1);
		cpx sum = 0;
		float sqrt_n = float(sqrt( nf ));
		//This += and starting the loop at 2 handles the inverse formula summing over d<n
		sum += f(nf) * (~f)( 1 );
		for( float i = 2; i < sqrt_n; ++i ) 
			{
			float d = nf / i;
			if( d == floor( d ) )
				{
				sum += f( i ) * (~f)( d );
				sum += f( d ) * (~f)( i );
				}
			}
		if( sqrt_n == floor( sqrt_n ) )
			sum += f( sqrt_n ) * (~f)( sqrt_n );

		return sum * (-1.0 / f(1));
		});
	}
cpx AF::ArithFunc::operator()( uint n ) const
	{ 
	if( n == 0 ) 
		{
		std::cout << "ArithFuncs start at one. This is math, not programming!" << std::endl;
		}
	return func( n );
	}
cpx AF::ArithFunc::atNthPrime( uint n ) const { return this->operator()(primes[n-1]); }
cpx AF::ArithFunc::LSeries( AF::cpx s ) const
	{
	cpx sum = 0;
	for( int n = 1; n <= AF::compSeriesLength; ++n )
		sum +=  this->operator()(n) * pow(n,-s);
	return sum;
	}
vector<cpx> AF::ArithFunc::render(uint length) const
	{
	vector<cpx> sequence(length);
	for(int i = 0; i< length; ++i)
		sequence[i] = this->operator()(i+1);
	return sequence;
	}
vector<float> AF::ArithFunc::renderf(uint length) const
	{
	vector<float> sequence(length);
	
	for(int i = 0; i < length; ++i)
		{
		sequence[i] = func(i+1).real();
		}
	return sequence;
	}

//This could be optimized by several methods (e.g. symmetry for big numbers)
ArithFunc AF::phi_k( AF::kTuple k )
	{
	return (ArithFunc) [=]( uint n )
		{
		uint total = 0;
		for( int i = 0; i < n; ++i )
			{
			bool coprime = true;
			for( auto j : k )
				{
				if( gcd( ( j + i ) % n, n ) > 1 ) 
					coprime = false;
				}
			if( coprime ) ++total;
			}
		return total;
		};
	}
ArithFunc AF::n_k( int k ) { return (ArithFunc) [=]( uint n ) { return pow(n,k); }; }
ArithFunc AF::sigma_k( int k ) { return (ArithFunc) [=]( uint n ){ return ( AF::n_k(k)*AF::u )(n); }; }
ArithFunc AF::u_k( AF::kTuple k ) { return ~AF::phi_k( k ) * AF::n; }
