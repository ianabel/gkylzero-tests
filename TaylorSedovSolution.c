
#include <math.h>
#include <stdio.h>

#include "TaylorSedov.h"

double NewtonRootSolve( double x0, double ( *f )( double ), double ( *fPrime )( double ), double fTol, double rTol, unsigned int nMaxIters )
{
	double x = x0;
	unsigned int i;
	for ( i = 0; i < nMaxIters; ++i )
	{
		double y = x;
		x = y - f( y )/fPrime( y );
		if ( fabs( f( x ) ) < fTol )
			break;
		if ( fabs( x - y ) / fabs( x ) < rTol )
			break;
	}
	if ( i == nMaxIters ) {
		printf( "Exceeded max iterations in newton solve\n" );
	}
	return x;
}

double gas_gamma;

// Exponents arising in the solution

double Nu1() {
	return -( 13.0*gas_gamma*gas_gamma - 7.0*gas_gamma + 12.0 )/( ( 3.0*gas_gamma - 1.0 )*( 2.0*gas_gamma + 1.0 ) );
};

double Nu2() {
	return 5.0*( gas_gamma - 1.0 )/( 2.0*gas_gamma + 1.0 );
};

double Nu3() {
	return 3.0/( 2.0*gas_gamma + 1.0 );
};

double Nu4() {
	return -Nu1( gas_gamma )/( 2.0 - gas_gamma );
};

double Nu5() {
	return -2.0/( 2.0 - gas_gamma );
};

// Implicit solution for the velocity function, V(xi)
// W = (gamma + 1/gamma - 1) * (gamma*V - 1)
double W(  double V ) {
	return ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * V - 1.0 );
}

double Vmin() {
	return 1.0/gas_gamma;
}

double Xi( double V ) {
	return Xi_W( W( V ) );
}

double Xi_W( double W ) {
	double V = Vmin() + W;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	// double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * V - 1.0 );
	double Y = ( gas_gamma + 1.0 )*V/2.0;

	return pow( Y, -2.0/5.0 ) * pow( X, Nu1( gas_gamma )/5.0 ) * pow( W, Nu2( gas_gamma )/5.0 );
}

/* Needed to transform xi integrals into V integrals and to do newton solves */

double dXidV( double V ) {
	return dXidV_W( W( V ) );
}

double dXidV_W( double W ) {
	double V = Vmin() + W;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	// double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * V - 1.0 );
	double Y = ( gas_gamma + 1.0 )*V/2.0;

	double XiValue = Xi_W( V, W );

	double dYdV = ( gas_gamma + 1.0 )/2.0;
	double dWdV = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*gas_gamma;
	double dXdV = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 1.0 - 3.0 * gas_gamma  );

	return ( dYdV*( -2.0/5.0 )/Y + dXdV*( Nu1( gas_gamma )/5.0 )/X + dWdV*( Nu2( gas_gamma )/5.0 )/W )*XiValue;
}


double VofXi( double XiValue )
{
	double VMin = 1.0/gas_gamma;
	double eps = 1e-10;
	if ( Xi( VMin + eps ) < XiValue ) {
		// Just newton solve with Xi
		return NewtonRootSolve( 
	}
}

// Given V, solution for G & Z functions

double G( double V ) {
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * V - 1.0 );
	double Y = ( gas_gamma + 1.0 )*V/2.0;
	return ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*pow( W, Nu3( gas_gamma ) ) * pow( X, Nu4( gas_gamma ) ) * pow( ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*( 1.0 - V ), Nu5( gas_gamma ) );
}

double Z( double V, double gas_gamma ) {
	return gas_gamma*( gas_gamma - 1.0 )*( 1.0 - V )*V*V/( 2.0 * ( gas_gamma * V - 1.0 ) );
}

int main( int, char* )
{
	// Plot solution from xi = 0 to 1
	gas_gamma = 5.0/3.0;
	return 0;
}



