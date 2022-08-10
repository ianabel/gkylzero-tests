
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "TaylorSedov.h"

double GSLNewtonSolve( double x0, gsl_function_fdf *func, double absTol, double relTol, unsigned int nMaxIters )
{
	const gsl_root_fdfsolver_type * solverType = gsl_root_fdfsolver_newton;
	gsl_root_fdfsolver * solver = gsl_root_fdfsolver_alloc ( solverType );
	
	int err = gsl_root_fdfsolver_set( solver, func, x0 );

	unsigned int i;
	double x = x0, xLast = x0;
	for ( i=0; i < nMaxIters; ++i ) {
		err = gsl_root_fdfsolver_iterate( solver );

		xLast = x;
		x = gsl_root_fdfsolver_root( solver );
		int status = gsl_root_test_delta( x, xLast, absTol, relTol );
		if ( status == GSL_SUCCESS )
			break;
	}

	if ( i == nMaxIters ) {
		printf( "Maximum iterations exceeded in Newton Solve\n" );
	}

	gsl_root_fdfsolver_free( solver );

	return x;
}

// Exponents arising in the solution

double Nu1( TaylorSedovProblem* p ) {
	double gas_gamma = p->gas_gamma;
	return -( 13.0*gas_gamma*gas_gamma - 7.0*gas_gamma + 12.0 )/( ( 3.0*gas_gamma - 1.0 )*( 2.0*gas_gamma + 1.0 ) );
};

double Nu2( TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	return 5.0*( gas_gamma - 1.0 )/( 2.0*gas_gamma + 1.0 );
};

double Nu3( TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	return 3.0/( 2.0*gas_gamma + 1.0 );
};

double Nu4( TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	return -Nu1( p )/( 2.0 - gas_gamma );
};

double Nu5( TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	return -2.0/( 2.0 - gas_gamma );
};

// Implicit solution for the velocity function, V(xi)
// W = (gamma + 1/gamma - 1) * (gamma*V - 1)
double W(  double V, TaylorSedovProblem* p ) {
	double gas_gamma = p->gas_gamma;
	return ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * V - 1.0 );
}

double Vmin( TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	return 1.0/gas_gamma;
}

double Xi_offset( double dV, TaylorSedovProblem* p ) {
	double gas_gamma = p->gas_gamma;
	double V = Vmin( p ) + dV;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * dV );
	double Y = ( gas_gamma + 1.0 )*V/2.0;

	return pow( Y, -2.0/5.0 ) * pow( X, Nu1( p )/5.0 ) * pow( W, Nu2( p )/5.0 );
}

double Xi( double V, TaylorSedovProblem* p ) {
	double gas_gamma = p->gas_gamma;
	return Xi_offset( V - Vmin( p ), p );
}

/* Needed to transform xi integrals into V integrals and to do newton solves */

double dXidV_offset( double dV, TaylorSedovProblem* p ) {
	double gas_gamma = p->gas_gamma;
	double V = Vmin( p ) + dV;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * dV );
	double Y = ( gas_gamma + 1.0 )*V/2.0;

	double XiValue = Xi_offset( dV, p );

	double dYdV = ( gas_gamma + 1.0 )/2.0;
	double dWdV = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*gas_gamma;
	double dXdV = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 1.0 - 3.0 * gas_gamma  );

	return ( dYdV*( -2.0/5.0 )/Y + dXdV*( Nu1( p )/5.0 )/X + dWdV*( Nu2( p )/5.0 )/W )*XiValue;
}

double dXidV( double V, TaylorSedovProblem* p ) {
	return dXidV_offset( V - Vmin( p ), p );
}

void XdX( double V, TaylorSedovProblem *p, double *xi, double *dxi ) 
{
	*xi = Xi( V, p );
	*dxi = dXidV( V, p );
}

double VofXi( double XiValue, TaylorSedovProblem* p )
{
	double gas_gamma = p->gas_gamma;
	double VMin = 1.0/gas_gamma;
	double VMax = 2.0/( gas_gamma + 1.0 );
	double eps = 1e-10;
	if ( Xi( VMin + eps, p ) < XiValue ) {
		gsl_function_fdf FDF = {
			.f = Xi,
			.df = dXidV,
			.fdf = XdX,
			.params = p
		};
		return GSLNewtonSolve( ( VMin + VMax )/2.0, &FDF, 1e-14, 1e-10, 50 );
	} else {
		printf( "Wargl\n" );
		return -1;
	}
}

// Given V, solution for G & Z functions

double G( double V, TaylorSedovProblem* p ) {
	double gas_gamma = p->gas_gamma;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * V - 1.0 );
	double Y = ( gas_gamma + 1.0 )*V/2.0;
	return ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*pow( W, Nu3( p ) ) * pow( X, Nu4( p ) ) * pow( ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*( 1.0 - V ), Nu5( p ) );
}

double Z( double V, TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	return gas_gamma*( gas_gamma - 1.0 )*( 1.0 - V )*V*V/( 2.0 * ( gas_gamma * V - 1.0 ) );
}

int main( int, char** )
{
	TaylorSedovProblem problem = {
		.rhoZero = 1.0,
		.InjectedEnergy= 1.0,
		.gas_gamma = 5.0/3.0
	};

	printf( "Xi(0.675) = %.8f\ndXidV(0.675) = %.8f", Xi( .675, &problem ), dXidV( .675, &problem ) );
	printf( "V(xi = 0.5) = %.8f\n", VofXi( 0.5, &problem ) );
	return 0;
}



