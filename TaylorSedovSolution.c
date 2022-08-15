
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

#include "TaylorSedov.h"

double GSLRootSolve( double x_l, double x_u, gsl_function *func, double absTol, double relTol, unsigned int nMaxIters )
{
	const gsl_root_fsolver_type * solverType = gsl_root_fsolver_brent;
	gsl_root_fsolver * solver = gsl_root_fsolver_alloc ( solverType );
	
	int err = gsl_root_fsolver_set( solver, func, x_l, x_u );

	unsigned int i;
	double x_lower = x_l;
	double x_upper = x_u;
	for ( i=0; i < nMaxIters; ++i ) {
		err = gsl_root_fsolver_iterate( solver );
		
		x_lower = gsl_root_fsolver_x_lower ( solver );
		x_upper = gsl_root_fsolver_x_upper ( solver );

		int status = gsl_root_test_interval( x_lower, x_upper, absTol, relTol );

		if ( status == GSL_SUCCESS )
			break;
	}


	double root = gsl_root_fsolver_root( solver );

	if ( i == nMaxIters ) {
		printf( "Maximum iterations exceeded in Newton Solve\n" );
	}

	gsl_root_fsolver_free( solver );

	return root;
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

double Vmax( TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	return 2.0/( gas_gamma + 1.0 );
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

double dXi5dV( double dV, TaylorSedovProblem* p )
{
	double gas_gamma = p->gas_gamma;
	double V = Vmin( p ) + dV;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * dV );
	double Y = ( gas_gamma + 1.0 )*V/2.0;

	double Xi5Value = pow( X, Nu1( p ) )*pow( W, Nu2( p ) )/( Y*Y );

	double dYdV = ( gas_gamma + 1.0 )/2.0;
	double dWdV = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*gas_gamma;
	double dXdV = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 1.0 - 3.0 * gas_gamma  );

	return ( -2.0 * dYdV/Y + Nu1( p )*dXdV/X + Nu2( p )*dWdV/W )*Xi5Value;

}

typedef struct VSolve_t {
	TaylorSedovProblem  *p;
	double XiValue;
} VSolve;

double Xi_for_V( double V_Offset, void *ptr )
{
	VSolve *vPtr = ( VSolve * )ptr;
	return Xi_offset( V_Offset, vPtr->p ) - vPtr->XiValue;
}

double VofXi( double XiValue, TaylorSedovProblem* p )
{
	double VMin = Vmin( p );
	double VMax = Vmax( p );

	if ( XiValue == 0.0 )
		return 0.0;
	if ( XiValue == 1.0 )
		return VMax - VMin;
	if ( XiValue < 0.0 || XiValue > 1.0 )
	{
		printf( "PANIC\n" );
		return -1;
	}
	
	VSolve param_data = {
		.p = p,
		.XiValue = XiValue
	};

	gsl_function xiGSL = {
		.function = Xi_for_V,
		.params = &param_data
	};

	double V_offset = GSLRootSolve( 0.0, VMax-VMin, &xiGSL, 1e-8, 1e-6, 75 );

	return V_offset;
	
}

// Given V, solution for G & Z functions

double G_offset( double V_offset, TaylorSedovProblem* p ) {
	double gas_gamma = p->gas_gamma;
	double V = Vmin( p ) + V_offset;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	double W = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) ) * ( gas_gamma * V_offset );
	double U = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*( 1.0 - V );
	return ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*pow( W, Nu3( p ) ) * pow( X, Nu4( p ) ) * pow( U, Nu5( p ) );
}

double Z_offset( double V_offset, TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	double V = Vmin( p ) + V_offset;
	return gas_gamma*( gas_gamma - 1.0 )*( 1.0 - V )*V*V/( 2.0 * ( gas_gamma * V_offset ) );
}

// Return Z(xi)*G(xi)*xi^2 for accuracy near 0
double ZG2_offset( double V_offset, TaylorSedovProblem *p ) {
	double gas_gamma = p->gas_gamma;
	double V = Vmin( p ) + V_offset;
	double X = ( ( gas_gamma + 1.0 )/( 7.0 - gas_gamma ) )*( 5.0 - ( 3.0*gas_gamma - 1.0 )*V );
	double U = ( ( gas_gamma + 1.0 )/( gas_gamma - 1.0 ) )*( 1.0 - V );
	double Y = ( gas_gamma + 1.0 )*V/2.0;

	// Z = 2*gamma * (gamma - 1)/(gamma+1)^2 * U * Y^2 /W
	// G = (gamma + 1)/(gamma - 1) * W^Nu3 X^Nu4 U^Nu5
	// Xi^2 = Y^(-4/5) X ^ (2Nu1/5) * W^(2Nu2/5)
	//
	//
	return ( 2* gas_gamma / ( gas_gamma + 1.0 ) ) * pow( Y, 2.0 - 4.0/5.0 ) * pow( X, Nu4( p ) + 0.4*Nu1( p ) ) * pow( U, Nu5( p ) + 1.0 );
}

/*
 * Normalised density, radial velocity, and pressure
 */

double RhoTilde( double xi, TaylorSedovProblem *p )
{
	double V_offset = VofXi( xi, p );
	return G_offset( V_offset, p );
}

double pTilde( double xi, TaylorSedovProblem *p )
{
	double V_offset = VofXi( xi, p );

	return ZG2_offset( V_offset, p )/p->gas_gamma;
}

double uTilde( double xi, TaylorSedovProblem *p )
{
	double V_offset = VofXi( xi, p );
	return xi*( V_offset + Vmin( p ) );
}

double AlphaIntegrand( double V_offset, void *ptr )
{
	TaylorSedovProblem *p = ( TaylorSedovProblem* )ptr;
	double gas_gamma = p->gas_gamma;
	double V = V_offset + Vmin( p );
	return G_offset( V_offset, p ) * ( V*V/2.0 + Z_offset( V_offset, p )/( gas_gamma * ( gas_gamma - 1.0 ) ) )
	           * ( dXi5dV( V_offset,p )/5.0 );
}

double Alpha( TaylorSedovProblem *p )
{
	gsl_integration_workspace *work = gsl_integration_workspace_alloc( 2048 );

	gsl_function AlphaGSL = {
		.function = AlphaIntegrand,
		.params = p
	};

	int status;
	double alphaIntegral,errorEstimate;

	status = gsl_integration_qags( &AlphaGSL, 0, Vmax( p ) - Vmin( p ), 1e-8, 1e-6, 50, work, &alphaIntegral, &errorEstimate );

	if ( status != GSL_SUCCESS )
	{
		printf( "Failure to integrate\n" );
		return -1;
	}

	gsl_integration_workspace_free( work );

	double alphaFive = 25.0 / ( 16.0 * M_PI * alphaIntegral );
	return pow( alphaFive, 1.0/5.0 );
}

void SetAlpha( TaylorSedovProblem *p )
{
	p->alpha = Alpha( p );
}

// Return Shock Location
double TaylorSedovR( double t, TaylorSedovProblem *p ) {
	return p->alpha * pow( p->InjectedEnergy * t * t / p->rhoZero, 0.2 );
}

double TaylorSedovRho( double t, double r, TaylorSedovProblem *p ) {
	double XiVal = r / TaylorSedovR( t, p );
	if ( XiVal <= 1.0 )
		return p->rhoZero * RhoTilde( XiVal, p );
	else 
		return p->rhoZero;
}

double TaylorSedovU( double t, double r, TaylorSedovProblem *p ) {
	double XiVal = r/ TaylorSedovR( t, p );
	if ( XiVal <= 1.0 ) {
		double Rdot = 2*TaylorSedovR( t, p )/( 5 * t );
		return Rdot * uTilde( XiVal, p );
	} else {
		return 0;
	}
}

double TaylorSedovP( double t, double r, TaylorSedovProblem *p ) {
	double XiVal = r/ TaylorSedovR( t, p );
	if ( XiVal <= 1.0 ) {
		double Rdot = 2*TaylorSedovR( t, p )/( 5 * t );
		return p->rhoZero * Rdot * Rdot * pTilde( XiVal, p );
	} else {
		return EXTERNAL_PRESSURE;
	}
}

