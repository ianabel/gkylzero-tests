
#include <stdio.h>
#include "TaylorSedov.h"

double Nu1( TaylorSedovProblem * );
double Nu2( TaylorSedovProblem * );
double Nu3( TaylorSedovProblem * );
double Nu4( TaylorSedovProblem * );
double Nu5( TaylorSedovProblem * );

int main( int, char ** )
{
	TaylorSedovProblem testProblem = {
		.gas_gamma = 5.0/3.0,
		.rhoZero = 1.0,
		.InjectedEnergy = 1.0
	};

	SetAlpha( &testProblem );

	printf( "Nu3 = %.8f\n", Nu3( &testProblem ) );
	printf( "Nu4 = %.8f\n", Nu4( &testProblem ) );
	printf( "Nu5 = %.8f\n", Nu5( &testProblem ) );
	return 0;
}

