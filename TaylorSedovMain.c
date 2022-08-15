
#include <stdio.h>
#include "TaylorSedov.h"

int main( int, char ** )
{
	TaylorSedovProblem testProblem = {
		.gas_gamma = 5.0/3.0,
		.rhoZero = 1.0,
		.InjectedEnergy = 1.0
	};

	SetAlpha( &testProblem );
	printf( "Alpha = %.8f\n", testProblem.alpha );
	printf( "Shock location at t = 0.1 is R = %.8f\n", TaylorSedovR( 0.1, &testProblem ) );
	return 0;
}

