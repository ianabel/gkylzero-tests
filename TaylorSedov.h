#ifndef TAYLORSEDOV_H
#define TAYLORSEDOV_H

#define R_IDX 0
#define THETA_IDX 1

#define PHI_IDX 2

#define rCoordinate( xn ) xn[ R_IDX ]
#define thetaCoordinate( xn ) xn[ THETA_IDX ]
#define phiCoordinate( xn ) xn[ PHI_IDX ]


#define RHO_IDX 0
#define RHO_UR_IDX 1 + R_IDX
#define RHO_UTHETA_IDX 1 + THETA_IDX
#define RHO_UPHI_IDX 1 + PHI_IDX
#define ENERGY_IDX 4

#define RHO( out ) out[ RHO_IDX ]
#define RHO_UR( out ) out[ RHO_UR_IDX ]
#define RHO_UTHETA( out ) out[ RHO_UTHETA_IDX ]
#define RHO_UPHI( out ) out[ RHO_UPHI_IDX ]
#define ENERGY( out ) out[ ENERGY_IDX ]

/*
 * This test sets up the spherically-symmettric
 * Taylor-Sedov blast wave into a constant-density medium
 */

typedef struct TaylorSedovProblem_t {
	double rhoZero; // External density
	double InjectedEnergy;
	double gas_gamma;
	double alpha;
} TaylorSedovProblem;

void SetAlpha( TaylorSedovProblem * );

double TaylorSedovR( double , TaylorSedovProblem * );
double TaylorSedovRho( double , double , TaylorSedovProblem * );
double TaylorSedovU( double, double, TaylorSedovProblem * );
double TaylorSedovP( double, double, TaylorSedovProblem * );


#endif // TAYLORSEDOV_H
