
#define R_IDX 0
#define THETA_IDX 1
#define PHI_IDX 2

#define rCoordinate( xn ) xn[ R_IDX ]
#define thetaCoordinate( xn ) xn[ THETA_IDX ]
#define phiCoordinate( xn ) xn[ PHI_IDX ]


#define RHO_IDX 0
#define RHO_UX_IDX 1
#define RHO_UY_IDX 2
#define RHO_UZ_IDX 3
#define INT_ENG_IDX 4

#define RHO( out ) out[ RHO_IDX ]
#define RHO_UX( out ) out[ RHO_UX_IDX ]
#define RHO_UY( out ) out[ RHO_UY_IDX ]
#define RHO_UZ( out ) out[ RHO_UZ_IDX ]
#define INT_ENG( out ) out[ INT_ENG_IDX ]

/*
 * This test sets up the spherically-symmettric
 * Taylor-Sedov blast wave into a constant-density medium
 */

struct TaylorSedovProblem_t {
	double rhoZero; // External density
	double InjectedEnergy;
} TaylorSedovProblem;


