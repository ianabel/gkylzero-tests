#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

#define _GNU_SOURCE
#include <fenv.h>

#include "TaylorSedov.h"

#ifndef M_PI
#define M_PI 3.14159265358979323844
#endif


typedef struct TaylorSedovTestCTX_t {
	TaylorSedovProblem *p;
	double r_min;
	double dr;
	double tZero;
} TaylorSedovTestCTX;

#ifdef USE_BC_FUNC
void evalTaylorSedovBC( double t_sim, int nc, const double * skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  TaylorSedovTestCTX *tsCTX = ( TaylorSedovTestCTX* )ctx;
  TaylorSedovProblem *p = ( ( TaylorSedovTestCTX* )ctx )->p;
  double r = tsCTX->r_min - tsCTX->dr/2.0;
  double t = t_sim + tsCTX->tZero;
  double u = TaylorSedovU( t, r, p );
  double rho = TaylorSedovRho( t, r, p );
  double P = TaylorSedovP( t, r, p );

  RHO( ghost )        = rho;
  RHO_UR( ghost )     = rho*u;
  RHO_UTHETA( ghost ) = 0.0;
  RHO_UPHI( ghost )   = 0.0;
  ENERGY( ghost )     = P/( p->gas_gamma - 1.0 ) + ( 1.0/2.0 )*rho*u*u;
}
#endif

void evalTaylorSedovInit(double t_sim, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  TaylorSedovTestCTX *tsCTX = ( TaylorSedovTestCTX* )ctx;
  TaylorSedovProblem *p = ( ( TaylorSedovTestCTX* )ctx )->p;

  double t = t_sim + tsCTX->tZero;
  double r = rCoordinate( xn );
  double u = TaylorSedovU( t, r, p );
  double rho = TaylorSedovRho( t, r, p );
  double P = TaylorSedovP( t, r, p );

  RHO( fout )        = rho;
  RHO_UR( fout )     = rho*u;
  RHO_UTHETA( fout ) = 0.0;
  RHO_UPHI( fout )   = 0.0;
  ENERGY( fout )     = P/( p->gas_gamma - 1.0 ) + ( 1.0/2.0 )*rho*u*u;
}

// map (r,theta) -> (x,y)
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = rCoordinate( xc ), theta = thetaCoordinate( xc );

  // Physical space is always (x,y,z)
  //
  xp[ 0 ] = r * cos( theta ); 

  xp[ 1 ] = r * sin( theta );

}

int main(int argc, char **argv)
{
  // feenableexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 64);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 2);

  double theta = 0.01; // wedge angle

  double r_min = 0.10;
  double r_max = 1.38;

  double dr = ( r_max - r_min )/NX;

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  TaylorSedovProblem testProblem = {
	  .dimension = 2,
	  .rhoZero = 1.0,
	  .InjectedEnergy = 1.0,
	  .gas_gamma = 5.0/3.0
  };

  TaylorSedovTestCTX context = {
	  .r_min = 0.10,
	  .dr = dr,
	  .p = &testProblem,
  };

  // Precompute the proportionality constant
  SetAlpha( &testProblem );

  // Work out InjectedEnergy by picking where the shockwave should be at t = tEnd
  double tEnd  = 0.1; // In seconds since the explosion
  double RZero = 0.4; 
  double REnd  = 1.2;

  // Work out what the initial time-since-explosion is
  // from the requirement that the shock end up at a given position at a given
  // physical time
  double tZero = pow( RZero/REnd, 5.0/2.0 ) * tEnd;
  context.tZero = tZero;

  // Now work out energy
  testProblem.InjectedEnergy = pow( RZero / testProblem.alpha, 5.0 ) * testProblem.rhoZero / ( tZero*tZero );


  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new( testProblem.gas_gamma );

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .evolve = 1,
    .ctx = &context,
    .init = evalTaylorSedovInit,

#ifdef USE_BC_FUNC
    .bc_lower_func = evalTaylorSedovBC,
    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_COPY },
#else
    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
#endif
    .bcy = { GKYL_SPECIES_WEDGE, GKYL_SPECIES_WEDGE },
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "euler_taylorsedov_test2d",

    .ndim = 2,
    // grid in computational space
    .lower = { r_min, -theta/2 },
    .upper = { r_max,  theta/2 },
    .cells = { NX, NY },

    .mapc2p = mapc2p, // mapping of computational to physical space

    .cfl_frac = 0.1,

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.1;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_wv_eqn_release(euler);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
