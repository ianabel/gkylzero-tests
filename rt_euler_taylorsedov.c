#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

#include "TaylorSedov.h"


void eval(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double r = rCoordinate( xn );

}

// map (r,theta) -> (x,y)
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = rCoordinate( xc ), theta = thetaCoordinate( xc ), phi = phiCoordinate( xc );

  // Physical space is always (x,y,z)
  //
  xp[ 0 ] = r * sin( theta ) * cos( phi ); 

  xp[ 1 ] = r * sin( theta ) * sin( phi );

  xp[ 2 ] = r * cos( theta );
}

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 64);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 2);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 2);

  double theta = 0.01; // wedge angle

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  TaylorSedovProblem testProblem = {
	  .rhoZero = 1.0;
	  .InjectedEnergy = 1.0;
	  .gas_gamma = 5.0/3.0;
  };

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new( testProblem.gas_gamma );

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .evolve = 1,
    .ctx = &testProblem,
    .init = evalTaylorSedovInit,
	 .bc_lower_func = evalTaylorSedovBC, // { evalTaylorSedovBC, NULL, NULL },

    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_COPY },
	 .bcy = { GKYL_SPECIES_WEDGE, GKYL_SPECIES_WEDGE },
	 // .bcz ignored because set to be periodic later
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "euler_taylorsedov_test",

    .ndim = 3,
    // grid in computational space
    .lower = { 0.1, -theta/2, 0.0 },
    .upper = { 10.0,  theta/2, 2.0*M_PI },
    .cells = { NX, NY, NZ },

    .mapc2p = mapc2p, // mapping of computational to physical space

    .cfl_frac = 0.9,

    .num_periodic_dirs = 1,
    .periodic_dirs = { 2 },
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
