/* ====================================================
 * CUTEst interface for ASA-CG             Feb 23, 2020
 *
 * W. Hager
 *
 * (Based on CUTEst cg_descent_main.c
 *  of W. Hager, Apr 5, 2014)
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cutest.h"
#include "asa_user.h"

  /* prototypes */
  double asa_value(asa_objective *asa);
  void asa_grad(asa_objective *asa);
  double asa_valgrad(asa_objective *asa);

  /* main program */
  int MAINENTRY(void) {
    char *fname = "OUTSDIF.d"; /* CUTEst data file */
    integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
    integer io_buffer = 11;    /* FORTRAN unit for internal i/o */
    integer iout = 6;          /* FORTRAN unit number for error output */
    integer ierr;              /* Exit flag from OPEN and CLOSE */
    integer status;            /* Exit flag from CUTEst tools */
    double  grad_tol = 1.e-6; /* required gradient tolerance */

    integer     M, N;
    doublereal *x, *bl, *bu;
    char       *pname;
    logical     efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;

    doublereal  calls[7], cpu[2];
    int         i, status_asa;

    FILE *spec;
    asa_stat Stats;
    asacg_parm cgParm;
    asa_parm asaParm;

    FILE *fp;
    struct timeval start, end;
    doublereal elapsedTime;

    /* Open problem description file OUTSDIF.d */
    ierr = 0;

    FORTRAN_open(&funit, fname, &ierr);
    if (ierr) {
      printf("Error opening file OUTSDIF.d.\nAborting.\n");
      exit(1);
    }

    /* Determine problem size */
    CUTEST_cdimen(&status, &funit, &N, &M);
    if (status) {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }
    if (M) {
      /* stop if the problem has constraints */
      printf(" ** the problem %s has %i constraints\n",
	     pname,  M);
      printf("    cg_descent is for unconstrained optimization\n");
      exit(-1);
    }

    /* Reserve memory for variables, bounds, and multipliers */
    /* and call appropriate initialization routine for CUTEst */
    MALLOC(x,  N, doublereal);
    MALLOC(bl, N, doublereal);
    MALLOC(bu, N, doublereal);
    CUTEST_usetup(&status, &funit, &iout, &io_buffer, &N,
		  x, bl, bu);
    if (status) {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }

    /* Get problem name */
    MALLOC(pname,  FSTRING_LEN+1, char);
    CUTEST_probname(&status, pname);
    if (status) {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }

    /* Make sure to null-terminate problem name */
    pname[FSTRING_LEN] = '\0';
    i = FSTRING_LEN - 1;
    while(i-- > 0 && pname[i] == ' ') {
      pname[i] = '\0';
    }

    /* Start wall timer */
    gettimeofday(&start, NULL);

    /* Set any parameter values here
       First read in the default parameter values */
    asa_cg_default (&cgParm);
    asa_default (&asaParm);

    /* if you want to change parameters, change them here: */
    cgParm.PrintParms = FALSE;
    cgParm.PrintLevel = 0;
    asaParm.PrintParms = FALSE;
    asaParm.PrintLevel = 0;
    asaParm.StopRule = FALSE;
    cgParm.maxit = 1000000;

    /* Call ASA_CG */
    status_asa = asa_cg(x, bl, bu, N, &Stats, &cgParm, &asaParm,
			grad_tol, asa_value, asa_grad, asa_valgrad, NULL, NULL);

    /* Stop timer */
    gettimeofday(&end, NULL);
    elapsedTime = end.tv_sec + end.tv_usec / 1e6 -
      start.tv_sec - start.tv_usec / 1e6; /* in seconds */


    /* Get CUTEst statistics */
    CUTEST_creport(&status, calls, cpu);
    if (status) {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }

    /* Print results to file */
    fp = fopen("asa_cg.all", "a");
    fprintf(fp, "%-10s %-10d  %-10li %-10.2f %-10.2f %-10d %-6d %-15.7e %-15.7e %-12.4f %-12.4f\n",
	    pname, N, Stats.cgiter, calls[0], calls[1],
	    N - Stats.nfree, status_asa, Stats.f, Stats.pgnorm, elapsedTime, cpu[1]);
    fclose(fp);

    fp = fopen("asa_cg.pp", "a");
    fprintf(fp, "%-10s %-10d %-15.7g %-12.4f %-12.4f\n",
	    pname, status_asa, calls[0], elapsedTime, cpu[1]);
    fclose(fp);

    printf(" *********************** CUTEst statistics ************************\n");
    printf(" Code used               : asa_cg\n");
    printf(" Problem                 : %-s\n", pname);
    printf(" # variables             = %-10d\n", N);
    printf(" # cg iterations         = %li\n", Stats.cgiter);
    printf(" # objective functions   = %-15.7g\n", calls[0]);
    printf(" # objective gradients   = %-15.7g\n", calls[1]);
    printf(" # objective Hessians    = %-15.7g\n", calls[2]);
    printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
    printf(" Exit code               = %-10d\n", status_asa);
    printf(" Final f                 = %-15.7g\n",Stats.f);
    printf(" Final ||g||             = %-15.7g\n",Stats.pgnorm);
    printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
    printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
    printf(" ******************************************************************\n");

    ierr = 0;
    FORTRAN_close(&funit, &ierr);
    if (ierr) {
      printf("Error closing %s on unit %d.\n", fname, funit);
      printf("Trying not to abort.\n");
    }

    /* Free workspace */
    FREE(pname);
    FREE(x); FREE(bl); FREE(bu);

    CUTEST_uterminate(&status);

    return 0;
  }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif

double asa_value(asa_objective *asa) {
  double f;
  integer status;
  integer N;
  N = asa->n;
  CUTEST_ufn(&status, &N, asa->x, &f);
  if ((status == 1) || (status == 2)) {
    printf("** CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }

  return (f);
}

void asa_grad(asa_objective *asa) {
  integer status;
  integer N;
  N = asa->n;
  CUTEST_ugr(&status, &N, asa->x, asa->g);
  if ((status == 1) || (status == 2)) {
    printf("** CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }
}

double asa_valgrad(asa_objective *asa) {
  logical grad;
  double f;
  integer N;
  integer status;
  grad = 1;
  N = asa->n;
  CUTEST_uofg(&status, &N, asa->x, &f, asa->g, &grad);
  if ((status == 1) || (status == 2)) {
    printf("** CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }
  return (f);
}
