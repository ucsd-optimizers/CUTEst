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

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cutest.h"
#include "asa_user.h"

  /* prototypes */
  double asa_value(asa_objective *asa);
  void asa_grad(asa_objective *asa);
  double asa_valgrad(asa_objective *asa);

  /* global variables */
  integer CUTEst_nvar;        /* number of variables */
  integer CUTEst_ncon;        /* number of constraints */

  /* main program */
  int MAINENTRY(void) {
    char *fname = "OUTSDIF.d"; /* CUTEst data file */
    integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
    integer io_buffer = 11;    /* FORTRAN unit for internal i/o */
    integer iout = 6;          /* FORTRAN unit number for error output */
    integer ierr;              /* Exit flag from OPEN and CLOSE */
    integer status;            /* Exit flag from CUTEst tools */
    double  grad_tol = 1.e-6; /* required gradient tolerance */

    doublereal *x, *bl, *bu;
    char       *pname;
    logical     efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;
    logical     constrained = FALSE_;

    doublereal  calls[7], cpu[2];
    int         i, status_asa;

    FILE *spec;
    asa_stat Stats;
    asacg_parm cgParm;
    asa_parm asaParm;

    /* Open problem description file OUTSDIF.d */
    ierr = 0;

    FORTRAN_open(&funit, fname, &ierr);
    if (ierr) {
      printf("Error opening file OUTSDIF.d.\nAborting.\n");
      exit(1);
    }

    /* Determine problem size */
    CUTEST_cdimen(&status, &funit, &CUTEst_nvar, &CUTEst_ncon);
    if (status) {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }
    /* Determine whether to call constrained or unconstrained tools */
    if (CUTEst_ncon) constrained = TRUE_;

    /* stop if the problem has constraints */
    if (constrained) {
      printf(" ** the problem %s has %i constraints\n",
	     pname,  CUTEst_ncon);
      printf("    cg_descent is for unconstrained optimization\n");
      exit(-1);
    }

    /* Reserve memory for variables, bounds, and multipliers */
    /* and call appropriate initialization routine for CUTEst */
    MALLOC(x,  CUTEst_nvar, doublereal);
    MALLOC(bl, CUTEst_nvar, doublereal);
    MALLOC(bu, CUTEst_nvar, doublereal);
    CUTEST_usetup(&status, &funit, &iout, &io_buffer, &CUTEst_nvar,
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

    /* Set any parameter values here
       First read in the default parameter values */
    asa_cg_default (&cgParm);
    asa_default (&asaParm);

    /* if you want to change parameters, change them here: */
    cgParm.PrintParms = TRUE;
    cgParm.PrintLevel = 0;
    asaParm.PrintParms = TRUE;
    asaParm.PrintLevel = 0;

    /* Call ASA_CG */
    status_asa = asa_cg(x, bl, bu, CUTEst_nvar, &Stats, &cgParm, &asaParm,
			grad_tol, asa_value, asa_grad, asa_valgrad, NULL, NULL);

    /* Get CUTEst statistics */
    CUTEST_creport(&status, calls, cpu);
    if (status) {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }

    printf(" *********************** CUTEst statistics ************************\n");
    printf(" Code used               : asa_cg\n");
    printf(" Problem                 : %-s\n", pname);
    printf(" # variables             = %-10d\n", CUTEst_nvar);
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

  CUTEST_ufn(&status, &CUTEst_nvar, asa->x, &f);
  if ((status == 1) || (status == 2)) {
    printf("** CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }

  return (f);
}

void asa_grad(asa_objective *asa) {
  integer status;
  CUTEST_ugr(&status, &CUTEst_nvar, asa->x, asa->g);
  if ((status == 1) || (status == 2)) {
    printf("** CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }
}

double asa_valgrad(asa_objective *asa) {
  logical grad;
  double f;
  integer status;
  grad = 1;
  CUTEST_uofg(&status, &CUTEst_nvar, asa->x, &f, asa->g, &grad);
  if ((status == 1) || (status == 2)) {
    printf("** CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }
  return (f);
}
