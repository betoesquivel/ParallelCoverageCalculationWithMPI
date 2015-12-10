//-----------------------------------------------------------------------------
// c c o v e r a g e   --   calculate the coverage of feature detectors
//-----------------------------------------------------------------------------
// COMPILE: mpicc -o coverage ParallelCoverage.c
// RUN: mpiexec -machinefile machinefile -n 90 ./coverage
//
// In one node within an MPI cluster with open_mpi installed where:
// - machinefile is the name of a file with all addresses of the nodes
// in the cluster each in a new line
// - 90 is the number of processes you wish to run in your cluster
//-----------------------------------------------------------------------------
// Base source code by: Dr. Adrian Clark
// Reg no. of student: 1500921
// Submission date: 10 December, 2015

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXFN   200
#define MAXLINE 2000
#define MAXLOCS 200000

int procs;
int rank;

int
load_locations (char *fn, double locs[][2], int first)
{
  char line[MAXLINE], *kwd;
  int nlocs;
  FILE *f;

  if ((f = fopen (fn, "r")) == NULL) {
    fprintf (stderr, "Cannot open %s!\n", fn);
    exit (EXIT_FAILURE);
  }
  fgets (line, MAXLINE, f);
  fgets (line, MAXLINE, f);
  nlocs = first;
  while (fgets (line, MAXLINE, f) != NULL) {
    kwd = strtok (line, " ");
    locs[nlocs][1] = atof (kwd);
    kwd = strtok ((char *) NULL, " ");
    locs[nlocs][0] = atof (kwd);
    nlocs += 1;
    if (nlocs >= MAXLOCS) {
      fprintf (stderr, "Not enough space to hold all locations: %d!\n", nlocs);
      exit (EXIT_FAILURE);
    }
  }
  fclose (f);
  return nlocs;
}

/*
 * Sets the value for a lower and upper bound (lb, ub respectively) based on
 * a number of processes procs in which to divide the data with data_length in
 * an array, given a section_length and the number of processes and rank in the
 * cluster.
 */
void
get_bounds(int *lb, int *ub, int section_length, int data_length)
{
  *lb = rank * section_length;
  *ub = (rank >= (procs-1)) ? data_length : (rank+1) * section_length;
}

/*
 * Distributes the calculation of coverage for the data points received.
 * Each cluster iterates over a section of the data points array.
 * The section is delimitted by lowerBoundary and upperBoundary, which
 * are calculated by the routine get_bounds and passed by reference.
 * After this, variables with the partial accumulated data to calculate
 * coverage (hmsum, np_overall, nc) are "reduced" through by adding
 * them altogether and writing the result in the following variables
 * in the root node (rank=0):
 * reduced_hmsum, reduced_np_overall, reduced_nc
 *
 * Finally these variables are used to calculate the coverage results
 * and store them in the variables received by reference:
 * result, npaths, ncoin.
 *
 */
void
coverage (double locs[][2], int nlocs, double *result, int *npaths, int *ncoin)
{
  int np_overall = 0, nc = 0, np, i, j;
  double hmsum = 0.0, dsum, hm, d, y1, x1, y2, x2;

  int section_length = nlocs/procs;
  int lowerBoundary = 0, upperBoundary = 0;
  get_bounds(&lowerBoundary, &upperBoundary, section_length, nlocs);

  for (i = lowerBoundary; i < upperBoundary; i++) {
    y1 = locs[i][0];
    x1 = locs[i][1];
    np = 0;
    dsum = 0.0;
    for (j = 0; j < nlocs; j++) {

      if (i != j) {

        y2 = locs[j][0] - y1;
        x2 = locs[j][1] - x1;
        d = sqrt (y2 * y2 + x2 * x2);
        if (d <= 0.0) {

          nc += 1;

        } else {

          dsum += 1.0 / d;
          np += 1;
          np_overall += 1;

        }

      }
    }
    hm = np / dsum;
    hmsum += 1.0 / hm;

  }

  int reduced_np_overall = 0, reduced_nc = 0;
  double reduced_hmsum = 0.0;

  /*
   * Make sure all processes are synchronized, because MPI_Reduce does not
   * guarantee synchronization.
   */
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce(&hmsum,&reduced_hmsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  MPI_Reduce(&np_overall,&reduced_np_overall,1,
      MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  MPI_Reduce(&nc,&reduced_nc,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  if (rank == 0) {
    *result = nlocs / reduced_hmsum;
    *npaths = reduced_np_overall;
    *ncoin  = reduced_nc;
  } else {
    *result = 1;
    *npaths = 1;
    *ncoin  = 1;
  }
}

double
mean_coverage1 (char *opname, char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  double locs[MAXLOCS][2];
  double csum = 0.0, cov;
  char fn[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn, "data/%s_%s%d.txt", opname, imsets[i], fno);
      nlocs = load_locations (fn, locs, 0);
      coverage (locs, nlocs, &cov, &np, &nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


double
mean_coverage2 (char *opname1, char *opname2, char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  double locs[MAXLOCS][2];
  double csum = 0.0, cov;
  char fn1[MAXFN], fn2[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn1, "data/%s_%s%d.txt", opname1, imsets[i], fno);
      nlocs = load_locations (fn1, locs, 0);
      sprintf (fn2, "data/%s_%s%d.txt", opname2, imsets[i], fno);
      nlocs = load_locations (fn2, locs, nlocs);
      coverage (locs, nlocs, &cov, &np, &nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


double
mean_coverage3 (char *opname1, char *opname2, char *opname3,
		char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  double locs[MAXLOCS][2];
  double csum = 0.0, cov;
  char fn1[MAXFN], fn2[MAXFN], fn3[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn1, "data/%s_%s%d.txt", opname1, imsets[i], fno);
      nlocs = load_locations (fn1, locs, 0);
      sprintf (fn2, "data/%s_%s%d.txt", opname2, imsets[i], fno);
      nlocs = load_locations (fn2, locs, nlocs);
      sprintf (fn3, "data/%s_%s%d.txt", opname3, imsets[i], fno);
      nlocs = load_locations (fn3, locs, nlocs);
      coverage (locs, nlocs, &cov, &np, &nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


double
mean_coverage4 (char *opname1, char *opname2, char *opname3, char *opname4,
		char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  double locs[MAXLOCS][2];
  double csum = 0.0, cov;
  char fn1[MAXFN], fn2[MAXFN], fn3[MAXFN], fn4[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn1, "data/%s_%s%d.txt", opname1, imsets[i], fno);
      nlocs = load_locations (fn1, locs, 0);
      sprintf (fn2, "data/%s_%s%d.txt", opname2, imsets[i], fno);
      nlocs = load_locations (fn2, locs, nlocs);
      sprintf (fn3, "data/%s_%s%d.txt", opname3, imsets[i], fno);
      nlocs = load_locations (fn3, locs, nlocs);
      sprintf (fn4, "data/%s_%s%d.txt", opname4, imsets[i], fno);
      nlocs = load_locations (fn4, locs, nlocs);
      coverage (locs, nlocs, &cov, &np, &nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


int
main (int argc, char **argv)
{
  // Init MPI and get information about the node in the cluster.
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  char *opnames[] = {"ebr", "ibr", "mser", "sfop"};
  int nopnames = 4;
  char *imsets[] = {"bark", "bikes", "boat", "graf", "leuv", "trees", "ubc",
		    "wall"};
  int nimsets = 8;
  int i1, i2, i3, i4;
  double mc;

  // Calculate the coverage for each individual operator.
  for (i1 = 0; i1 < nopnames; i1++) {
    mc = mean_coverage1 (opnames[i1], imsets, nimsets);
    if (rank == 0) {
      printf ("%s: %f\n", opnames[i1], mc);
    }
  }

  // Calculate the coverage for each combination of pairs of operators.
  for (i1 = 0; i1 < nopnames; i1++) {
    for (i2 = i1+1; i2 < nopnames; i2++) {
      mc = mean_coverage2 (opnames[i1], opnames[i2], imsets, nimsets);
      if (rank == 0) {
        printf ("%s + %s: %f\n", opnames[i1], opnames[i2], mc);
      }
    }
  }

  // Calculate the coverage for each combination of triplets of operators.
  for (i1 = 0; i1 < nopnames; i1++) {
    for (i2 = i1+1; i2 < nopnames; i2++) {
      for (i3 = i2+1; i3 < nopnames; i3++) {
        mc = mean_coverage3 (opnames[i1], opnames[i2], opnames[i3],
			     imsets, nimsets);
        if (rank == 0) {
          printf ("%s + %s + %s: %f\n", opnames[i1], opnames[i2],
              opnames[i3], mc);
        }
      }
    }
  }

  // Calculate the coverage for each combination of quadruplets of operators.
  for (i1 = 0; i1 < nopnames; i1++) {
    for (i2 = i1+1; i2 < nopnames; i2++) {
      for (i3 = i2+1; i3 < nopnames; i3++) {
        for (i4 = i3+1; i4 < nopnames; i4++) {
          mc = mean_coverage4 (opnames[i1], opnames[i2], opnames[i3],
                   opnames[i4],  imsets, nimsets);
          if (rank == 0) {
            printf ("%s + %s + %s + %s: %f\n", opnames[i1], opnames[i2],
            opnames[i3], opnames[i4], mc);
          }
        }
      }
    }
  }

  // End MPI Connection
  MPI_Finalize();

  return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
// End of ccoverage.c
//-----------------------------------------------------------------------------
