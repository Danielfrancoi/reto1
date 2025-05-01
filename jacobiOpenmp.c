#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h> // Inclusión de la biblioteca OpenMP

void jacobi(int nsweeps, int n, double *u, double *f) {
  int i, sweep;
  double h = 1.0 / n;
  double h2 = h * h;
  double *utmp = (double *)malloc((n + 1) * sizeof(double));

  /* Fill boundary conditions into utmp */
  utmp[0] = u[0];
  utmp[n] = u[n];

  for (sweep = 0; sweep < nsweeps; sweep += 2) {
    /* Old data in u; new data in utmp */
    #pragma omp parallel for shared(u, utmp, f, h2) private(i)
    for (i = 1; i < n; ++i)
      utmp[i] = (u[i - 1] + u[i + 1] + h2 * f[i]) / 2;

    /* Old data in utmp; new data in u */
    #pragma omp parallel for shared(u, utmp, f, h2) private(i)
    for (i = 1; i < n; ++i)
      u[i] = (utmp[i - 1] + utmp[i + 1] + h2 * f[i]) / 2;
  }

  free(utmp);
}

void write_solution(int n, double *u, const char *fname) {
  int i;
  double h = 1.0 / n;
  FILE *fp = fopen(fname, "w+");
  for (i = 0; i <= n; ++i)
    fprintf(fp, "%g %g\n", i * h, u[i]);
  fclose(fp);
}

int main(int argc, char **argv) {
  int i, n, nsteps, num_threads;
  double *u, *f, h;
  char *fname;
  struct timeval tstart, tend;

  // Process input arguments
  n = (argc > 1) ? atoi(argv[1]) : 100;
  nsteps = (argc > 2) ? atoi(argv[2]) : 100;
  fname = (argc > 3) ? argv[3] : NULL;
  num_threads = (argc > 4) ? atoi(argv[4]) : 4; //omp_get_max_threads();
  h = 1.0 / n;

  // Set the number of threads for OpenMP
  omp_set_num_threads(num_threads);

  // Allocate and initialize arrays
  u = (double *)malloc((n + 1) * sizeof(double));
  f = (double *)malloc((n + 1) * sizeof(double));

  // Initialize the solution array with zeros
  memset(u, 0, (n + 1) * sizeof(double));

  // Set the right-hand side of the equation
  #pragma omp parallel for shared(f, h) private(i)
  for (i = 0; i <= n; ++i)
    f[i] = i * h;

  // Run the solver
  gettimeofday(&tstart, NULL);
  jacobi(nsteps, n, u, f);
  gettimeofday(&tend, NULL);

  // Calculate the execution time
  double exec_time;
  exec_time =
      (tend.tv_sec - tstart.tv_sec) + 1e-6 * (tend.tv_usec - tstart.tv_usec);

  // Output the execution time and number of threads used
  printf("%d,%f,%d\n", n, exec_time, num_threads);

  if (fname)
    // Write the solution to the specified file
    write_solution(n, u, fname);

  // Free the arrays
  free(f);
  free(u);
  return 0;
}