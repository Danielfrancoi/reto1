/*
 * MPI-parallel Jacobi solver for a 1D Poisson-like problem
 * Based on the sequential jacobiSecuencial.c, extended with MPI communication.
 *
 * Compile with: mpicc -O2 -o jacobi_mpi jacobi_mpi.c
 * Run with: mpirun -np <num_procs> ./jacobi_mpi <n> <nsteps> [output_filename]
 *
 *   - n:       Number of intervals (grid points = n+1, boundary at 0 and n)
 *   - nsteps:  Number of Jacobi iterations (each iteration does two half-sweeps)
 *   - output_filename (optional): if provided, rank 0 will write the final solution
 *                                 to this file as two columns: x u(x)
 *
 * Note: Dirichlet boundary conditions u(0)=u(n)=0 are assumed.
 *       The right-hand side f[i] = (i * h), for i=0..n, where h=1.0/n.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    int rank, size;
    int n, nsteps;
    char *fname = NULL;
    double h, h2;
    int *counts = NULL, *displs = NULL;
    int total_interior, base, rem;
    int local_n, g_start, g_end;
    double *u_local = NULL, *utmp_local = NULL, *f_local = NULL;
    double tstart, tend, exec_time;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Process command-line arguments (all ranks see the same argv) */
    n      = (argc > 1) ? atoi(argv[1]) : 100;
    nsteps = (argc > 2) ? atoi(argv[2]) : 100;
    if (argc > 3) {
        fname = argv[3];
    }

    /* 
     * We require at least one interior point per rank:
     * total_interior = n - 1. If size > total_interior, abort.
     */
    total_interior = n - 1;
    if (rank == 0) {
        if (size > total_interior) {
            fprintf(stderr, "Error: number of processes (%d) must be <= n-1 (%d)\n",
                    size, total_interior);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (size > total_interior) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /* 
     * Compute the decomposition of the interior indices [1..n-1].
     * counts[r] = number of interior points on rank r
     * displs[r] = offset (in interior index space, zero-based) where rank r begins
     */
    counts = (int *)malloc(size * sizeof(int));
    displs = (int *)malloc(size * sizeof(int));
    base = total_interior / size;
    rem  = total_interior % size;
    if (rank == 0) {
        int offset = 0;
        for (int r = 0; r < size; ++r) {
            counts[r] = base + (r < rem ? 1 : 0);
            displs[r] = offset;
            offset += counts[r];
        }
    }
    /* Broadcast counts[] and displs[] to all ranks */
    MPI_Bcast(counts, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);

    /* Each rank determines its local number of interior points and global start index */
    local_n = counts[rank];
    g_start = 1 + displs[rank];          /* global index of u_local[1] */
    g_end   = g_start + local_n - 1;      /* global index of u_local[local_n] */

    /* Grid spacing */
    h  = 1.0 / n;
    h2 = h * h;

    /* Allocate local arrays of length local_n + 2 (including two ghost cells) */
    u_local    = (double *)malloc((local_n + 2) * sizeof(double));
    utmp_local = (double *)malloc((local_n + 2) * sizeof(double));
    f_local    = (double *)malloc((local_n + 2) * sizeof(double));

    /* Initialize interior u-values to zero, set ghost boundaries to zero if at domain edge */
    for (int i = 1; i <= local_n; ++i) {
        u_local[i]    = 0.0;
        utmp_local[i] = 0.0;
    }
    if (rank == 0) {
        u_local[0]    = 0.0;   /* left Dirichlet boundary */
        utmp_local[0] = 0.0;
    }
    if (rank == size - 1) {
        u_local[local_n + 1]    = 0.0;   /* right Dirichlet boundary */
        utmp_local[local_n + 1] = 0.0;
    }

    /* Initialize right-hand side f_local for interior points */
    /* f[i] = (global_index) * h */
    for (int i = 1; i <= local_n; ++i) {
        int gi = g_start + (i - 1);
        f_local[i] = gi * h;
    }

    /* Synchronize before starting timing */
    MPI_Barrier(MPI_COMM_WORLD);
    tstart = MPI_Wtime();

    /* ----- Jacobi iterations (each iteration does two half-sweeps) ----- */
    for (int sweep = 0; sweep < nsteps; sweep += 2) {
        /* --- First half-sweep: compute utmp from u_local --- */

        /* Exchange u_local boundary values with neighbors to fill ghost cells */
        /* Send u_local[1] to rank-1, receive into u_local[0] from rank-1 */
        if (rank > 0) {
            MPI_Sendrecv(&u_local[1], 1, MPI_DOUBLE, rank - 1, 0,
                         &u_local[0], 1, MPI_DOUBLE, rank - 1, 1,
                         MPI_COMM_WORLD, &status);
        } else {
            u_local[0] = 0.0;  /* left boundary */
        }
        /* Send u_local[local_n] to rank+1, receive into u_local[local_n+1] from rank+1 */
        if (rank < size - 1) {
            MPI_Sendrecv(&u_local[local_n], 1, MPI_DOUBLE, rank + 1, 1,
                         &u_local[local_n + 1], 1, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, &status);
        } else {
            u_local[local_n + 1] = 0.0;  /* right boundary */
        }

        /* Compute utmp_local for interior points */
        for (int i = 1; i <= local_n; ++i) {
            utmp_local[i] = 0.5 * (u_local[i - 1] + u_local[i + 1] + h2 * f_local[i]);
        }
        /* Set boundary ghost in utmp_local if on domain edge */
        if (rank == 0) {
            utmp_local[0] = 0.0;
        }
        if (rank == size - 1) {
            utmp_local[local_n + 1] = 0.0;
        }

        /* --- Second half-sweep: compute u_local from utmp_local --- */

        /* Exchange utmp_local boundary values with neighbors */
        if (rank > 0) {
            MPI_Sendrecv(&utmp_local[1], 1, MPI_DOUBLE, rank - 1, 0,
                         &utmp_local[0], 1, MPI_DOUBLE, rank - 1, 1,
                         MPI_COMM_WORLD, &status);
        } else {
            utmp_local[0] = 0.0;
        }
        if (rank < size - 1) {
            MPI_Sendrecv(&utmp_local[local_n], 1, MPI_DOUBLE, rank + 1, 1,
                         &utmp_local[local_n + 1], 1, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, &status);
        } else {
            utmp_local[local_n + 1] = 0.0;
        }

        /* Compute updated u_local for interior points */
        for (int i = 1; i <= local_n; ++i) {
            u_local[i] = 0.5 * (utmp_local[i - 1] + utmp_local[i + 1] + h2 * f_local[i]);
        }
        /* Boundary ghosts in u_local remain set for next iteration */
    }
    /* ------------------------------------------------------------------- */

    /* Synchronize and end timing */
    MPI_Barrier(MPI_COMM_WORLD);
    tend = MPI_Wtime();
    exec_time = tend - tstart;

    /* Rank 0 prints execution time (matching sequential: "n,exec_time") */
    if (rank == 0) {
        printf("%d,%f\n", n, exec_time);
    }

    /* If an output filename is provided, gather solution and write it */
    if (fname != NULL) {
        /* Gather only the interior u-values from all ranks to rank 0 */
        double *u_interior = NULL;
        int *recvcounts = NULL;
        int *recvdispls = NULL;

        if (rank == 0) {
            u_interior  = (double *)malloc((n - 1) * sizeof(double));
            recvcounts  = (int *)malloc(size * sizeof(int));
            recvdispls  = (int *)malloc(size * sizeof(int));
            for (int r = 0; r < size; ++r) {
                recvcounts[r] = counts[r];
                recvdispls[r] = displs[r];
            }
        }

        /* Each rank sends its u_local[1..local_n] */
        MPI_Gatherv(
            &u_local[1],           /* send buffer */
            local_n,               /* send count */
            MPI_DOUBLE,
            u_interior,            /* recv buffer on rank 0 */
            recvcounts,            /* recv counts */
            recvdispls,            /* recv displacements */
            MPI_DOUBLE,
            0,                     /* root */
            MPI_COMM_WORLD
        );

        if (rank == 0) {
            /* Reconstruct full solution including boundaries */
            double *u_global = (double *)malloc((n + 1) * sizeof(double));
            u_global[0] = 0.0;
            for (int i = 1; i < n; ++i) {
                u_global[i] = u_interior[i - 1];
            }
            u_global[n] = 0.0;

            /* Write (x, u) to file */
            FILE *fp = fopen(fname, "w");
            if (fp == NULL) {
                fprintf(stderr, "Error: could not open file %s for writing\n", fname);
            } else {
                for (int i = 0; i <= n; ++i) {
                    double x = i * h;
                    fprintf(fp, "%g %g\n", x, u_global[i]);
                }
                fclose(fp);
            }

            free(u_global);
            free(u_interior);
            free(recvcounts);
            free(recvdispls);
        }
    }

    /* Clean up */
    free(u_local);
    free(utmp_local);
    free(f_local);
    free(counts);
    free(displs);

    MPI_Finalize();
    return 0;
}
