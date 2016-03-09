#include <math.h>
#include <mpi.h>     /* for setting up an MPI_Datatype for PRNG state */
#include <stdbool.h>
#include <stdint.h>  /* for `uint32_t` */
#include <string.h>  /* for `memcpy` */

/* Pi and its square root aren't defined in C99, so if they're not already
   defined, define them now. */
#ifndef PI
#define PI 3.14159265358979
#endif
#ifndef SQRT_PI
#define SQRT_PI 1.77245385090552
#endif

/* Pseudorandom number generator state for PRNG functions. */
typedef struct PRNG {
	uint32_t a;
	uint32_t b;
	uint32_t c;
	uint32_t d;
	uint32_t index;
	uint32_t state[16];
} PRNG;

void init_prng(PRNG* state, uint32_t* seed);
long double rnd(PRNG* state, long double maximum);
void nrv_pair(PRNG* state, double sd, double* nrv1, double* nrv2);
unsigned int prv(PRNG* state, double lambda);
MPI_Datatype prng_mpi_type(void);
void copy_prng(PRNG* dest, const PRNG* src);
bool sanitize_prng(PRNG* state);

double u_no_flow_icdf(double chi, double q);
double u_icdf(double U, double chi, double q);
