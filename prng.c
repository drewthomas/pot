#include "prng.h"

/* Initialize a PRNG state with a given 512-bit seed. */
void init_prng(PRNG* state, uint32_t* seed)
{
	memcpy(state->state, seed, 16 * sizeof(uint32_t));
	state->index = 0;
}

/* Generate a uniformly distributed random variate between 0 & `maximum`. */
long double rnd(PRNG* state, long double maximum)
{
	/* The internals: this is a WELL512 PRNG. */
	state->a = state->state[state->index];
	state->c = state->state[(state->index + 13) & 15];
	state->b = state->a ^ state->c ^ (state->a << 16) ^ (state->c << 15);
	state->c = state->state[(state->index + 9) & 15];
	state->c ^= (state->c >> 11);
	state->a = state->state[state->index];
	state->state[state->index] = state->b ^ state->c;
	state->d = state->a ^ ((state->a << 5) & 0xda442d24UL);
	state->index = (state->index + 15) & 15;
	state->a = state->state[state->index];
	state->state[state->index] = state->a ^ state->b ^ state->d
	                             ^ (state->a << 2) ^ (state->b << 18)
	                             ^ (state->c << 28);
	return maximum * (state->state[state->index] / (long double) UINT32_MAX);
}

/* Generate a pair of normal random variates with mean zero and store them
   in `nrv1` and `nrv2`. */
void nrv_pair(PRNG* state, double sd, double* nrv1, double* nrv2)
{
	double r;
	double theta;
	long double urv1;
	long double urv2 = rnd(state, 1.0L);

	/* Calculate the variates with the Box-Muller transform's polar form. */
	do {
		urv1 = rnd(state, 1.0L);
	} while (urv1 == 0.0);  /* `rnd` very occasionally returns zero! */
	r = sqrt(-2.0 * (double) logl(urv1));
	theta = 2.0 * PI * urv2;
	*nrv1 = sd * r * cos(theta);
	*nrv2 = sd * r * sin(theta);
}

/* Generate a Poisson-distributed random variate with mean `lambda`. */
unsigned int prv(PRNG* state, double lambda)
{
	double exp_lambda = exp(-lambda);
	unsigned int k = 0;
	double ratio_term = 1.0;
	double uni = rnd(state, 1.0L);

	while (ratio_term) {
		uni -= ratio_term * exp_lambda;
		if (uni < 0.0) {
			break;
		}
		k++;
		ratio_term *= lambda / k;
	}

	if (ratio_term) {
		return k;
	} else {
		return 0;
	}
}

/* Define an MPI data type for the PRNG data structure.
   (Don't forget to free it later!) */
/* FIXME: at some stage I should probably rewrite this to treat the `struct`
   as a bunch of different fields rather than a simple bunch of bytes. The
   bag-of-bytes approach could go wrong if a state file created on one
   platform is loaded on another. */
MPI_Datatype prng_mpi_type(void)
{
/*
	PRNG example;
	MPI_Aint MPIPRNG_disps[6];
	MPI_Datatype MPIPRNG_subtypes[] = {
		MPI_UINT32_T,
		MPI_UINT32_T,
		MPI_UINT32_T,
		MPI_UINT32_T,
		MPI_UINT32_T,
		MPI_UINT32_T
	};

	MPI_Get_address(&(example.a), &(MPIPRNG_disps[0]));
	MPI_Get_address(&(example.b), &(MPIPRNG_disps[0]));
	MPI_Get_address(&(example.c), &(MPIPRNG_disps[0]));
	MPI_Get_address(&(example.d), &(MPIPRNG_disps[0]));
	MPI_Get_address(&(example.index), &(MPIPRNG_disps[0]));
	MPI_Get_address(&(example.state), &(MPIPRNG_disps[0]));
*/
	MPI_Datatype MPIPRNG;
	
	MPI_Type_contiguous(sizeof(PRNG), MPI_BYTE, &MPIPRNG);
	MPI_Type_commit(&MPIPRNG);

	return MPIPRNG;
}

/* Copy the contents of one PRNG data structure into another. */
void copy_prng(PRNG* dest, const PRNG* src)
{
	*dest = *src;
	memcpy(dest->state, src->state, 16 * sizeof(uint32_t));
}

/* Check whether a PRNG state seems invalid, and if so try
   resetting/sanitizing it to a sane state. */
bool sanitize_prng(PRNG* state)
{
	if (state->index >= 16) {
		state->index = state->index % 16;
		return true;
	}

	return false;
}

/* Compute the `q`th quantile of the inverse cumulative distribution
   function of a particle's speed at infinity, according to Hutchinson
   (2003)'s eqs. 18 & 19 in the special case of no flow.

   19/05/2014: I've compared this function's output against SCEPTIC's
   `finvtfunc` and Mathematica's for `chi` = 0; this function matches
   Mathematica's results except at extreme `q` values (where the inverse
   CDF is nearly flat) and differs in the 2nd/3rd decimal place from
   SCEPTIC's results. Since the Mathematica results are derived semi-
   analytically and this function pretty much matches those, the
   SCEPTIC function appears to be the one that's (slightly) wrong here! */
double u_no_flow_icdf(double chi, double q)
{
	double omc = 1.0 - chi;
	double p;
	double u = 1.29951;  /* use (`chi` = 0) median as starting guess */
	/* FIXME: if/when I start using nonzero `chi` for this function,
	   use an approximate expression for the median as a func. of `chi`
	   to initialize `u`. */
	double usq;

	/* The Newton-Raphson method does the inversion; it needs 4.75
	   iterations on average when `chi` is zero. */
	do {
		usq = u * u;
		p = 1.0 - ((usq + omc) * exp(-usq) / omc) - q;
		u -= p / (2.0 * u * exp(-usq) * (chi - usq) / -omc);
	} while ((p * p) > 1e-16);  /* squaring because p may be negative */

	return u;
}

/* Compute the `q`th quantile of the inverse cumulative distribution
   function of a particle's speed at infinity according to Hutchinson
   (2003)'s eq. 19 given nonzero flow (i.e. the normalized velocity `U` > 0).

   Warning: this function's caller must ensure `U` is greater than zero.

   19/05/2014: I've tested this function too against Mathematica and
   SCEPTIC for `chi` = 0 and it appears to work (so long as one bears in
   mind that SCEPTIC's equivalent of `U` is `Uc`, not `vd`). */
/* FIXME: when `chi` = 0, this function works when `U` is less than 26.48.
   I should either make it more robust for larger `U`, or make it explicitly
   detect when `U` is >= 26.48 and fail with a warning. */
double u_icdf(double U, double chi, double q)
{
	const double Usq = U * U;
	const double psi = SQRT_PI * (1.0 + (2.0 * (Usq - chi))) / 2.0;
	const double p_den = 2.0 * ((U * exp(-Usq)) + (psi * erf(U)));

	double dp_den;
	double dp_num;
	double p;
	double p_num;
	double temp;
	double u;
	double usq;

	/* Try to find a good starting estimate for `u` for the Newton-Raphson
	   iterations to come. */
	if (U > 2.1) {
		/* Use the median in the large `U` limit (when `chi` = 0). */
		u = U + (1.0 / U);
	} else {
		/* Use a quadratic (in `U`) approximation for the median `u`
		   (assuming `chi` = 0). */
		u = 1.2586147 + (0.221 * U) + (0.1903312 * U * U);
	}

	/* Precompute the u-independent denominator of the unwieldy fraction
	   in the definition of the cumulative dist. func.'s derivative. */
	dp_den = 2.0 * SQRT_PI * (U + (psi * exp(Usq) * erf(U)));

	/* Do the inversion with the Newton-Raphson method. */
	do {
		usq = u * u;
		temp = - (U + u) * (U + u);
		p_num = (U - u) * exp(temp);
		p_num += (U + u) * exp((4.0 * U * u) + temp);
		p_num += psi * (erf(U - u) + erf(U + u));
		p = 1.0 - (p_num / p_den) - q;
		dp_num = (SQRT_PI * (1.0 + (2.0 * (Usq - usq)))) - (2.0 * psi);
		dp_num *= exp(-u * ((2.0 * U) + u)) - exp((2.0 * U * u) - usq);
		u -= p / (dp_num / dp_den);
	} while ((p * p) > 1e-16);  /* squaring because p may be negative */

	return u;
}

/* Reference:
     I. H. Hutchinson (2003). Plasma Physics and Controlled Fusion, 45(8),
     pp. 1477-1500. */
