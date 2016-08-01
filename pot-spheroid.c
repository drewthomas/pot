#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw.h>

#include <GL/glu.h>


#ifdef FPE
	/* This enables the use of `feenableexcept` and `fedisableexcept`,
	   which make debugging easier by allowing the program to unmask
	   floating point exceptions. Note that this is not ANSI C. */
	#define _GNU_SOURCE
	#include <fenv.h>
#endif

#include "prng.h"

/* If `DISABLE_GRAIN` is defined, the grain doesn't collect particles;
   instead particles pass through it without being absorbed.
   However, note that if the grain's given an initial charge, it still
   interacts with the other particles electromagnetically. */
/*
#define DISABLE_GRAIN
*/

/* The number of radial, azimuthal, and zenithal bins to use if computing
   spatial profiles of particle number and/or electric potential. */
#define CELLS_R 50
#define CELLS_AZI 4
#define CELLS_ZEN 10

/* If `GL` is defined, show a graphical display of the simulation. */
#define GL
/*
*/

/* If `GL` is defined, and the user requests tracking of particle
   trajectories, `TRAJ_LEN` is the "length" or number of data points to
   store from each trajectory, and `TRAJ_ALPHA` is the alpha value used when
   drawing the trajectories. `TRAJ_ALPHA` may be undefined, in which case
   the trails are simply rendered in opaque blue without any alpha value.
   (If `GL` isn't defined `TRAJ_LEN` & `TRAJ_ALPHA` have no effect.) */
#define TRAJ_LEN 5000
#define TRAJ_ALPHA 0.4

/* The x, y & z components respectively of the externally-applied,
   time-constant magnetic field (in tesla). */
#define EXTERNAL_B_X 0.0
#define EXTERNAL_B_Y 0.0
#define EXTERNAL_B_Z 0.0

/* The ratio of the ion mass to the electron mass. */
#define MIME 1836.153

/* Default output directory to use for output files. */
#define OUTPUT_DIR "./output/"

/* If `RECORD_MACRO_DATA` is defined, record macroscopic (i.e. whole-system)
   data in the file at the path `RECORD_MACRO_DATA`. */
#define RECORD_MACRO_DATA (OUTPUT_DIR "pot-macroscopic.dat")

/* If `RECORD_PHI_DATA` is defined, periodically record electric potential
   data in the file at the path `RECORD_PHI_DATA`. 
#define RECORD_PHI_DATA (OUTPUT_DIR "pot-phi.dat")

*/

/* If `RECORD_V_R_DATA` is defined, periodically record electrons' & ions'
   average radial velocity data at the path `RECORD_V_R_DATA`. */
/*
#define RECORD_V_R_DATA (OUTPUT_DIR "pot-v-r.dat")
*/

/* If `RECORD_N_DATA` is defined, periodically record particle density data
   in the file at the path `RECORD_N_DATA`. */
/*
#define RECORD_N_DATA (OUTPUT_DIR "pot-n.dat")
*/

/* If `RECORD_V_EACH_DATA` is defined, periodically record every particle's
   velocity at the path `RECORD_V_EACH_DATA`. */
/*
#define RECORD_V_EACH_DATA (OUTPUT_DIR "pot-v-each.dat")
*/

/* The opening angle criterion parameter (a.k.a. theta). */
#define OAP 1.0

/* If `DIPOLES` is defined, compute and use electric dipole moments for
   tree cells, instead of just monopole moments.
   If enabled, this should make electric field & potential calculations
   a little more accurate at the expense of some CPU time and memory. */
/*
#define DIPOLES
*/

/* If `SPHERICAL` is defined, the simulation domain is a sphere instead of
   a cube. */
#define SPHERICAL
/*
*/

/* If `PROLATE` is defined, the dust is a prolate spheroid instead of a
   sphere - JH. */
#define PROLATE
/*
*/

/* If `OBLATE` is defined, the dust is an oblate spheroid instead of a
   sphere - JH.
#define OBLATE

*/

/* If `OUT_BOUNCE` is defined, the program bounces particles off the walls
   of the simulation when they're about to leave the simulation area, with
   the value of `OUT_BOUNCE` being the coefficient of restitution.
   Otherwise, the program reinjects particles that depart from the
   simulation area at the opposite wall with random velocity. */
/*
#define OUT_BOUNCE 0.95
*/

/* By default the numeric integrator uses the velocity Verlet algorithm
   to step particles forward. But if `EULER_RICHARDSON` is defined, it
   uses the Euler-Richardson algorithm instead. (N.B. that Euler-Richardson
   is a form of RK2.) */
/*
#define EULER_RICHARDSON
*/

/* Alternatively, if `BORIS` is defined, the numeric integrator steps
   particles forward with the Boris algorithm. */
#define BORIS
/*
*/

/* If `DEBUG_FIELDS` is defined, record the estimated electric fields on
   a subset of about `DEBUG_FIELDS` particles every few ion time steps on
   standard output. */
/*
#define DEBUG_FIELDS 500
*/

/* If `DEBUG_KEPLER` is defined, try to simulate a simple two-body
   electron-ion stable orbit to test the numeric integrator. */
/*
#define DEBUG_KEPLER
*/

/* If `DEBUG_RUTHERFORD` is defined, simulate Rutherford scattering of
   an electron by an ion to test the numeric integrator. */
/*
#define DEBUG_RUTHERFORD
*/

/* If `DEBUG_LANGMUIR` is defined, simulate a Langmuir wave by forming a
   proverbial electron slab at the start of the simulation (by displacing
   a subset of electrons in the x-direction). */
/*
#define DEBUG_LANGMUIR
*/

/* If `POIS_POS_COLLECTION` is defined, positive particles that collide with
   the grain are reinjected, but do not affect the grain's charge. Instead,pot-rin.pbs.bash
   the grain gains positive charge according to a separate, abstract Poisson
   process.
   Important note: if this is defined, `RECORD_MACRO_DATA` must also be
   defined, because the Poisson-based positive particle collection depends
   on the ion temperature that `write_macroscopic_data_to_file` computes. */
/*
#define POIS_POS_COLLECTION
*/

/* If `SCEPTIC_REINJECTION` is defined, and the simulation domain is
   spherical (i.e. if `SPHERICAL` is also defined), try to mimic how Ian
   Hutchinson's SCEPTIC reinjects particles. */
#define SCEPTIC_REINJECTION

/* If `DUST_FIELD_ONLY` is defined, the code does not simulate
   particle-particle interactions; the particles experience an electric
   field and potential from only the central dust grain. */
/*
#define DUST_FIELD_ONLY
*/

/* A path from which to read random bytes for a PRNG seed. */
#define RANDOMNESS_SOURCE "/dev/urandom"

/* Maximum path length allowed when walking particle tree.
   In 3-D, 400 allows for a (big!) particle tree of depth 50. */
#define MAX_PATH_LEN 400

/* Default values for runtime options the user can set as arguments. */
#define DEFAULT_N 150001     /* number of particles */
#ifdef SPHERICAL
	#define DEFAULT_R 654   /* simulation area radius in pixels */
#else
	#define DEFAULT_L 2000   /* simulation area length in pixels */
#endif
#define DEFAULT_DT 1e-12     /* time step length in seconds */
#define DEFAULT_MPP 4e-7     /* pixel size in metres */
#define DEFAULT_RUN_FOR 0    /* no. of time steps (0 => run indefinitely) */
#define DEFAULT_SOFT 4e-10   /* softening parameter in metres */
#define DEFAULT_T_E 11600    /* electron injection temperature in K */
#define DEFAULT_T_I 11600    /* ion injection temperature in K */
#define DEFAULT_X_DRIFT 0.0  /* drift velocity in x-direction */
#define DEFAULT_A 2.5e-6     /* grain radius in metres */
#ifdef GL
	#define DEFAULT_TRACK 0  /* number of particle trajectories to show */
#endif
#define DEFAULT_TIM_LIM 0.0  /* runtime limit (in minutes; 0 => no limit) */
#define DEFAULT_SETTLE 0.0   /* initial settling/equilibration time
	now zero as particles have correct initial spatial dist. - JH */
#ifdef PROLATE
	#define DEFAULT_ASPECT_RATIO 5.0  /* aspect ratio of grain - JH */
#endif
#ifdef OBLATE
	#define DEFAULT_ASPECT_RATIO 0.5  /* aspect ratio of grain - JH */
#endif

/* If a header file hasn't already defined pi, define pi here because C99
   doesn't require a definition of pi. */
#ifndef PI
	#define PI 3.14159265358979323846264338328L
#endif

/* If `POIS_POS_COLLECTION` is defined, `RECORD_MACRO_DATA` had better
   be defined too, for the reason given above. */
#ifdef POIS_POS_COLLECTION
	#ifndef RECORD_MACRO_DATA
		#error "Poisson positive collection enabled but not macroscopic data"
	#endif
#endif

/* If `SCEPTIC_REINJECTION` is defined, `SPHERICAL` should be defined too.
   Trying to mimic SCEPTIC's reinjection method makes no sense for a cubic
   simulation domain. (Also, if using SCEPTIC-like reinjection, define
   `USING_PHI_AT`, since SCEPTIC-like reinjection requires `phi_at`. */
#ifdef SCEPTIC_REINJECTION
	#ifndef SPHERICAL
		#error "SCEPTIC-like reinjection requires spherical domain"
	#endif
	#define USING_PHI_AT
#endif

/* The particle motion-solving numeric integrator can only use one algorithm
   at a time, so the `EULER_RICHARDSON` and `BORIS` definitions are mutually
   exclusive. */
#ifdef EULER_RICHARDSON
	#ifdef BORIS
		#error "can't use both Euler-Richardson & Boris integrators"
	#endif
#endif

/* A spheroid cannot be both prolate and oblate, hence `PROLATE` and `OBLATE`
   are mutually exclusive. */
#ifdef PROLATE
	#ifdef OBLATE
		#error "spheroid cannot be prolate and oblate simultaneously"
	#endif
#endif

/* Recording electric potential information requires definition of the
   `phi_at` function, so if recording electric pot'l information, define
   the preprocessor token that triggers `phi_at`'s definition. */
#ifdef RECORD_PHI_DATA
	#ifndef USING_PHI_AT
		#define USING_PHI_AT
	#endif
#endif

/* Define a macro to compute Pythagorean distance given three orthogonal
   displacements.
   I used to use two nested `hypotl` calls since they're supposed to handle
   underflow & overflow nicely. But then I found that doing so in the inner
   functions raised the runtime by 50%-60% compared to computing distances
   in the less rigorous way used below. And, on reflection, overflows won't
   ever happen because the distances involved are invariably small, and
   underflows shouldn't matter much either. */
#define pythag(x, y, z) sqrt((x)*(x) + (y)*(y) + (z)*(z))

/* Define a `puts`-like macro for stderr. */
#define e_puts(s) fputs(s "\n", stderr)

/* Define a macro to tell whether `x` is a sensible floating-point number. */
#define isnt_weird(x) (isnormal(x) || (fpclassify(x) == FP_ZERO))

#define KB (1.38065e-23)        /* Boltzmann's constant */
#define EPSILON0 (8.85419e-12)  /* electric constant/vacuum permittivity */
#define MU0 (4 * PI * 1e-7)     /* magnetic constant/vacuum permeability */
#define QE (1.6021766e-19)      /* elementary/electron charge */
#define ME (9.109383e-31)       /* electron mass */

/* The different particle species in the simulation plasma.
   The Species enum's values serve as array indices for `CHARGE` & `MASS`. */
enum Species {
	DUST = 0,
	ION = 1,
	ELECTRON = 2
};
const int CHARGE[] = { 0, 1, -1 };
const long double MASS[] = { 1.67e-12, MIME * ME, ME };

/* A particle. (Usually, but not always, a point particle.) */
typedef struct Particle {
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	enum Species species;
} Particle;

/* A node (whether leaf, root, or neither) in a particle tree. */
typedef struct Tree {
	double x1;  /* x-coordinate of left face */
	double y1;  /* y-coordinate of top face */
	double z1;  /* z-coordinate of front face */
	double x2;  /* x-coordinate of right face */
	double y2;  /* y-coordinate of bottom face */
	double z2;  /* z-coordinate of back face */
	unsigned int part;  /* particle contained in this cell/node, if any */
	double cx;  /* centre of charge, x-coordinate */
	double cy;  /* centre of charge, y-coordinate */
	double cz;  /* centre of charge, z-coordinate */
	int q;  /* net charge */
	unsigned int aq;  /* "absolute charge", i.e. sum of charge magnitudes */
	#ifdef DIPOLES
	double px;  /* x-comp. of electric dipole moment about charge centre */
	double py;  /* y-comp. of electric dipole moment about charge centre */
	double pz;  /* z-comp. of electric dipole moment about charge centre */
	#endif
	unsigned int gxlylz;
	unsigned int gxlygz;
	unsigned int lxlylz;
	unsigned int lxlygz;
	unsigned int lxgylz;
	unsigned int lxgygz;
	unsigned int gxgylz;
	unsigned int gxgygz;
} Tree;

/* Settings for a simulation run. */
typedef struct Run {
	unsigned long int N;  /* maximum number of particles */
	#ifdef SPHERICAL
	unsigned int r;  /* window/simulation area radius (in pixels) */
	#else
	unsigned int L;  /* window/simulation area length (in pixels) */
	#endif
	long double dt;             /* time step to use (in seconds) */
	double mpp;                 /* metres per pixel, i.e. pixel size */
	unsigned long int run_for;  /* number of time steps to run for */
	double soft;                /* softening parameter */
	double x_drift;             /* relative x-direction drift velocity */
	double a;                   /* grain radius (in metres) */
	#ifdef GL
	unsigned short track;  /* number of particle trajectories to display */
	#endif
	double settle;   /* settling time before enabling grain absorption */
	char* sw_path;   /* path to which to write final simulation state */
	char* sr_path;   /* path from which to read initial simulation state */
	double tim_lim;  /* wall-time limit */
	double Te;       /* initial & reinjection electron temperature */
	double Ti;       /* initial & reinjection ion temperature */
	#if defined(PROLATE) || defined(OBLATE)
	double A;	 /* aspect ratio of grain - JH*/
	#endif
} Run;

/* A single process' state and interprocess communication data. */
typedef struct Proc {
	int rank;  /* MPI proc. rank, also used as unique MPI proc. ID */
	int num_procs;        /* total number of processes running */
	Particle* parts;      /* particle array */
	Particle* new_parts;  /* temporary storage for updated particle array */
	PRNG ps;              /* this process' PRNG state */
	unsigned int** pits;  /* particle indices to step */
	MPI_Comm work_comm;   /* MPI communicator for workers to communicate */
	long double epe[3];   /* each species' total electrostatic pot'l energy */
	FILE* fp_macro;       /* macroscopic summary file ptr. (rank 0 only) */
	int dust_grain_q;     /* net elementary charges on dust grain */
	long double dust_grain_m;  /* dust grain mass */
	long double* a_x1;  /* mid-time step acceleration, x component */
	long double* a_y1;  /* mid-time step acceleration, y component */
	long double* a_z1;  /* mid-time step acceleration, z component */
	Tree* root;             /* root node of particle tree */
	unsigned int tr_alloc;  /* allocated number of nodes for p'cle tree */
	unsigned int tss[3];    /* time steps to use for each species */
	#ifdef GL
	float** traj_x;  /* particle trajectories' x components (in pixels) */
	float** traj_y;  /* particle trajectories' y components (in pixels) */
	float** traj_z;  /* particle trajectories' z components (in pixels) */
	unsigned char* traj_reset;  /* which trajectories needs resetting */
	#endif
	double* prev_x;  /* particles' x-positions at previous time step */
	double* prev_y;  /* particles' y-positions at previous time step */
	double* prev_z;  /* particles' z-positions at previous time step */
	unsigned int* starts;  /* start of p'cles moved this step (rank 0 only) */
	#ifdef POIS_POS_COLLECTION
	double Ti;              /* cached ion temperature */
	double pos_rate_const;  /* constant part of ion collection rate */
	#endif
	unsigned long int t_steps;  /* number of time steps elapsed */
	unsigned int num_workers;   /* number of worker processes running */
	unsigned int tot_absor_e;   /* total electrons absorbed (rank 0 only) */
	#ifdef SCEPTIC_REINJECTION
	double chi_b;  /* normalized electric potential measured by this proc. */
	#endif
} Proc;

#include "coord.c"  /* easier to just suck coord.c in here, trust me */

#ifdef SPHERICAL
#ifdef SCEPTIC_REINJECTION

/* Solve the orbit integral for the angle `alpha` between an initially
   distant particle's initial velocity and the position where it would enter
   the simulation domain, as Ian Hutchinson's SCEPTIC does. The arguments
   `s`, `chi_b` & `u2` have the same meaning as in Hutchinson's work (the
   ratio of the impact parameter `b` to the simulation domain radius, the
   normalized edge electric potential, and the square of the normalized
   velocity at infinity respectively); `Rol` denotes the ratio of the domain
   radius to the relevant (Debye) shielding length. The integration
   algorithm is a crudely adaptive Simpson's rule; usually this works well
   because the orbit integrand is usually a smooth, well-defined curve.

   But not always! Certain parameter values represent particles which
   would encounter effective potential barriers on their way to the
   simulation domain. For such particles the orbit integrand is imaginary
   at some radii and `alpha` is undefined and physically meaningless.
   In such situations this function is liable to return NaN; it is the
   caller's responsibility to detect that (or to avoid the entire situation
   by never calling this function with parameter values that contradict its
   implicit OML model).

   (This function integrates a simpler, less exact integrand than SCEPTIC.
    When computing the integral, SCEPTIC incorporates a more complicated
    version of the Debye-Hueckel potential profile which accounts for ion
    depletion, whereas this function uses the ordinary D-H profile,
    neglecting ion depletion. This makes little difference for small
    grains because for small grains the profiles differ only very slightly
    in absolute terms.) */
double integrate_for_alpha(double s, double chi_b, double u2, double Rol)
{
	const double chu2 = chi_b / u2;
	const double max_lin_dev = 0.1;
	const double s2 = s * s;
	const double y_one = s / sqrt(1.0 - chu2 - s2);  /* `xi` = 1 integ'nd */

	double alpha = 0.0;
	double xi[3] = { 0.0, 0.5, 1.0 };  /* current `xi` interval */
	double y[3] = { s, -9.9, y_one };  /* integrand values at `xi` */

	/* Handle the special case of head-on motion towards the
	   simulation domain's centre (i.e. zero impact parameter). */
	if (s == 0.0) {
		return 0.0;
	}

	do {
		/* Try defining the current integration interval as starting
		   at `xi[0]` and stretching all the way to 1. */
		xi[2] = 1.0;
		y[2] = y_one;
		xi[1] = (xi[0] + xi[2]) / 2.0;
		y[1] = 1.0 - (chu2 * xi[1] * exp(Rol * (1.0 - (1.0 / xi[1]))));
		y[1] = 1.0 / sqrt((y[1] / s2) - (xi[1] * xi[1]));

		/* If the integrand (`y`) appears to be linear over this interval,
		   excellent -- don't bother entering the loop. Otherwise, enter
		   the loop and iteratively shrink the interval until it's small
		   enough that `y` is close to linear throughout. */
		while (fabs(y[2] - y[1] - (y[1] - y[0])) > max_lin_dev) {
			/* `y` deviates substantially from linearity over the
			   current `xi` interval, so shrink the interval to half
			   its size and recompute `y` at the new mid-point and end
		       point, before checking again for near-linearity. */
			xi[2] = xi[1];
			y[2] = y[1];
			xi[1] = (xi[0] + xi[2]) / 2.0;
			y[1] = 1.0 - (chu2 * xi[1] * exp(Rol * (1.0 - (1.0 / xi[1]))));
			y[1] = 1.0 / sqrt((y[1] / s2) - (xi[1] * xi[1]));
		}

		/* `y` is approximately linear over the interval `xi`,
		   so apply Simpsons's rule over `xi`, incrementing `alpha`.
		   Then shift the interval's start to where the interval's end
		   currently is. */
		alpha += (xi[2] - xi[0]) * (y[0] + (4.0 * y[1]) + y[2]) / 6.0;
		xi[0] = xi[2];
		y[0] = y[2];
	} while (xi[2] < 1.0);

	/* Display a warning about this integral if it evaluated to NaN
	   (implying bad parameters for this function). */
	if (isnan(alpha)) {
		fprintf(stderr, "integrate_for_alpha got NaN: %.8g %.9g %g %g\n",
		        s, u2, chi_b, Rol);
	}

	return alpha;
}

/* Reinject a particle at the simulation domain edge using a reinjection
   algorithm like that in Ian Hutchinson's SCEPTIC. */
void sceptically_reinject(Particle* part, Run ru, PRNG* state, double chi_b)
{
	double a[2];   /* general-purpose rotation axis in the x-y plane */
	double alpha;  /* Hutchinson's alpha angle */
	double b;      /* normalized impact parameter */
	double c;      /* cosine of Hutchinson's theta angle */
	double co_z;   /* cosine of angle zeta, used below */
	double phi;    /* azimuthal angle in spherical coordinate system */
	double psi;    /* auxiliary rotation angle */
	double q;      /* u.r.v. quantile for sampling from CDFs */
	double Rol;    /* ratio of simul'n domain radius to shielding length */
	double U;      /* normalized drift velocity */
	double u_mag;  /* initial speed at infinity */
	double u2;     /* square of initial speed at infinity */
	double u[3];   /* particle's velocity at infinity (Cart. coords.) */
	double v[3];   /* reinj. particle's veloc. at boundary (sph. coords.) */
	double v_n;    /* characteristic normalization velocity */

	/* Compute the normalization velocity (square root of 2 k_B T / m).
	   Also, if `part` represents an electron, renormalize `chi_b` and
	   flip its sign so that it makes sense for an electron. */
	v_n = 2.0 * KB / MASS[part->species];
	if (part->species == ELECTRON) {
		v_n *= ru.Te;
		chi_b *= -ru.Ti / ru.Te;
	} else {
		v_n *= ru.Ti;
	}
	v_n = sqrt(v_n);

	/* Choose a random initial speed at infinity `u_mag`, ensuring its
	   square is no less than `chi_b` (otherwise the implied impact
	   parameter `b`, to be sampled below, would be imaginary --
	   corresponding to a particle which lacks the energy to reach the
	   simulation domain under the OML assumption). Then choose a random
	   cosine `c` of theta, theta being the angle between the drift
	   velocity and the initial velocity at infinity. */
	if (ru.x_drift) {
		U = ru.x_drift / v_n;  /* constant, so compute before the loop */
		do {
			u_mag = u_icdf(U, chi_b, rnd(state, 1.0));
			u2 = u_mag * u_mag;
		} while (u2 < chi_b);
		q = rnd(state, 1.0);
		c = log(q + (1.0 - q) * exp(-4.0 * u_mag * U)) / (2.0 * u_mag * U);
		c += 1.0;
	} else {
		do {
			u_mag = u_no_flow_icdf(chi_b, rnd(state, 1.0));
			u2 = u_mag * u_mag;
		} while (u2 < chi_b);
		c = rnd(state, 2.0) - 1.0;
	}

	/* With `u_mag` & `c` now in hand, it's time to randomly decide
	   the three components of the initial velocity `u`.
	   First, initialize `u` to be entirely in the z-direction by setting
	   the z-component to `u_mag`. Then rotate `u` about the y-axis
	   through an angle of pi/2 minus arccos `c` radians, to bring it
	   onto the cone of possible velocities defined as having angle theta
	   about the x-axis (the drift direction). Finally, to randomize `u`'s
	   precise direction, rotate it again, this time about the x-axis,
	   through a uniformly distributed random angle (call it `psi`). */
	u[0] = 0.0;
	u[1] = 0.0;
	u[2] = u_mag;
	a[0] = 0.0;
	a[1] = 1.0;
	rotate_vec_about_xy_dir(a, c, sqrt(1.0 - (c * c)),
	                        &(u[0]), &(u[1]), &(u[2]), 0.0);
	a[0] = 1.0;
	a[1] = 0.0;
	psi = rnd(state, 2.0 * PI);
	rotate_vec_about_xy_dir(a, sin(psi), cos(psi),
	                        &(u[0]), &(u[1]), &(u[2]), 0.0);

	/* Choose a random impact parameter `b` at infinity.
	   Bad things would happen here were `chi_b` > `u2`, but the code
	   sampling `u_mag` above should ensure `chi_b` <= `u2`. */
	b = sqrt(rnd(state, 1.0 - (chi_b / u2)));

	/* Compute Hutchinson's alpha angle from `b`, `chi_b` and `u2`. */
	Rol = QE * QE * (ru.N - 1) / 2.0;
	Rol /= (4.0 * PI / 3.0) * ru.mpp * ru.r * EPSILON0 * KB * ru.Te;
	Rol = sqrt(Rol);
	alpha = integrate_for_alpha(b, chi_b, u2, Rol);
	if (isnan(alpha)) {
		/* `alpha` isn't well defined; given the sampled `b`, `chi_b`,
		   `u2` and potential profile used by `integrate_for_alpha`,
		   this particle would apparently have encountered an effective
		   potential barrier on its way to the simulation domain.
		   So try to reinject it again. Since this problem is unlikely
		   to occur many times in a row, this function reattempts
		   reinjection by calling itself; since reinjection attempts
		   usually succeed, this should never overflow the stack. */
		if (part->species == ELECTRON) {
			/* Undo the electron-specific renormalization of `chi_b` done
			   at the start of this function before passing it to this
			   function again. */
			chi_b /= -ru.Ti / ru.Te;
		}
		sceptically_reinject(part, ru, state, chi_b);
		return;
	}

	/* Time to generate the particle's position at the simulation edge
	   based on its `u` & `b`. Begin by putting the particle at the edge at
	   an angle `alpha` from the z-axis, since that's easy. */
	phi = rnd(state, 2.0 * PI);
	part->x = ru.mpp * ru.r * sin(alpha) * cos(phi);
	part->y = ru.mpp * ru.r * sin(alpha) * sin(phi);
	part->z = ru.mpp * ru.r * cos(alpha);

	/* But the particle still needs to be moved so that it's at an
	   an angle `alpha` from the vector `u`, not from the z-axis! So
	   derive a rotation axis `a` by taking the cross product of
	   `u` and a unit z-axis vector. Then normalize `a`. */
	a[0] = u[1];
	a[1] = -u[0];
	a[0] /= sqrt((u[0] * u[0]) + (u[1] * u[1]));
	a[1] /= sqrt((u[0] * u[0]) + (u[1] * u[1]));

	/* Deduce (the cosine of) the angle zeta through which to rotate the
	   particle about `a`, then rotate the particle through zeta. Also,
	   increment each of the particle's coordinates by the simulation
	   radius to bring the particle inside the actual simulation domain
	   (which has an origin different to the ideal mathematical origin of
	   (0, 0, 0)). */
	co_z = u[2] / u_mag;
	rotate_vec_about_xy_dir(a, sin(acos(co_z)), co_z, &(part->x),
	                        &(part->y), &(part->z), ru.mpp * ru.r);

	/* Now the particle's velocity must be set. As before, set the velocity
	   based on the pretence that `u` is along the z-axis. In this case the
	   reinjection velocity has no azimuthal component (by c.o.a.m.), which
	   beautifully simplifies the velocity calculation; the zenith comp't
	   follows directly from angular momentum conservation and the radial
	   comp't follows from Pythagoras' theorem and energy conservation. */
	v[1] = 0.0;
	v[2] = u_mag * b;
	v[0] = -sqrt(u2 - chi_b - (v[2] * v[2]));
	velocity_in_cartesian_coords(v, phi, alpha,
		                         &(part->vx), &(part->vy), &(part->vz));

	/* Like the position, the velocity vector must also be rotated
	   through zeta about `a`, since `u` is (generally) not actually
	   parallel to the z-axis. */
	rotate_vec_about_xy_dir(a, sin(acos(co_z)), co_z, &(part->vx),
	                        &(part->vy), &(part->vz), 0.0);

	/* The calculations in this function have used normalized units.
	   Unnormalize the velocity so that it's back in the expected units
	   (metres per second). */
	part->vx *= v_n;
	part->vy *= v_n;
	part->vz *= v_n;
}

#endif /* #ifdef SCEPTIC_REINJECTION */
#endif /* #ifdef SPHERICAL */

/* Inject/reinject a particle with a random position & velocity. */
void inject(Particle* part, Run ru, bool reinjection, Proc* pro)
{
	/* Standard deviation of the velocity component distribution. */
	double sigma;

	/* This process' PRNG's state. */
	PRNG* state = &(pro->ps);

	#ifdef SPHERICAL
	/* Spherical coordinates at which to reinject a particle. */
	double phi;
	double r;
	double theta;
	/* Spherical components of the (re)injection velocity. */
	double v[3];  /* 0 is v_r, 1 is v_phi, and 2 is v_theta */
	#endif

	/* Deduce the appropriate standard deviation for this particle's
	   velocity component sampling distribution from the temperature and
	   particle's mass. */
	if (part->species == ION) {
		sigma = sqrt(KB * ru.Ti / MASS[ION]);
	} else {
		sigma = sqrt(KB * ru.Te / MASS[ELECTRON]);
	}
	#ifndef SCEPTIC_REINJECTION
	if (!reinjection) {
		/* When initializing a particle at a simulation's start, lower the
		   effective temperature, because at equilibrium the simulation
		   invariably settles to a temperature that's only 65%-75% of the
		   injection temperature. */
		sigma *= 0.85;  /* 0.85 is approx. sqrt(0.72) */
	}
	#endif

	#ifdef SPHERICAL

	/* Decide how far from the centre to (re)inject the particle, based on
	   whether this is an initialization or a reinjection. */
	if (!reinjection) {
		/* Put the particle at a random distance from the centre (although
		   not so close to the centre that it's within the grain).
		   The **CUBE** root transformation should ensure a uniform
		   spatial distribution. */
		/* Bug fix - cube root transform required, not square root, so
		   settling time is now redundant! JH (spotted by Drew) */
		do {
			r = ru.mpp * cbrt(rnd(state, ru.r * ru.r * ru.r));
		} while (r <= ru.a);
	} else {
		/* Put the particle at the simulation's outer boundary. */
		r = ru.mpp * ru.r;
	}

	#ifdef SCEPTIC_REINJECTION
	if (reinjection) {
		/* Reinject this particle with a SCEPTIC-like algorithm. */
		sceptically_reinject(part, ru, state, pro->chi_b);
		return;  /* position & velocity now set, nothing left to do here */
	} else {
		/* Inject this particle at a random position in the domain. */
		phi = rnd(state, 2.0 * PI);
		theta = acos(rnd(state, 2.0) - 1.0);
		part->x = (ru.mpp * ru.r) + (r * sin(theta) * cos(phi));
		part->y = (ru.mpp * ru.r) + (r * sin(theta) * sin(phi));
		part->z = (ru.mpp * ru.r) + (r * cos(theta));
	}
	#else
	/* (Re)place the particle at a random position at distance `r`
	   from the simulation's centre. */
	phi = rnd(state, 2.0 * PI);
	theta = acos(rnd(state, 2.0) - 1.0);
	part->x = (ru.mpp * ru.r) + (r * sin(theta) * cos(phi));
	part->y = (ru.mpp * ru.r) + (r * sin(theta) * sin(phi));
	part->z = (ru.mpp * ru.r) + (r * cos(theta));
	#endif

	/* Set random velocity components, then add the drift velocity. */
	if (!reinjection) {
		nrv_pair(state, sigma, &(part->vx), &(part->vy));
		nrv_pair(state, sigma, &(part->vz), &(part->vz));
	} else {
		nrv_pair(state, sigma, &(v[0]), &(v[1]));
		nrv_pair(state, sigma, &(v[2]), &(v[2]));
		if (v[0] > 0.0) {
			/* No point injecting a particle at the outer boundary
			   with positive radial velocity. */
			v[0] = -v[0];
		}
		#ifdef DEBUG_RUTHERFORD
		if (part->species == ELECTRON) {
			/* For Rutherford scattering testing it's desirable to have the
			   electron pass reasonably close to the ion, so boost the radial
			   velocity component and shrink the non-radial components. */
			v[0] *= 1.5;
			v[1] /= 8.0;
			v[2] /= 8.0;
		}
		#endif
		/* Translate the injection velocity in spherical coordinates into
		   the desired Cartesian coordinates. Luckily I already had this
		   worked out for John and me's spherical absorber problem! */
		velocity_in_cartesian_coords(v, phi, theta,
		                             &(part->vx), &(part->vy), &(part->vz));
	}
	part->vx += ru.x_drift;

	#else

	/* Set random velocity components. */
	nrv_pair(state, sigma, &(part->vx), &(part->vy));
	nrv_pair(state, sigma, &(part->vz), &(part->vz));

	if (!reinjection) {
		/* At the simulation's start, place particles uniformly randomly. */
		part->x = ru.mpp * rnd(state, ru.L);
		part->y = ru.mpp * rnd(state, ru.L);
		part->z = ru.mpp * rnd(state, ru.L);
	} else {
		/* Reinject a particle that leaves the simulation. How it gets
		   reinjected depends on whether it left through a wall or was
		   instead absorbed by the dust grain. */
		if ((part->x >= 0.0) && (part->x <= (ru.mpp * ru.L))
		    && (part->y >= 0.0) && (part->y <= (ru.mpp * ru.L))
		    && (part->z >= 0.0) && (part->z <= (ru.mpp * ru.L))) {
			/* The particle needs to be reinjected but is inside the
			   simulation area. This implies the dust grain absorbed it.
			   Reinject the particle at a random wall. FIXME: check that
			   the velocity after reinjection isn't into the wall.
			   Yes, this PRNG call will very occasionally return 6, and
			   this will introduce a tiny bias in favour of reinjection at
			   a particular wall. But it is surely negligible. */
			part->x = ru.mpp * rnd(state, ru.L);
			part->y = ru.mpp * rnd(state, ru.L);
			part->z = ru.mpp * rnd(state, ru.L);
			switch ((int) rnd(state, 6.0)) {
			case 0:  part->x = 0.0;           break;
			case 1:  part->x = ru.mpp * ru.L; break;
			case 2:  part->y = 0.0;           break;
			case 3:  part->y = ru.mpp * ru.L; break;
			case 4:  part->z = 0.0;           break;
			default: part->z = ru.mpp * ru.L; break;
			}
		} else {
			/* The grain left by striking a wall. Reinject it at the
			   opposite wall with its other coordinates intact. */
			if (part->x < 0.0) {
				part->x = ru.mpp * ru.L;
				if (part->vx > 0.0) {
					part->vx = -part->vx;
				}
			} else if (part->x > (ru.mpp * ru.L)) {
				part->x = 0.0;
				if (part->vx < 0.0) {
					part->vx = -part->vx;
				}
			}
			if (part->y < 0.0) {
				part->y = ru.mpp * ru.L;
				if (part->vy > 0.0) {
					part->vy = -part->vy;
				}
			} else if (part->y > (ru.mpp * ru.L)) {
				part->y = 0.0;
				if (part->vy < 0.0) {
					part->vy = -part->vy;
				}
			}
			if (part->z < 0.0) {
				part->z = ru.mpp * ru.L;
				if (part->vz > 0.0) {
					part->vz = -part->vz;
				}
			} else if (part->z > (ru.mpp * ru.L)) {
				part->z = 0.0;
				if (part->vz < 0.0) {
					part->vz = -part->vz;
				}
			}
		}
	}

	/* Add the drift velocity. */
	part->vx += ru.x_drift;

	#endif
}

#ifdef DEBUG_RUTHERFORD
double rutherford_impact_parameter(Particle* ion, Particle* electron)
{
	double s[3];     /* vector from electron to ion */
	double v_dot_s;  /* dot product of electron velocity and `s` */
	double v2;       /* square of electron's speed */

	s[0] = ion->x - electron->x;
	s[1] = ion->y - electron->y;
	s[2] = ion->z - electron->z;
	v_dot_s = (electron->vx * s[0]) + (electron->vy * s[1])
	          + (electron->vz * s[2]);
	v2 = (electron->vx * electron->vx) + (electron->vy * electron->vy)
	     + (electron->vz * electron->vz);
	return sqrt(((s[0] * s[0]) + (s[1] * s[1]) + (s[2] * s[2]))
	            - (v_dot_s * v_dot_s / v2));
}

/* Estimate the electron's squared speed infinitely far from the ion
   based on its current squared speed and distance from the ion.
   pot does not use this function, because it introduces a systematic
   deviation from the Rutherford scattering prediction because the
   scattering angle is still estimated in the imprecise way (based on
   the electron's initial & final velocities at the simulation edge). */
double true_v0_squared(Particle* ion, Particle* electron)
{
	double r;
	double v2;

	r = pythag(electron->x - ion->x, electron->y - ion->y,
	           electron->z - ion->z);
	v2 = (electron->vx * electron->vx) + (electron->vy * electron->vy)
	     + (electron->vz * electron->vz);
	return v2 + (QE * QE / (2.0 * PI * EPSILON0 * MASS[ELECTRON] * r));
}
#endif

/* Initialize the particles by setting random initial positions &
   velocities and assigning the particles appropriate charges & masses. */
void init_particles(Proc* pro, Run ru)
{
	unsigned int i;
	unsigned int N_i;
	unsigned int p_idx;
	Particle* parts = pro->parts;
	unsigned int proc;

	/* Initialize the dust grain. */
	#ifdef SPHERICAL
	parts[0].x = ru.r * ru.mpp;
	#else
	parts[0].x = ru.L * ru.mpp / 2.0;
	#endif
	parts[0].y = parts[0].x;
	parts[0].z = parts[0].x;
	parts[0].vx = 0.0;
	parts[0].vy = 0.0;
	parts[0].vz = 0.0;
	parts[0].species = DUST;

	/* To begin with, set every non-dust particle as an electron. */
	for (p_idx = 1; p_idx < ru.N; p_idx++) {
		parts[p_idx].species = ELECTRON;
	}

	/* Now, change electrons to ions one by one until there're the desired
	   number of ions. (As many ions as electrons if possible, otherwise
	   one fewer ion than there are electrons.) Moreover, do so in a way
	   that puts all of a process' ions in one contiguous block, and all of
	   of a process' electrons in one contiguous block. */
	i = 0;
	N_i = (ru.N - 1) / 2;
	while (N_i) {
		for (proc = 0; proc < pro->num_workers; proc++) {
			p_idx = pro->pits[proc][0] + i;
			if (!p_idx) {
				/* Don't set this particle as an ion -- it's the grain! */
				continue;
			}
			parts[p_idx].species = ION;
			N_i--;
			if (!N_i) {
				/* No more electrons need be reset as ions. */
				break;
			}
		}
		i++;
	}

	/* Inject the particles. */
	for (p_idx = 1; p_idx < ru.N; p_idx++) {
		inject(&(parts[p_idx]), ru, false, pro);
	}

	#ifdef DEBUG_KEPLER
		/* Position the electron & ion for the 2-body Kepler problem test. */
		if (ru.N != 3) {
			e_puts("Number of particles must be set to 3 "
			       "for two-body debugging.");
			abort();
		}
		#ifdef SPHERICAL
		parts[2].species = ELECTRON;
		parts[2].x = ((ru.r + 0.5) * ru.mpp);
		parts[2].y = ((ru.r - 0.5) * ru.mpp);
		parts[2].z = 0.5 * ru.r * ru.mpp;
		parts[2].vx = 0.0;
		parts[2].vy = QE * sqrt(CHARGE[ION] / (4.0 * PI * EPSILON0 * MASS[ELECTRON] * 0.5 * ru.r * ru.mpp));
		parts[2].vz = 0.0;
		parts[1].species = ION;
		parts[1].x = parts[2].x;
		parts[1].y = parts[2].y;
		parts[1].z = ru.r * ru.mpp;
		parts[1].vx = 0.0;
		parts[1].vy = -(MASS[ELECTRON] / MASS[ION]) * parts[2].vy;
		parts[1].vz = 0.0;
		#else
		#error "for Kepler tests, use a spherical domain"
		#endif
	#endif

	#ifdef DEBUG_RUTHERFORD
		/* Position the ion for Rutherford scattering tests, and inject
		   the electron somewhere at the simulation edge. */
		if (ru.N != 3) {
			e_puts("Number of particles must be set to 3 "
			       "for Rutherford scattering debugging.");
			abort();
		}
		parts[1].species = ION;
		#ifdef SPHERICAL
		parts[1].x = ((ru.r + 0.5) * ru.mpp);
		parts[1].y = ((ru.r - 0.5) * ru.mpp);
		parts[1].z = ru.r * ru.mpp;
		#else
		#error "for Rutherford scattering, use a spherical domain"
		#endif
		parts[1].vx = 0.0;
		parts[1].vy = 0.0;
		parts[1].vz = 0.0;
		parts[2].species = ELECTRON;
		inject(&(parts[2]), ru, true, pro);
		printf("%g %g %g %g ",
		       rutherford_impact_parameter(&(parts[1]), &(parts[2])),
		       parts[2].vx, parts[2].vy, parts[2].vz);
	#endif
}

#ifdef FPE

/* When `FPE` is defined, `main` unmasks some floating point exceptions
   to allow detection of floating point bugs. However, annoyingly, some
   OpenGL calls like to trigger irrelevant FPEs. This function disables
   the FPEs `main` enabled so some OpenGL calls can be made safely.
   N.B.: this isn't ANSI C! */
void disable_fpes(void)
{
	fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
}

/* Turn most FP exceptions back on again by disabling exception masking. */
void enable_fpes(void)
{
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
}

#endif

/* If the GUI's enabled, clear the window. */
#ifdef GL
void blank_window(void)
{
	#ifdef FPE
	disable_fpes();
	#endif
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	#ifdef FPE
	enable_fpes();
	#endif
}
#endif

/* Clear the window and draw the particles at their positions. */
#ifdef GL
void draw_particles(GLUquadric* qr, const Proc* pro, Run ru)
{
	unsigned int i;
	unsigned int j;
	unsigned int p_idx;
	float x;
	float y;
	float z;

	#ifdef FPE
	disable_fpes();
	#endif

	for (i = 1; i < ru.N; i++) {

		glPushMatrix();
		glTranslatef(pro->parts[i].x / ru.mpp,
		             pro->parts[i].y / ru.mpp,
		             pro->parts[i].z / ru.mpp);
		if (pro->parts[i].species == ION) {
			glColor3f(1.0, 0.0, 0.0);
		} else {
			glColor3f(0.2, 0.2, 1.0);
		}
		gluSphere(qr, 3.0, 8, 8);
		glPopMatrix();

		/* If particle `i` is an electron, and there are not many particles,
		   and no particle trajectories are being tracked, draw a translucent
		   blue line between the centre of particle `i` and the dust grain's
		   centre. */
		if ((!ru.track) && (ru.N < 500) && (pro->parts[i].species == ELECTRON)) {
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glColor4f(0.0, 0.0, 0.8, 0.15);
			glBegin(GL_LINES);
			glVertex3f(pro->parts[i].x / ru.mpp,
			           pro->parts[i].y / ru.mpp,
			           pro->parts[i].z / ru.mpp);
			glVertex3f(pro->parts[0].x / ru.mpp,
			           pro->parts[0].y / ru.mpp,
			           pro->parts[0].z / ru.mpp);
			glEnd();
			glDisable(GL_BLEND);
		}

		glLineWidth(1.0);
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_LINE_STRIP);
		glVertex3f(0.0, 0.0, 0.0);
		#ifdef SPHERICAL
		glVertex3f(2.0 * ru.r, 0.0, 0.0);
		#else
		glVertex3f(ru.L, 0.0, 0.0);
		#endif
		glEnd();
		glColor3f(0.0, 1.0, 0.0);
		glBegin(GL_LINE_STRIP);
		glVertex3f(0.0, 0.0, 0.0);
		#ifdef SPHERICAL
		glVertex3f(0.0, 2.0 * ru.r, 0.0);
		#else
		glVertex3f(0.0, ru.L, 0.0);
		#endif
		glEnd();
		glColor3f(0.0, 0.0, 1.0);
		glBegin(GL_LINE_STRIP);
		glVertex3f(0.0, 0.0, 0.0);
		#ifdef SPHERICAL
		glVertex3f(0.0, 0.0, 2.0 * ru.r);
		#else
		glVertex3f(0.0, 0.0, ru.L);
		#endif
		glEnd();
		glLineWidth(1.0);

	}

	/* Draw recorded particle trajectories as series of straight lines. */
	for (i = 0; i < ru.track; i++) {

		p_idx = (i * ((ru.N - 1) / ru.track)) + 1;

		#ifdef TRAJ_ALPHA
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if (pro->parts[p_idx].species == ION) {
			glColor4f(1.0, 0.3, 0.3, TRAJ_ALPHA);
		} else {
			glColor4f(0.4, 0.4, 1.0, TRAJ_ALPHA);
		}
		#else
		if (pro->parts[p_idx].species == ION) {
			glColor3f(1.0, 0.3, 0.3);
		} else {
			glColor3f(0.4, 0.4, 1.0);
		}
		#endif

		glBegin(GL_LINE_STRIP);
		for (j = 1; j <= TRAJ_LEN; j++) {
			/* Notice that trajectories are plotted backwards, i.e.
			   starting with the most recent point and drawing lines
			   connecting progressively older points until the list of
			   points runs out. */
			x = pro->traj_x[i][TRAJ_LEN - j];
			y = pro->traj_y[i][TRAJ_LEN - j];
			z = pro->traj_z[i][TRAJ_LEN - j];
			if ((x == 0.0) && (y == 0.0) && (z == 0.0)) {
				break;
			}
			glVertex3f(x, y, z);
		}
		glEnd();

		#ifdef TRAJ_ALPHA
		glDisable(GL_BLEND);
		#endif

	}

	/* Enable translucency for the imminent grain & boundary rendering. */
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/* Render the dust grain in green (opaque & bright green if
	   particle collection's enabled, translucent green if not). */
	#ifdef DISABLE_GRAIN
	glColor4f(0.0, 1.0, 0.0, 0.13);
	#endif
	glPushMatrix();
	glTranslatef(pro->parts[0].x / ru.mpp,
	             pro->parts[0].y / ru.mpp,
	             pro->parts[0].z / ru.mpp);
	#ifndef DISABLE_GRAIN
	glColor3f(0.0, 1.0, 0.0);
	#endif
	#if defined(PROLATE) || defined(OBLATE)
	/* If `PROLATE`/`OBLATE` is defined, plot spheroid, else sphere - JH */
	glScalef(pow( ru.A, 2.0/3.0), pow( ru.A, -1.0/3.0), pow( ru.A, -1.0/3.0));
	#endif
	gluSphere(qr, ru.a / ru.mpp, 99, 99);
	glPopMatrix();
	#ifdef DISABLE_GRAIN
	#endif

	#ifdef SPHERICAL
	/* Render the simulation domain's boundary as a translucent sphere. */
	glColor4f(1.0, 1.0, 1.0, 0.15);
	glPushMatrix();
	glTranslatef(ru.r, ru.r, ru.r);
	gluSphere(qr, ru.r, 299, 299);
	glPopMatrix();
	#else
	/* Render the simulation domain's boundary as a translucent cube with
	   one open face. */
	glDisable(GL_CULL_FACE);  /* some cube faces go missing otherwise */
	glPushMatrix();
	glTranslatef(0.0, 0.0, ru.L);
	glColor4f(0.05, 0.05, 0.05, 0.4);
	glRectf(0.0, 0.0, ru.L, ru.L);
	glPopMatrix();
	glPushMatrix();
	glRotatef(270.0, 0.0, 1.0, 0.0);
	glColor4f(0.15, 0.15, 0.15, 0.4);
	glRectf(0.0, 0.0, ru.L, ru.L);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(ru.L, 0.0, 0.0);
	glRotatef(270.0, 0.0, 1.0, 0.0);
	glColor4f(0.15, 0.15, 0.15, 0.4);
	glRectf(0.0, 0.0, ru.L, ru.L);
	glPopMatrix();
	glPushMatrix();
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glColor4f(0.1, 0.1, 0.1, 0.4);
	glRectf(0.0, 0.0, ru.L, ru.L);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.0, ru.L, 0.0);
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glColor4f(0.1, 0.1, 0.1, 0.4);
	glRectf(0.0, 0.0, ru.L, ru.L);
	glPopMatrix();
	glEnable(GL_CULL_FACE);  /* turn culling back on so spheres look good */
	#endif

	glDisable(GL_BLEND);  /* turn transparency back off */

	#ifdef FPE
	enable_fpes();
	#endif
}
#endif

void init_node_children(Tree* node_list, unsigned int node_no)
{
	Tree* children[8];
	unsigned int i;
	double mid_x;
	double mid_y;
	double mid_z;
	Tree* tr = &(node_list[node_no]);

	children[0] = &(node_list[tr->gxlylz]);
	children[1] = &(node_list[tr->gxlygz]);
	children[2] = &(node_list[tr->lxlylz]);
	children[3] = &(node_list[tr->lxlygz]);
	children[4] = &(node_list[tr->lxgylz]);
	children[5] = &(node_list[tr->lxgygz]);
	children[6] = &(node_list[tr->gxgylz]);
	children[7] = &(node_list[tr->gxgygz]);

	for (i = 0; i < 8; i++) {
		children[i]->x1 = tr->x1;
		children[i]->x2 = tr->x2;
		children[i]->y1 = tr->y1;
		children[i]->y2 = tr->y2;
		children[i]->z1 = tr->z1;
		children[i]->z2 = tr->z2;
		children[i]->part = UINT_MAX;
		children[i]->gxlylz = 0;
		children[i]->gxlygz = 0;
		children[i]->lxlylz = 0;
		children[i]->lxlygz = 0;
		children[i]->lxgylz = 0;
		children[i]->lxgygz = 0;
		children[i]->gxgylz = 0;
		children[i]->gxgygz = 0;
	}

	mid_x = (tr->x1 + tr->x2) / 2;
	mid_y = (tr->y1 + tr->y2) / 2;
	mid_z = (tr->z1 + tr->z2) / 2;
	children[0]->x1 = mid_x;
	children[0]->y2 = mid_y;
	children[0]->z2 = mid_z;
	children[1]->x1 = mid_x;
	children[1]->y2 = mid_y;
	children[1]->z1 = mid_z;
	children[2]->x2 = mid_x;
	children[2]->y2 = mid_y;
	children[2]->z2 = mid_z;
	children[3]->x2 = mid_x;
	children[3]->y2 = mid_y;
	children[3]->z1 = mid_z;
	children[4]->x2 = mid_x;
	children[4]->y1 = mid_y;
	children[4]->z2 = mid_z;
	children[5]->x2 = mid_x;
	children[5]->y1 = mid_y;
	children[5]->z1 = mid_z;
	children[6]->x1 = mid_x;
	children[6]->y1 = mid_y;
	children[6]->z2 = mid_z;
	children[7]->x1 = mid_x;
	children[7]->y1 = mid_y;
	children[7]->z1 = mid_z;
}

void update_node_moments(Tree* node, Tree* chi)
{
	#ifdef DIPOLES
	double old_cx = node->cx;
	double old_cy = node->cy;
	double old_cz = node->cz;
	#endif
	node->cx = (node->aq * node->cx + chi->aq * chi->cx) / (node->aq + chi->aq);
	node->cy = (node->aq * node->cy + chi->aq * chi->cy) / (node->aq + chi->aq);
	node->cz = (node->aq * node->cz + chi->aq * chi->cz) / (node->aq + chi->aq);
	#ifdef DIPOLES
	node->px -= (node->cx - old_cx) * node->q;
	node->py -= (node->cy - old_cy) * node->q;
	node->pz -= (node->cz - old_cz) * node->q;
	node->px += chi->px - (node->cx - chi->cx) * chi->q;
	node->py += chi->py - (node->cy - chi->cy) * chi->q;
	node->pz += chi->pz - (node->cz - chi->cz) * chi->q;
	#endif
	node->q += chi->q;
	node->aq += chi->aq;
}

void compute_node_values(Proc* pro, unsigned int node_no)
{
	Tree* chi;
	Tree* node = &(pro->root[node_no]);
	Particle* p;

	if (node->part != UINT_MAX) {
		p = &(pro->parts[node->part]);
		if (p->species == DUST) {
			node->q = pro->dust_grain_q;
		} else {
			node->q = CHARGE[p->species];
		}
		node->aq = abs(node->q);
		node->cx = p->x;
		node->cy = p->y;
		node->cz = p->z;
		#ifdef DIPOLES
		node->px = 0.0;
		node->py = 0.0;
		node->pz = 0.0;
		#endif
	} else {
		node->q = 0;
		node->aq = 0;
		node->cx = 0.0;  /* precise value here doesn't matter */
		node->cy = 0.0;  /* precise value here doesn't matter */
		node->cz = 0.0;  /* precise value here doesn't matter */
		#ifdef DIPOLES
		node->px = 0.0;  /* but this must be nil */
		node->py = 0.0;  /* but this must be nil */
		node->pz = 0.0;  /* but this must be nil */
		#endif
		if (node->gxlylz) {
			compute_node_values(pro, node->gxlylz);
			chi = &(pro->root[node->gxlylz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
		if (node->gxlygz) {
			compute_node_values(pro, node->gxlygz);
			chi = &(pro->root[node->gxlygz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
		if (node->lxlylz) {
			compute_node_values(pro, node->lxlylz);
			chi = &(pro->root[node->lxlylz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
		if (node->lxlygz) {
			compute_node_values(pro, node->lxlygz);
			chi = &(pro->root[node->lxlygz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
		if (node->lxgylz) {
			compute_node_values(pro, node->lxgylz);
			chi = &(pro->root[node->lxgylz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
		if (node->lxgygz) {
			compute_node_values(pro, node->lxgygz);
			chi = &(pro->root[node->lxgygz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
		if (node->gxgylz) {
			compute_node_values(pro, node->gxgylz);
			chi = &(pro->root[node->gxgylz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
		if (node->gxgygz) {
			compute_node_values(pro, node->gxgygz);
			chi = &(pro->root[node->gxgygz]);
			if (chi->aq) {
				update_node_moments(node, chi);
			}
		}
	}
}

/* Allocate more memory for storing a tree's node information, handling
   the necessary book-keeping for the caller, `build_tree`. */
bool allocate_more_tree_memory(Proc* pro, Tree** root, Tree** tr, unsigned int cur)
{
	Tree* new_root;

	pro->tr_alloc = (pro->tr_alloc * 1.1) + 9;
	new_root = realloc(pro->root, sizeof(Tree) * pro->tr_alloc);
	if (new_root == NULL) {
		fprintf(stderr, "%i:%lu: can't allocate %lu bytes "
				"for particle tree.\n", pro->rank, pro->t_steps,
				sizeof(Tree) * pro->tr_alloc);
		return false;
	} else {
		fprintf(stderr, "%i:%lu: allocated more particle tree memory "
				"(%u nodes' worth).\n", pro->rank, pro->t_steps,
				pro->tr_alloc);
	}

	pro->root = new_root;
	*root = pro->root;

	/* `realloc` might have moved the block of memory to which
	   `pro->root` was pointing; reset `tr` to allow for this
	   possibility. */
	*tr = &(pro->root[cur]);

	return true;
}

bool is_particle_position_nan(const Proc* pro, unsigned int part_no)
{
	const Particle* part = &(pro->parts[part_no]);

	if (isnan(part->x) || isnan(part->y) || isnan(part->z)) {
		fprintf(stderr, "%i:%lu: particle %u has position (%g, %g, %g).\n",
				pro->rank, pro->t_steps, part_no,
				part->x, part->y, part->z);
		return true;
	}

	return false;
}

/* Build a particle tree for this process' particles, particle by particle.
   (This may throw itself into an infinite loop if two particles have
    identical positions, which is physically unrealistic but can happen if
    the PRNG isn't seeded well. To save time, this function doesn't always
    test for multiple particles at the same position; the caller is trusted
    to ensure that no two particles have exactly the same positions.)
   (Oh, yeah, infinite loops have been known to occur in the function if
    a particle's position is NaN, which should never happen, but might if
    e.g. the PRNG is buggy.) */
Tree* build_tree(Proc* pro, Run ru)
{
	unsigned int cur;
	double mid_x;
	double mid_y;
	double mid_z;
	unsigned int next;
	Particle* part;
	unsigned int part_no;
	Tree* root = pro->root;
	Tree* tr;
	unsigned int used = 0;  /* number of allocated nodes used in tree */

	/* Store the root particle tree node as the 1st entry in the array
	   of particle tree nodes. */
	root->x1 = 0.0;
	root->y1 = 0.0;
	root->z1 = 0.0;
	#ifdef SPHERICAL
	root->x2 = 2.0 * ru.r * ru.mpp;
	#else
	root->x2 = ru.L;
	#endif
	root->y2 = root->x2;
	root->z2 = root->y2;
	root->part = UINT_MAX;
	root->cx = 0.0;
	root->cy = 0.0;
	root->cz = 0.0;
	root->gxlylz = 0;
	root->gxlygz = 0;
	root->lxlylz = 0;
	root->lxlygz = 0;
	root->lxgylz = 0;
	root->lxgygz = 0;
	root->gxgylz = 0;
	root->gxgygz = 0;

	#if defined(PROLATE) || defined(OBLATE)
	part_no = 1;	/* Excluding dust (part_no = 0) from tree - JH */
	#else
	part_no = 0;
	#endif

	while (part_no < ru.N) {

		part = &(pro->parts[part_no]);

		/* Check the particle's inside the top-level cell; if it's not
		   then abort the run; it'd be impossible to put it into a cell
		   and trying to do so would cause an infinite loop here. */
		if ((part->x < 0.0) || (part->x > root->x2)
		    || (part->y < 0.0) || (part->y > root->y2)
		    || (part->z < 0.0) || (part->z > root->z2)) {
			/* FIXME: if this error message is triggered, every worker
			   process displays it (because each worker independently
			   builds the particle tree itself)! */
			fprintf(stderr, "%i: error: particle %u at (%g, %g, %g) "
			        "with velocity (%g, %g, %g) not in simulation area.\n",
			        pro->rank, part_no, part->x, part->y, part->z,
			        part->vx, part->vy, part->vz);
			return NULL;
		}

		/* Traverse the tree until the leaf node that would contain
		   particle `part_no` is reached. */
		next = 0;
		do {
			cur = next;
			tr = &(pro->root[cur]);
			mid_x = (tr->x1 + tr->x2) / 2.0;
			mid_y = (tr->y1 + tr->y2) / 2.0;
			mid_z = (tr->z1 + tr->z2) / 2.0;
			if (part->x < mid_x) {
				if (part->y < mid_y) {
					if (part->z < mid_z) {
						next = tr->lxlylz;
					} else {
						next = tr->lxlygz;
					}
				} else {
					if (part->z < mid_z) {
						next = tr->lxgylz;
					} else {
						next = tr->lxgygz;
					}
				}
			} else {
				if (part->y < mid_y) {
					if (part->z < mid_z) {
						next = tr->gxlylz;
					} else {
						next = tr->gxlygz;
					}
				} else {
					if (part->z < mid_z) {
						next = tr->gxgylz;
					} else {
						next = tr->gxgygz;
					}
				}
			}
		} while (next);

		/* If this node is empty, hooray! The current particle can simply
		   be assigned to it. Otherwise, start partitioning this node. */
		if (tr->part == UINT_MAX) {
			tr->part = part_no;
			part_no++;  /* move on to the next particle */
		} else {
			/* This node already has a particle assigned to it. Partition it
			   into new nodes, move the already assigned particle to the
			   appropriate new child node, then go around the loop again to
			   try to find an empty home for particle `part_no` again. */

			if (pro->tr_alloc < (used + 9)) {
				/* Uh oh, there won't be enough memory allocated for the new
				   batch of nodes. Allocate more memory first. */
				if (!allocate_more_tree_memory(pro, &root, &tr, cur)) {
					return NULL;
				}
				/* FIXME: experimental error check to try to catch particles
				   at the same position. Take this out if it doesn't ever
				   help. What it does: checks whether this tree node's
				   particle has exactly the same position as `part`.
				   If so, fail. */
				if ((part->x == pro->parts[tr->part].x)
				    && (part->y == pro->parts[tr->part].y)
				    && (part->z == pro->parts[tr->part].z)) {
					fprintf(stderr, "%i:%lu: particles %u & %u both have "
					        "the position (%g, %g, %g).\n",
					        pro->rank, pro->t_steps, part_no, tr->part,
					        part->x, part->y, part->z);
					return NULL;
				}
				/* FIXME: more experimental error checks; see whether
				   this particle and/or the tree node particle have NaN
				   positions. */
				if (is_particle_position_nan(pro, part_no)
				    || is_particle_position_nan(pro, tr->part)) {
					return NULL;
				}
			}

			tr->gxlylz = ++used;
			tr->gxlygz = ++used;
			tr->lxlylz = ++used;
			tr->lxlygz = ++used;
			tr->lxgylz = ++used;
			tr->lxgygz = ++used;
			tr->gxgylz = ++used;
			tr->gxgygz = ++used;
			init_node_children(pro->root, cur);
			part = &(pro->parts[tr->part]);
			if (part->x < mid_x) {
				if (part->y < mid_y) {
					if (part->z < mid_z) {
						pro->root[tr->lxlylz].part = tr->part;
					} else {
						pro->root[tr->lxlygz].part = tr->part;
					}
				} else {
					if (part->z < mid_z) {
						pro->root[tr->lxgylz].part = tr->part;
					} else {
						pro->root[tr->lxgygz].part = tr->part;
					}
				}
			} else {
				if (part->y < mid_y) {
					if (part->z < mid_z) {
						pro->root[tr->gxlylz].part = tr->part;
					} else {
						pro->root[tr->gxlygz].part = tr->part;
					}
				} else {
					if (part->z < mid_z) {
						pro->root[tr->gxgylz].part = tr->part;
					} else {
						pro->root[tr->gxgygz].part = tr->part;
					}
				}
			}
			tr->part = UINT_MAX;
		}

	}

	compute_node_values(pro, 0);

	return pro->root;
}

/* Display an error message and immediately halt; traversing the tree calls
   for too long a path through the tree. */
void particle_tree_path_too_long(int rank, const char* reason)
{
	fprintf(stderr, "%i: particle tree path length exceeded %u "
	        "when estimating %s.\n", rank, MAX_PATH_LEN, reason);
	abort();
}

/* Funny preprocessor if/or/or statement. The function which determines field
   from grain must be called if particle-particle interactions are off or
   if the dust produces a special, e.g. spheroidal, potential - JH */
#if defined(DUST_FIELD_ONLY) || defined(PROLATE) || defined(OBLATE)

/* Compute the electric field a particle feels from the central dust grain.
   (This function modifies its `E` argument, and needs the particle's array
    index `part` to check whether it's zero, i.e. whether the calculation
    is of the grain's self-interaction, which should be nil!) */
void field_particle_feels_from_grain(unsigned int part, const Proc* pro,
		const Run* ru, double x, double y, double z, long double* E)
{
	#if defined(PROLATE) || defined (OBLATE)
	double rho2;
	double rt;
	double E_n;
	double beta;
	double theta;
	double dx;
	double dy;
	double dz;
	#ifdef PROLATE
		double sinh2;
		double c = ( ru->a * sqrt(ru->A*ru->A - 1.0) ) *
			pow(ru->A, -1.0/3.0);
		double phi_const = (QE * pro->dust_grain_q * acosh(ru->A) ) / 
			(4.0 * PI * EPSILON0 * c * 
			log(ru->A - sqrt(ru->A * ru->A -1.0)) );
	#endif
	#ifdef OBLATE
		double cosh2;
		double c = ( ru->a * sqrt(1.0 - ru->A*ru->A) ) *
			pow(ru->A, -1.0/3.0);
		double phi_const = (QE * pro->dust_grain_q * acos(ru->A) ) / 
			(8.0 * PI * EPSILON0 * c * 
			(atan( (1.0 - sqrt(1.0 - ru->A * ru->A)) / ru->A )
				- 0.25*PI ) );
	#endif
	const Particle* p = &(pro->parts[0]);

	if (!part) {
		/* Particle `part` is the central dust grain! */
		E[0] = 0.0;
		E[1] = 0.0;
		E[2] = 0.0;
		return;
	}

	dx = x - p->x;
	dy = y - p->y;
	dz = z - p->z;
	rho2 = (dy * dy) + (dz * dz);
	#ifdef PROLATE
		rt = sqrt( rho2*rho2 + dx*dx*dx*dx + c*c*c*c + 2.0*rho2*dx*dx
			+ 2.0*rho2*c*c - 2.0*dx*dx*c*c );
		sinh2 = (rho2 + dx*dx - c*c + rt) / (2.0*c*c);
		E_n = -phi_const / sqrt( sinh2 * rt );
		beta = atan( (-sinh2*dx) / ( (sinh2+1.0) * sqrt(rho2) ) );
	#endif
	#ifdef OBLATE
		rt = sqrt( rho2*rho2 + dx*dx*dx*dx + c*c*c*c + 2.0*rho2*dx*dx
			- 2.0*rho2*c*c + 2.0*dx*dx*c*c );
		cosh2 = (rho2 + dx*dx + c*c + rt) / (2.0*c*c);
		E_n = -phi_const / sqrt( cosh2 * rt );
		beta = atan( (-cosh2*dx) / ( (cosh2-1.0) * sqrt(rho2) ) );
	#endif
	theta = atan( fabs(dz)/fabs(dy) );

	E[0] = -E_n * sin(beta);
	E[1] = E_n * cos(beta) * cos(theta) * fabs(dy) / dy;
	E[2] = E_n * cos(beta) * sin(theta) * fabs(dz) / dz;

	#else

	double dist2;
	double dx;
	double dy;
	double dz;
	long double E_scale = QE * pro->dust_grain_q;
	const Particle* p = &(pro->parts[0]);

	if (!part) {
		/* Particle `part` is the central dust grain! */
		E[0] = 0.0;
		E[1] = 0.0;
		E[2] = 0.0;
		return;
	}

	dx = x - p->x;
	dy = y - p->y;
	dz = z - p->z;
	dist2 = (dx * dx) + (dy * dy) + (dz * dz);
	dist2 += ru->soft * ru->soft;  /* FIXME: is this necessary here? */
	E_scale /= ((4 * PI * EPSILON0) * dist2 * sqrt(dist2));
	E[0] = E_scale * dx;
	E[1] = E_scale * dy;
	E[2] = E_scale * dz;
	#endif
}

#endif

#if !defined(DUST_FIELD_ONLY) || defined(PROLATE) || defined(OBLATE)

/* Compute the electric field a particle feels from the particles in a
   particle (sub)tree. (This function modifies its `E` argument.)
   (This function needs the array index of the particle so it can exclude
    that particle's own field from the calculation, and it needs an explicit
    position because the numerical solver might need the fields the particle
    would experience at a hypothetical, different position.)

   If `PROLATE` or `OBLATE` is defined then the dust field is evaluated 
   outside the tree with a spheroidal potential - JH */
void field_particle_feels_at(unsigned int part, const Proc* pro, const Run* ru,
		double x, double y, double z, long double* E, Tree* root)
{
	double dist;
	double dist2;
	double dx;
	double dy;
	double dz;
	long double E_scale;
	Tree* node;
	#ifdef DIPOLES
	double r_dot_p_d2;
	#endif
	Particle* p;
	unsigned int path[MAX_PATH_LEN];
	unsigned int path_len = 1;

	E[0] = 0.0;
	E[1] = 0.0;
	E[2] = 0.0;

	#if defined(PROLATE) || defined(OBLATE)
	/* If spheroid, initialise E with grain's field - JH */
	field_particle_feels_from_grain(part, pro, ru, x, y, z, E);
	#endif

	path[0] = 0;

	do {
		node = &(root[path[path_len-1]]);
		path_len--;

		if (node->part == UINT_MAX) {
			/* This cell contains multiple particles (or none at all!),
			   so try using the simpler particle-cell interaction. */
			dx = x - node->cx;
			dy = y - node->cy;
			dz = z - node->cz;
			dist2 = (dx * dx) + (dy * dy) + (dz * dz);
			if (((node->x2 - node->x1) * (node->x2 - node->x1)) < (OAP * OAP * dist2)) {
				/* Opening angle criterion satisfied: calculate field from
				   cell without breaking it down further. */
				#ifdef DIPOLES
				if (node->aq) {
				#else
				if (node->q) {
				#endif
					dist = sqrt(dist2);
					E_scale = QE / (4 * PI * EPSILON0 * dist2 * dist);
					#ifdef DIPOLES
					r_dot_p_d2 = (dx * node->px) + (dy * node->py) + (dz * node->pz);
					r_dot_p_d2 /= dist2;
					E[0] += E_scale * ((node->q * dx) - node->px + (3 * dx * r_dot_p_d2));
					E[1] += E_scale * ((node->q * dy) - node->py + (3 * dy * r_dot_p_d2));
					E[2] += E_scale * ((node->q * dz) - node->pz + (3 * dz * r_dot_p_d2));
					#else
					E_scale *= node->q;
					E[0] += E_scale * dx;
					E[1] += E_scale * dy;
					E[2] += E_scale * dz;
					#endif
				}
			} else {
				/* Opening angle criterion not satisfied: add this cell's
				   subcells to the list of cells to visit. */
				if ((path_len + 8) > MAX_PATH_LEN) {
					particle_tree_path_too_long(pro->rank, "fields");
				}
				if (node->gxlylz) {
					/* This cell points to a subcell in its `gxlylz` member.
					   Assume this cell's in a consistent state, i.e. that
					   it's safe to assume that this cell's other 7 subcell
					   array index members are correctly set. This saves 7
					   more if statement checks. */
					path[path_len++] = node->gxlylz;
					path[path_len++] = node->gxlygz;
					path[path_len++] = node->lxlylz;
					path[path_len++] = node->lxlygz;
					path[path_len++] = node->lxgylz;
					path[path_len++] = node->lxgygz;
					path[path_len++] = node->gxgylz;
					path[path_len++] = node->gxgygz;
				}
			}
		} else if (node->part != part) {
			/* This cell contains one particle and hence no subcells, so
			   directly compute the field from that particle. */
			p = &(pro->parts[node->part]);
			dx = x - p->x;
			dy = y - p->y;
			dz = z - p->z;
			dist2 = (dx * dx) + (dy * dy) + (dz * dz);
			if (dist2 != 0.0) {
				/* Note that the other particle only influences this particle
				   with a field if they're at different positions (which
				   would ideally always be the case). If the particles have
				   literally exactly the same position (which can happen with
				   an unlucky PRNG) they don't notice each other...like ships
				   in the night. */
				dist2 += ru->soft * ru->soft;
				dist = sqrt(dist2);
				if (p->species == DUST) {
					E_scale = QE * pro->dust_grain_q / ((4 * PI * EPSILON0) * dist2 * dist);
				} else {
					E_scale = QE * CHARGE[p->species] / ((4 * PI * EPSILON0) * dist2 * dist);
				}
				E[0] += E_scale * dx;
				E[1] += E_scale * dy;
				E[2] += E_scale * dz;
			} else {
				fprintf(stderr, "%i: multiple particles at the location"
				        " (%g, %g, %g)!\n", pro->rank, x, y, z);
				abort();
			}
		}

	} while (path_len);

}

#endif

#ifdef USING_PHI_AT
/* Compute the net electrostatic potential at a point from the particles
   in a particle (sub)tree.
   FIXME: it feels gross to repeat so much code from `phi_particle_feels`
   but I haven't immediately thought of a nice way to avoid the repetition
   without sacrificing performance. */
double phi_at(double x, double y, double z, const Proc* pro, const Run* ru, Tree* root)
{
	#if defined(PROLATE) || defined(OBLATE)
	#ifdef PROLATE
		double c = ( ru->a * sqrt(ru->A*ru->A - 1.0) ) / pow(ru->A, 1.0/3.0);
		double phi_const = (QE * pro->dust_grain_q * acosh(ru->A) ) / 
			(4.0*PI*EPSILON0*c*log(ru->A - sqrt(ru->A * ru->A -1.0)) );
		double rho2 = (y * y) + (z * z);
		double rt = sqrt( rho2*rho2 + x*x*x*x + c*c*c*c + 2.0*rho2*x*x
			+ 2.0*rho2*c*c - 2.0*x*x*c*c );
		double phi = phi_const * log(
			( sqrt(rho2 + x*x + c*c + rt) - c*sqrt(2.0) ) / 
			sqrt(rho2 + x*x - c*c + rt) );
	#endif
	#ifdef OBLATE
		double c = ( ru->a * sqrt(1.0 - ru->A*ru->A) ) *
			pow(ru->A, -1.0/3.0);
		double phi_const = (QE * pro->dust_grain_q * acos(ru->A) ) / 
			(4.0 * PI * EPSILON0 * c * 
			(atan( (1.0 - sqrt(1.0 - ru->A * ru->A)) / ru->A )
				- 0.25*PI ) );
		double rho2 = (y * y) + (z * z);
		double rt = sqrt( rho2*rho2 + x*x*x*x + c*c*c*c + 2.0*rho2*x*x
			- 2.0*rho2*c*c + 2.0*x*x*c*c );
		double phi = phi_const * ( atan(
			( sqrt(rho2 + x*x + c*c + rt) - c*sqrt(2.0) ) / 
			sqrt(rho2 + x*x - c*c + rt) ) - 0.25*PI );
	#endif
	#else
	double phi = 0.0;
	#endif

	double dist;
	double dist2;
	double dx;
	double dy;
	double dz;
	Tree* node;
	Particle* p;
	unsigned int path[MAX_PATH_LEN];
	unsigned int path_len = 1;
	#ifdef DIPOLES
	double r_dot_p;
	#endif

	path[0] = 0;

	do {
		node = &(root[path[path_len-1]]);

		if (node->part != UINT_MAX) {
			/* This cell contains one particle and hence no subcells, so
			   directly compute that particle's potential. */
			p = &(pro->parts[node->part]);
			dist2 = ((x - p->x) * (x - p->x)) + ((y - p->y) * (y - p->y)) + ((z - p->z) * (z - p->z));
			dist = sqrt(dist2 + (ru->soft * ru->soft));
			if (p->species == DUST) {
				phi += QE * pro->dust_grain_q / (4 * PI * EPSILON0 * dist);
			} else {
				phi += QE * CHARGE[p->species] / (4 * PI * EPSILON0 * dist);
			}
			path_len--;
		} else {
			dx = x - node->cx;
			dy = y - node->cy;
			dz = z - node->cz;
			dist2 = (dx * dx) + (dy * dy) + (dz * dz);
			if (((node->x2 - node->x1) * (node->x2 - node->x1)) < (OAP * OAP * dist2)) {
				/* Opening angle criterion satisfied: calculate this
				   cell's potential without breaking it down further. */
				#ifdef DIPOLES
				if (node->aq) {
				#else
				if (node->q) {
				#endif
					#ifdef DIPOLES
					dist = sqrt(dist2);
					r_dot_p = (dx * node->px) + (dy * node->py) + (dz * node->pz);
					phi += (QE / (4 * PI * EPSILON0)) * ((node->q / dist) + (r_dot_p / (dist * dist2)));
					#else
					phi += QE * node->q / (4 * PI * EPSILON0 * sqrt(dist2));
					#endif
				}
				path_len--;
			} else {
				/* Opening angle criterion not satisfied: add this cell's
				   subcells to the list of cells to visit. */
				path_len--;
				if ((path_len + 8) > MAX_PATH_LEN) {
					particle_tree_path_too_long(pro->rank, "phi");
				}
				if (node->gxlylz) {
					/* This cell points to a subcell in its `gxlylz` member.
					   Assume this cell's in a consistent state, i.e. that
					   it's safe to assume that this cell's other 7 subcell
					   array index members are correctly set. This saves 7
					   more if statement checks. */
					path[path_len++] = node->gxlylz;
					path[path_len++] = node->gxlygz;
					path[path_len++] = node->lxlylz;
					path[path_len++] = node->lxlygz;
					path[path_len++] = node->lxgylz;
					path[path_len++] = node->lxgygz;
					path[path_len++] = node->gxgylz;
					path[path_len++] = node->gxgygz;
				}
			}
		}

	} while (path_len);

	return phi;
}
#endif

#if defined(DUST_FIELD_ONLY) || defined(PROLATE) || defined(OBLATE)
/* Compute the electric potential a particle experiences from the central
   dust grain. */
long double phi_particle_feels_from_grain(unsigned int part, const Proc* pro, const Run* ru)
{
	double dx;
	double dy;
	double dz;
	Particle* p = &(pro->parts[0]);

	if (!part) {
		return 0.0L;
	}

	dx = pro->parts[part].x - p->x;
	dy = pro->parts[part].y - p->y;
	dz = pro->parts[part].z - p->z;

	#ifdef PROLATE
		double c = ru->a * sqrt(ru->A*ru->A - 1.0) / pow(ru->A, 1.0/3.0);
		double phi_const = (QE * pro->dust_grain_q * acosh(ru->A) ) / 
			(4.0*PI*EPSILON0*c*log(ru->A - sqrt(ru->A * ru->A -1.0)) );
		double rho2 = (dy * dy) + (dz * dz);
		double rt = sqrt( rho2*rho2 + dx*dx*dx*dx + c*c*c*c + 2.0*rho2*dx*dx
			+ 2.0*rho2*c*c - 2.0*dx*dx*c*c );
		return phi_const * log(
			( sqrt(rho2 + dx*dx + c*c + rt) - c*sqrt(2.0) ) / 
			sqrt(rho2 + dx*dx - c*c + rt) );
	#elif defined(OBLATE)
		double c = ( ru->a * sqrt(1.0 - ru->A*ru->A) ) *
			pow(ru->A, -1.0/3.0);
		double phi_const = (QE * pro->dust_grain_q * acos(ru->A) ) / 
			(4.0 * PI * EPSILON0 * c * 
			(atan( (1.0 - sqrt(1.0 - ru->A * ru->A)) / ru->A )
				- 0.25*PI ) );
		double rho2 = (dy * dy) + (dz * dz);
		double rt = sqrt( rho2*rho2 + dx*dx*dx*dx + c*c*c*c + 2.0*rho2*dx*dx
			- 2.0*rho2*c*c + 2.0*dx*dx*c*c );
		return phi_const * ( atan(
			( sqrt(rho2 + dx*dx + c*c + rt) - c*sqrt(2.0) ) / 
			sqrt(rho2 + dx*dx - c*c + rt) ) - 0.25*PI );
	#else
	double dist2 = (dx * dx) + (dy * dy) + (dz * dz);
	dist2 += (ru->soft * ru->soft);
	return QE * pro->dust_grain_q / (4 * PI * EPSILON0 * sqrt(dist2));
	#endif
}

#endif
#if !defined(DUST_FIELD_ONLY) || defined(PROLATE)

/* Compute the electric potential a particle experiences from the
   particles in a particle (sub)tree.
   (This function needs the particle's array index to exclude that
   particle's own potential from the calculation.) */
long double phi_particle_feels(unsigned int part, const Proc* pro, const Run* ru, Tree* root)
{
	double dist;
	double dist2;
	double dx;
	double dy;
	double dz;
	Tree* node;
	Particle* p;
	unsigned int path[MAX_PATH_LEN];
	unsigned int path_len = 1;
	long double phi = 0.0;
	#ifdef DIPOLES
	double r_dot_p;
	#endif
	#if defined(PROLATE) || defined(OBLATE)
	phi = phi_particle_feels_from_grain(part, pro, ru);
	#endif

	path[0] = 0;

	do {
		node = &(root[path[path_len-1]]);
		p = &(pro->parts[node->part]);
		path_len--;

		if (node->part == UINT_MAX) {
			dx = pro->parts[part].x - node->cx;
			dy = pro->parts[part].y - node->cy;
			dz = pro->parts[part].z - node->cz;
			dist2 = (dx * dx) + (dy * dy) + (dz * dz);
			if (((node->x2 - node->x1) * (node->x2 - node->x1)) < (OAP * OAP * dist2)) {
				/* Opening angle criterion satisfied: calculate this
				   cell's potential without breaking it down further. */
				#ifdef DIPOLES
				if (node->aq) {
				#else
				if (node->q) {
				#endif
					#ifdef DIPOLES
					dist = sqrt(dist2);
					r_dot_p = (dx * node->px) + (dy * node->py) + (dz * node->pz);
					phi += QE * ((node->q / dist) + (r_dot_p / (dist * dist2))) / (4 * PI * EPSILON0);
					#else
					phi += QE * node->q / (4 * PI * EPSILON0 * sqrt(dist2));
					#endif
				}
			} else {
				/* Opening angle criterion not satisfied: add this cell's
				   subcells to the list of cells to visit. */
				if ((path_len + 8) > MAX_PATH_LEN) {
					particle_tree_path_too_long(pro->rank, "phi");
				}
				if (node->gxlylz) {
					/* This cell points to a subcell in its `gxlylz` member.
					   Assume this cell's in a consistent state, i.e. that
					   it's safe to assume that this cell's other 7 subcell
					   array index members are correctly set. This saves 7
					   more if statement checks. */
					path[path_len++] = node->gxlylz;
					path[path_len++] = node->gxlygz;
					path[path_len++] = node->lxlylz;
					path[path_len++] = node->lxlygz;
					path[path_len++] = node->lxgylz;
					path[path_len++] = node->lxgygz;
					path[path_len++] = node->gxgylz;
					path[path_len++] = node->gxgygz;
				}
			}
		} else if (node->part != part) {
			/* This cell contains one particle and hence no subcells, so
			   directly compute that particle's potential. */
			dx = pro->parts[part].x - p->x;
			dy = pro->parts[part].y - p->y;
			dz = pro->parts[part].z - p->z;
			dist2 = (dx * dx) + (dy * dy) + (dz * dz);
			dist = sqrt(dist2 + (ru->soft * ru->soft));
			if (p->species == DUST) {
				phi += QE * pro->dust_grain_q / (4 * PI * EPSILON0 * dist);
			} else {
				phi += QE * CHARGE[p->species] / (4 * PI * EPSILON0 * dist);
			}
		}

	} while (path_len);

	return phi;
}

#endif

#ifdef OUT_BOUNCE
int check_particle_in_sim_area(Particle* part, Run ru, double prev_x, double prev_y, double prev_z)
#else
int check_particle_in_sim_area(Particle* part, Run ru, Proc* pro)
#endif
{
	bool reset_trajectory = false;

	#ifdef SPHERICAL

		/* The simulation domain's radius in absolute units
		   instead of pixels. */
		double rad = ru.r * ru.mpp;

		#ifdef OUT_BOUNCE

		double r[3];  /* particle's position, relative to the simulation
						 area's centre, in spherical coordinates */
		double v[3];  /* particle's velocity in spherical coordinates */

		/* Check whether the particle is outside the simulation area;
		   if so, try to bounce it off the wall back into the area. If
		   the particle has an inconsistent position-velocity combination
		   (i.e. it was embedded in the wall but travelling outwards)
		   return an error code representing this. */
		r[0] = pythag(part->x - rad, part->y - rad, part->z - rad);
		if (r[0] > rad) {
			r[1] = atan2l(part->y - rad, part->x - rad);
			r[2] = acosl((part->z - rad) / r[0]);
			v[0] = part->vx * cos(r[1]) * sin(r[2])
				   + part->vy * sin(r[1]) * sin(r[2])
				   + part->vz * cos(r[2]);
			v[1] = -part->vx * sin(r[1]) + part->vy * cos(r[1]);
			v[2] = part->vx * cos(r[1]) * cos(r[2])
				   + part->vy * sin(r[1]) * cos(r[2]) - part->vz * sin(r[2]);
			if (v[0] > 0.0) {
				v[0] *= -OUT_BOUNCE;
				v[1] *= OUT_BOUNCE;
				v[2] *= OUT_BOUNCE;
				velocity_in_cartesian_coords(v, r[1], r[2], &(part->vx),
				                             &(part->vy), &(part->vz));
			} else {
				return -1;
			}
			part->x = prev_x;
			part->y = prev_y;
			part->z = prev_z;
		}

		#else
		/* If the particle went outside the simulation area, reinject it
		   at a random point on the outer boundary. */
		if (pythag(part->x - rad, part->y - rad, part->z - rad) > rad) {
			inject(part, ru, true, pro);
			reset_trajectory = true;
		}
		#endif

	#else

		#ifdef OUT_BOUNCE

		/* If the particle's outside the simulation area, bounce it off the
		   the wall back into the area. If the particle has a seemingly
		   inconsistent position-velocity combination (e.g. it was in
		   the left wall travelling rightwards!) return an error code
		   representing this. Note that the particle could have such
		   inconsistent motion twice over (e.g. if it was in the top & left
		   walls travelling downward & right) but then only the first such
		   inconsistency is detected by this function before it returns. */
		if (part->x < 0.0) {
			part->x = prev_x;
			if (part->vx < 0.0) {
				part->vx *= -OUT_BOUNCE;
			} else {
				return -1;
			}
		} else if (part->x > (ru.L * ru.mpp)) {
			part->x = prev_x;
			if (part->vx > 0.0) {
				part->vx *= -OUT_BOUNCE;
			} else {
				return -2;
			}
		}
		if (part->y < 0.0) {
			part->y = prev_y;
			if (part->vy < 0.0) {
				part->vy *= -OUT_BOUNCE;
			} else {
				return -3;
			}
		}
		if (part->y > (ru.L * ru.mpp)) {
			part->y = prev_y;
			if (part->vy > 0.0) {
				part->vy *= -OUT_BOUNCE;
			} else {
				return -4;
			}
		}
		if (part->z < 0.0) {
			part->z = prev_z;
			if (part->vz < 0.0) {
				part->vz *= -OUT_BOUNCE;
			} else {
				return -5;
			}
		}
		if (part->z > (ru.L * ru.mpp)) {
			part->z = prev_z;
			if (part->vz > 0.0) {
				part->vz *= -OUT_BOUNCE;
			} else {
				return -6;
			}
		}

		#else

		/* If the particle went outside the simulation area, just reinject
		   it at the opposite wall with a new random velocity. */
		if ((part->x < 0.0) || (part->x > (ru.L * ru.mpp))
		    || (part->y < 0.0) || (part->y > (ru.L * ru.mpp))
		    || (part->z < 0.0) || (part->z > (ru.L * ru.mpp))) {
			inject(part, ru, true, pro);
			reset_trajectory = true;
	    }

		#endif

	#endif

	return (int) reset_trajectory;
}

#ifdef GL
/* Reset the record of a tracked particle's trajectory.
   Doing this when a particle goes off-screen or gets absorbed prevents
   the display of weird-looking tracks that abruptly jump across the
   simulation area. */
void reset_trajectory(Proc* pro, const Run* ru, unsigned int traj)
{
	if (traj >= ru->track) {
		/* This trajectory array index is too big to be correct. */
		return;
	}
	memset(pro->traj_x[traj], 0, sizeof(float) * TRAJ_LEN);
	memset(pro->traj_y[traj], 0, sizeof(float) * TRAJ_LEN);
	memset(pro->traj_z[traj], 0, sizeof(float) * TRAJ_LEN);
}
#endif

#ifdef GL
/* Return the index in the trajectory array where the particle with index
   `part_idx` in the particle array has its trajectory stored. */
unsigned int part_idx_to_traj_idx(const Run* ru, unsigned int part_idx)
{
	unsigned int traj_idx;

	if (ru->track) {
		if (!((part_idx - 1) % ((ru->N - 1) / ru->track))) {
			traj_idx = (part_idx - 1) / ((ru->N - 1) / ru->track);
			if (traj_idx >= ru->track) {
				return UINT_MAX;
			} else {
				return traj_idx;
			}
		}
	}

	return UINT_MAX;
}
#endif

#ifdef DEBUG_FIELDS
void assess_field_accuracy(unsigned int i, const Particle* p, long double* E, Particle* parts, const Run* ru, unsigned long int time_steps)
{
	long double dist;
	long double dist2;
	long double dx;
	long double dy;
	long double dz;
	long double E_pps[3];
	long double E_scale;
	unsigned int j;

	E_pps[0] = 0.0;
	E_pps[1] = 0.0;
	E_pps[2] = 0.0;
	#ifdef DISABLE_GRAIN
	for (j = 1; j < ru->N; j++) {
	#else
	for (j = 0; j < ru->N; j++) {
	#endif
		if (j == i) {
			continue;
		}
		dx = p->x - parts[j].x;
		dy = p->y - parts[j].y;
		dz = p->z - parts[j].z;
		dist2 = (dx * dx) + (dy * dy) + (dz * dz) + (ru->soft * ru->soft);
/*
		dist2 = (dx * dx) + (dy * dy) + (dz * dz);
*/
		dist = sqrt(dist2);
		/* FIXME: not right to use `CHARGE` for the dust grain, really. */
		E_scale = QE * CHARGE[parts[j].species] / (4 * PI * EPSILON0 * dist2 * dist);
		E_pps[0] += E_scale * dx;
		E_pps[1] += E_scale * dy;
		E_pps[2] += E_scale * dz;
	}
	printf("%Lg %u %i %Lg %Lg %Lg %Lg %Lg %Lg\n",
	       ru->dt * time_steps, i, p->species,
	       E[0], E_pps[0], E[1], E_pps[1], E[2], E_pps[2]);
}
#endif

#ifdef SCEPTIC_REINJECTION
/* Sample the normalized electric potential chi_b at a particular location
   on the simulation boundary, decided by this process' MPI rank.
   (Since the SCEPTIC-like reinjection algorithm assumes spherical symmetry
   of the electric potential anyway, there seems to be little point in
   sampling the electric potential from many angles, so simply sample in a
   circle in the x-y plane at equally-spaced azimuthal angles, one point for
   each process.) */
#ifdef DUST_FIELD_ONLY
double measure_chi_b(const Proc* pro, Run ru)
#else
double measure_chi_b(const Proc* pro, Run ru, Tree* root)
#endif
{
	#ifdef DUST_FIELD_ONLY
	double r2;
	#else
	double azi;  /* azimuthal angle */
	const double x_grain = pro->parts[0].x;
	const double y_grain = pro->parts[0].y;
	const double z_grain = pro->parts[0].z;
	const double zen = PI / 2.0;  /* zenith angle */
	#endif

	double phi;  /* unnormalized electric potential */
	const double r = ru.r * ru.mpp;  /* simulation radius in SI units */

	#ifdef DUST_FIELD_ONLY
	r2 = (r * r) + (ru.soft * ru.soft);
	phi = QE * pro->dust_grain_q / (4 * PI * EPSILON0 * sqrt(r2));
	#else
	azi = 2.0 * PI * pro->rank / pro->num_workers;
	phi = phi_at(x_grain + (r * cos(azi) * sin(zen)),
		         y_grain + (r * sin(azi) * sin(zen)),
		         z_grain + (r * cos(zen)),
		         pro, &ru, root);
	#endif

	return phi / (KB * ru.Ti / (CHARGE[ION] * QE));
}
#endif

/* Does a particle that's just been stepped cross the grain's surface?
   Two checks answer this.
   (1) Does the particle's new position put it inside the grain?
   (2) Did the particle cross the grain somewhere along its interpolated
       trajectory between its new position and its previous position?  */
bool does_particle_cross_grain(const Particle* p, const double* d_pos, const Particle* gr, double a
	#if defined(PROLATE) || defined(OBLATE)
	, double A
	#endif
	)
{
	#if defined(PROLATE) || defined(OBLATE)
	/* Working in x-rho plane, spheroid surface is given by
	   f = x^2 +rho^2*A^2 = a^2*A^(4/3). After check (1), use y = m*x + c
	   and z = M*x + C to write alpha in terms of x only. Then df/dx = 0
	   yields the `closest approach` of the particle for check (2) - JH */

	double x_min;
	double x_prev;
	double y_min;
	double z_min;
	double s[3];  /* grain's pos'n vector minus particle's pos'n vec. */

	s[0] = gr->x - p->x;
	s[1] = gr->y - p->y;
	s[2] = gr->z - p->z;

	/* Check (1) */ 
	if ( ( s[0]*s[0] + (s[1]*s[1] + s[2]*s[2])*A*A ) <= 
				(a * a * pow(A, 4.0/3.0) ) ) {
		return true;
	}

	/* Check (2) */
	x_min = ( d_pos[1]*(d_pos[1]*s[0] - d_pos[0]*s[1]) +
			d_pos[2]*(d_pos[2]*s[0] - d_pos[0]*s[2]) ) /
			( (d_pos[0]*d_pos[0]/(A*A)) + d_pos[1]*d_pos[1] +
			d_pos[2]*d_pos[2] );
	x_prev = s[0] - d_pos[0];

	if ( ( (x_min > s[0]) && (x_min > x_prev) ) || 
		( (x_min < s[0]) && (x_min < x_prev) ) ) {
		/* The point of closest approach isn't between the particle's
		   current & previous positions, so it doesn't count. */
		return false;
	}

	/* The point of closest approach is between the particle's previous &
	   current positions, but is it close enough to breach the grain? */
	y_min = d_pos[1]*(x_min - s[0])/d_pos[0] + s[1];
	z_min = d_pos[2]*(x_min - s[0])/d_pos[0] + s[2];

	return ( (x_min*x_min) + (y_min*y_min + z_min*z_min) * A * A )
				<= (a * a * pow( A, 4.0/3.0) );

	#else
	/* Drew works his magic with some vector algebra - JH */
	double alpha;
	double d_dot_s;  /* dot product of `d_pos` and `s` */
	double d_pos2;   /* squared length of `d_pos` */
	double dist2;
	double s[3];  /* grain's pos'n vector minus particle's pos'n vec. */
	double s2;    /* squared length of `s` */

	s[0] = gr->x - p->x;
	s[1] = gr->y - p->y;
	s[2] = gr->z - p->z;
	s2 = (s[0] * s[0]) + (s[1] * s[1]) + (s[2] * s[2]);

	/* Check 1: does the particle's new position put it inside the grain? */
	if (s2 <= (a * a)) {
		return true;
	}

	/* Begin check 2: did the particle breach the grain during the step? */
	d_pos2 = (d_pos[0] * d_pos[0]) + (d_pos[1] * d_pos[1])
	         + (d_pos[2] * d_pos[2]);
	d_dot_s = (d_pos[0] * s[0]) + (d_pos[1] * s[1]) + (d_pos[2] * s[2]);

	/* To solve this problem, extend the line connecting the particle's
	   current position and old position (the vector `d_pos`) to infinity
	   and see where that line is closest to the grain. */
	alpha = d_dot_s / d_pos2;
	if ((alpha < -1.0) || (alpha > 0.0)) {
		/* The point of closest approach isn't between the particle's
		   current & previous positions, so it doesn't count. */
		return false;
	}

	/* The point of closest approach is between the particle's previous &
	   current positions, but is it close enough to breach the grain? */
	dist2 = s2 - (alpha * d_dot_s);
	return dist2 <= (a * a);

	#endif
}

#ifdef DUST_FIELD_ONLY
void move_particles(Proc* pro, Run ru, MPI_Datatype MPIParticle, bool compute_epe, unsigned long int time_steps)
#else
void move_particles(Proc* pro, Run ru, Tree* root, MPI_Datatype MPIParticle, bool compute_epe, unsigned long int time_steps)
#endif
{
	const double B[3] = { EXTERNAL_B_X, EXTERNAL_B_Y, EXTERNAL_B_Z };
	long double dt;
	long double E[3];
	unsigned int i;
	unsigned int ims;  /* `i` minus `start`, used as array index in loops */
	Particle* p;
	Particle* p_n;
	long double qom;   /* "Q over m", i.e. charge divided by mass */
	unsigned int* p_idxs;
	MPI_Status recv_status;
	unsigned int start;

	#ifndef DISABLE_GRAIN
	unsigned int absorbed_e = 0;   /* electrons absorbed by grain */
	long double absorbed_m = 0.0;  /* mass absorbed by grain */
	int absorbed_q = 0;     /* elementary charges absorbed by grain */
	double delta_pos[3];    /* used to interpolate p'cle trajectory */
	long double dpx = 0.0;  /* grain x-momentum gained by absorption */
	long double dpy = 0.0;  /* grain y-momentum gained by absorption */
	long double dpz = 0.0;  /* grain z-momentum gained by absorption */
	unsigned int net_abs_e = 0;   /* electrons absorbed (proc. 0 only) */
	long double net_abs_m = 0.0;  /* net grain mass gained (pr. 0 only) */
	int net_abs_q = 0.0;          /* net charge grain abs'd (pr. 0 only) */
	long double net_dpx = 0.0;  /* net x-momentum grain gained (pr. 0 only) */
	long double net_dpy = 0.0;  /* net y-momentum grain gained (pr. 0 only) */
	long double net_dpz = 0.0;  /* net z-momentum grain gained (pr. 0 only) */
	#endif

	#ifdef DEBUG_RUTHERFORD
	double old_v[3];  /* electron's velocity at start of time step */
	#endif

	#ifdef POIS_POS_COLLECTION
	double phi_normed;  /* grain surface potential norm'd by T_i */
	double kb_ti;       /* cached ion temperature as an energy */
	double lambda;      /* rate of Poissonian ion collection */
	#endif

	#ifdef BORIS
	/* Precompute `B`'s magnitude, `B` as a direction vector, and the
	   delta phi rotation angle for the `B`-based rotation sub-step. */
	const double B_mag = pythag(B[0], B[1], B[2]);
	const double b[3] = { B[0] / B_mag, B[1] / B_mag, B[2] / B_mag };
	double dephi;
	#else
	/* Acceleration vector for the second half-step, used by the non-Boris
	   particle motion integrators. */
	long double a_x2;
	long double a_y2;
	long double a_z2;
	#endif

	/* None of the surplus processors should be executing this function,
	   but just in case they are, make sure they stop at this point. */
	if (((unsigned int) pro->rank) >= ru.N) {
		fprintf(stderr, "%i: called move_particles function unnecessarily.\n",
		        pro->rank);
		return;
	}

	/* If execution reaches this point, this process isn't idle -- it's
	   doing work. So it's now confirmed to be safe to dereference
	   `pro->pits` and set the convenience variable `p_idxs`, which is the
	   starting & ending particle array index for the particles this
	   process in is responsible for. */
	p_idxs = pro->pits[pro->rank];

	/* Copy the current particle array into the to-be-updated particle
	   array. */
	memcpy(pro->new_parts, pro->parts, sizeof(Particle) * ru.N);

	#ifdef SCEPTIC_REINJECTION
	/* SCEPTIC-like reinjection requires an estimate of the normalized
	   electric potential at the simulation boundary `chi_b`. Periodically
	   make such an estimate by averaging together samples taken by
	   each process. */
	if (!(pro->t_steps % 20)) {
		#ifdef DUST_FIELD_ONLY
		pro->chi_b = measure_chi_b(pro, ru);
		#else
		pro->chi_b = measure_chi_b(pro, ru, root);
		#endif
		MPI_Allreduce(MPI_IN_PLACE, &(pro->chi_b), 1,
		              MPI_DOUBLE, MPI_SUM, pro->work_comm);
		pro->chi_b /= pro->num_workers;
	}
	#endif

	/* Some particles probably don't need to be moved on this time step;
	   deduce which ones these are and move the array index to start
	   stepping particles from the right place in this process' array.
	   Note the implicit assumption that particles are arranged heaviest
	   first in the array.
	   FIXME?: as is, this implicitly assumes that if one species is due
	   to be stepped, all lighter species need to be stepped too, which
	   isn't true. Fortunately, this will merely improve the quality of
	   the numeric integration slightly while eating a little more CPU. */
	for (start = p_idxs[0]; start <= p_idxs[1]; start++) {
		if (!(time_steps % pro->tss[pro->parts[start].species])) {
			break;
		}
	}

	/* The worker processes tell their `start` values to the master, so
	   the master can deduce how many particles each process will actually
	   be moving (and hence how much updated particle data it's going
	   to receive later in this function). */
	MPI_Gather(&start, 1, MPI_UNSIGNED, pro->starts, 1, MPI_UNSIGNED,
	           0, pro->work_comm);

	/* The particles' initial positions may need to be referred to later to
	   discern whether they crossed the grain's surface during a time step.
	   So copy the particles' initial positions into temporary storage. */
	for (i = start; i <= p_idxs[1]; i++) {
		p = &(pro->parts[i]);
		pro->prev_x[i - start] = p->x;
		pro->prev_y[i - start] = p->y;
		pro->prev_z[i - start] = p->z;
	}

	#ifdef GL
	/* To begin with, none of the particles' trajectories need resetting,
	   so set the entire boolean array of which particles need resetting
	   to false.
	   (Of course the array is of unsigned char as MPI doesn't have a
	   boolean type. Oh well, at least I don't need sizeof here.) */
	memset(pro->traj_reset, 0, ru.track);
	#endif

	/* If testing Rutherford scattering, record the electron's velocity at
	   the beginning of this time step to allow calculation of its
	   deflection angle if it leaves the simulation region. */
	#ifdef DEBUG_RUTHERFORD
	old_v[0] = pro->parts[2].vx;
	old_v[1] = pro->parts[2].vy;
	old_v[2] = pro->parts[2].vz;
	#endif

	/* Start solving each particle's motion for the next time step. */
	for (i = start; i <= p_idxs[1]; i++) {
		p = &(pro->parts[i]);
		dt = ru.dt * pro->tss[p->species];
		#ifdef BORIS
		/* The Boris particle stepping algorithm jumps in with a half-step
		   of free motion before the field is even calculated. */
		p->x += 0.5 * dt * p->vx;
		p->y += 0.5 * dt * p->vy;
		p->z += 0.5 * dt * p->vz;
		#endif
		#ifdef DUST_FIELD_ONLY
		field_particle_feels_from_grain(i, pro, &ru, p->x, p->y, p->z, E);
		#else
		field_particle_feels_at(i, pro, &ru, p->x, p->y, p->z, E, root);
		#endif
		#ifdef DEBUG_FIELDS
		if ((!(time_steps % pro->tss[pro->parts[ION].species])) && i
		    && !((i-1) % ((ru.N - 1) / DEBUG_FIELDS))) {
			/* Assess the accuracy of this E field calculation by emitting
			   the exact E field and estimated E field for this particle. */
			assess_field_accuracy(i, p, E, pro->parts, &ru, time_steps);
		}
		#endif
		ims = i - start;
		if (p->species == DUST) {
			qom = QE * pro->dust_grain_q;
		} else {
			qom = QE * CHARGE[p->species];
		}
		qom /= MASS[p->species];
		#ifdef BORIS
			/* Compute the particle's acceleration under `E` & `B`,
			   separating this out into a half-acceleration, a `B`-based
			   rotation, and a final half-acceleration. Then apply the
			   final free-motion half-step. */
			pro->a_x1[ims] = 0.5 * qom * E[0];
			pro->a_y1[ims] = 0.5 * qom * E[1];
			pro->a_z1[ims] = 0.5 * qom * E[2];
			pro->new_parts[i].vx += dt * pro->a_x1[ims];
			pro->new_parts[i].vy += dt * pro->a_y1[ims];
			pro->new_parts[i].vz += dt * pro->a_z1[ims];
			if (B_mag) {
				/* Do the `B`-based Larmor orbital rotation. With a little
				   luck the compiler's smart enough to optimize this away
				   entirely when there's no magnetic field. */
				/* FIXME: because `dephi` is invariably small, could maybe
			   	   use fast approximations for `dephi`'s sine & cosine? */
				dephi = 2.0 * atanl(0.5 * dt * qom * B_mag);
/*
if (!(i % 99)) {
fprintf(stderr, "%i %lu | %4u %Lg %g | %g %g %g | %g %g %g -> %g", pro->rank, pro->t_steps, i, qom, dephi, pro->new_parts[i].vx, pro->new_parts[i].vy, pro->new_parts[i].vz, B[0], B[1], B[2], B_mag);
}
*/
				rotate_vec_about_dir(b, sin(dephi), cos(dephi), &(pro->new_parts[i].vx), &(pro->new_parts[i].vy), &(pro->new_parts[i].vz), 0.0);
			}
			pro->new_parts[i].vx += dt * pro->a_x1[ims];
			pro->new_parts[i].vy += dt * pro->a_y1[ims];
			pro->new_parts[i].vz += dt * pro->a_z1[ims];
			pro->new_parts[i].x = p->x + (0.5 * dt * pro->new_parts[i].vx);
			pro->new_parts[i].y = p->y + (0.5 * dt * pro->new_parts[i].vy);
			pro->new_parts[i].z = p->z + (0.5 * dt * pro->new_parts[i].vz);
		#else
			/* This velocity half-step is identical for both the
			   Euler-Richardson and velocity Verlet algorithms. */
			pro->a_x1[ims] = qom * (E[0] + p->vy * B[2] - p->vz * B[1]);
			pro->a_y1[ims] = qom * (E[1] + p->vz * B[0] - p->vx * B[2]);
			pro->a_z1[ims] = qom * (E[2] + p->vx * B[1] - p->vy * B[0]);
			pro->new_parts[i].vx += 0.5 * dt * pro->a_x1[ims];
			pro->new_parts[i].vy += 0.5 * dt * pro->a_y1[ims];
			pro->new_parts[i].vz += 0.5 * dt * pro->a_z1[ims];
			/* Now do the algorithm-specific bits of the relevant
			   non-Boris stepping algorithm. */
			#ifdef EULER_RICHARDSON
			pro->new_parts[i].x += 0.5 * dt * p->vx;
			pro->new_parts[i].y += 0.5 * dt * p->vy;
			pro->new_parts[i].z += 0.5 * dt * p->vz;
			#else
			pro->new_parts[i].x += dt * (p->vx + (0.5 * dt * pro->a_x1[ims]));
			pro->new_parts[i].y += dt * (p->vy + (0.5 * dt * pro->a_y1[ims]));
			pro->new_parts[i].z += dt * (p->vz + (0.5 * dt * pro->a_z1[ims]));
			#endif
		#endif  /* end preproc. dir've for Boris vs. non-Boris algorithms */
	}

	#ifndef BORIS
		/* The non-Boris particle stepping algorithms require a second force
		   evaluation halfway into the time step. For these, the master
		   process here collects all of the processes' half-time results
		   for the motion stepping calculation, then broadcasts them all
		   to the workers so they're ready for the second half. */
		if (pro->rank) {
			/* Worker processes send their intermediate particle states back
			   to the master. */
			MPI_Send(&(pro->new_parts[start]), p_idxs[1] - start + 1,
			         MPIParticle, 0, 0, pro->work_comm);
		} else {
			/* Master process overwrites part of its old particle state array
			   with its intermediate results, then overwrites the rest of the
			   old array with the other processes' intermediate results. */
			memcpy(&(pro->parts[p_idxs[0]]),
			       &(pro->new_parts[p_idxs[0]]),
			       sizeof(Particle) * (p_idxs[1] - p_idxs[0] + 1));
			for (i = 1; i < pro->num_workers; i++) {
				/* FIXME: should I be checking the ret. val. of this? */
				MPI_Recv(&(pro->parts[pro->starts[i]]),
/* FIXME
				         pro->pits[i][1] - pro->pits[i][0] + 1, MPIParticle,
*/
				         pro->pits[i][1] - pro->starts[i] + 1, MPIParticle,
				         i, 0, pro->work_comm, MPI_STATUS_IGNORE);
			}
		}
		#ifdef DUST_FIELD_ONLY
		/* Broadcast only the dust grain's state data, as only the grain's
		   electric field/potential is ultimately used in the calculations. */
		if (MPI_Bcast(pro->parts, 1, MPIParticle, 0, pro->work_comm) != MPI_SUCCESS) {
			fprintf(stderr, "%i: failed to send/receive grain state data.\n",
			        pro->rank);
			return;
		}
		#else
		/* Broadcast every particle's information. */
		if (MPI_Bcast(pro->parts, ru.N, MPIParticle, 0, pro->work_comm) != MPI_SUCCESS) {
			fprintf(stderr, "%i: failed to send/receive particle state data.\n",
			        pro->rank);
			return;
		}
		#endif
	#endif  /* end preproc. directive for non-Boris stepping algorithms */

	/* Finish solving particle motion for the next time step (for
	   non-Boris particle motion integrators). Once done, check that
	   all particles remain in the simulation area. */
	for (i = start; i <= p_idxs[1]; i++) {

		p = &(pro->parts[i]);
		p_n = &(pro->new_parts[i]);
		ims = i - start;

		#ifndef BORIS
			/* The unlucky non-Boris algorithms (Euler-Richardson and
			   velocity Verlet) have to carry out a second force evaluation
			   to finish stepping particles for this time step. */
			#ifdef DUST_FIELD_ONLY
			field_particle_feels_from_grain(i, pro, &ru, p->x, p->y, p->z, E);
			#else
			field_particle_feels_at(i, pro, &ru, p->x, p->y, p->z, E, root);
			#endif
			if (p->species == DUST) {
				qom = QE * pro->dust_grain_q;
			} else {
				qom = QE * CHARGE[p->species];
			}
			qom /= MASS[p->species];
			a_x2 = qom * (E[0] + p->vy * B[2] - p->vz * B[1]);
			a_y2 = qom * (E[1] + p->vz * B[0] - p->vx * B[2]);
			a_z2 = qom * (E[2] + p->vx * B[1] - p->vy * B[0]);
			dt = ru.dt * pro->tss[p->species];
			#ifdef EULER_RICHARDSON
			p_n->x = (p->x - 0.5 * dt * (p->vx - 0.5 * dt * pro->a_x1[ims])) + dt * p->vx;
			p_n->y = (p->y - 0.5 * dt * (p->vy - 0.5 * dt * pro->a_y1[ims])) + dt * p->vy;
			p_n->z = (p->z - 0.5 * dt * (p->vz - 0.5 * dt * pro->a_z1[ims])) + dt * p->vz;
			p_n->vx = (p->vx - 0.5 * dt * pro->a_x1[ims]) + dt * a_x2;
			p_n->vy = (p->vy - 0.5 * dt * pro->a_y1[ims]) + dt * a_y2;
			p_n->vz = (p->vz - 0.5 * dt * pro->a_z1[ims]) + dt * a_z2;
			#else
			p_n->vx += 0.5 * dt * a_x2;
			p_n->vy += 0.5 * dt * a_y2;
			p_n->vz += 0.5 * dt * a_z2;
			#endif
		#endif  /* end preproc. dir've for Boris vs. non-Boris algorithms */

		#ifndef DISABLE_GRAIN
		/* If necessary, check whether this particle's been absorbed by
		   the dust grain. */
		if ((i > 0) && ((time_steps * ru.dt) > ru.settle)) {
			/* If this particle's close enough to be absorbed by the
			   dust grain (assumed to be particle 0), it gets absorbed
			   and reinjected. */
			delta_pos[0] = p_n->x - pro->prev_x[ims];
			delta_pos[1] = p_n->y - pro->prev_y[ims];
			delta_pos[2] = p_n->z - pro->prev_z[ims];
			if (does_particle_cross_grain(p_n, delta_pos,
							&(pro->new_parts[0]), ru.a
							#if defined(PROLATE) || defined(OBLATE)
							,ru.A
							#endif
							)) {
				/* If the absorbed particle's an electron, increment the
				   count of absorbed electrons. */
				if (p_n->species == ELECTRON) {
					absorbed_e++;
				}
				/* Add the absorbed particle's mass, charge, and linear
				   momentum to the dust grain's. Exception: if the
				   particle's positive and `POIS_POS_COLLECTION` is
				   defined, don't add the particle's charge to the grain. */
				absorbed_m += MASS[p_n->species];
				absorbed_q += CHARGE[p_n->species];
				#ifdef POIS_POS_COLLECTION
				if (CHARGE[p_n->species] > 0) {
					/* Undo the addition of this positive particle's charge
					   to the grain's charge. */
					absorbed_q -= CHARGE[p_n->species];
				}
				#endif
				dpx += MASS[p_n->species] * p_n->vx;
				dpy += MASS[p_n->species] * p_n->vy;
				dpz += MASS[p_n->species] * p_n->vz;
				/* Reinject this absorbed particle and (if necessary)
				   reset its recorded trajectory. */
				inject(p_n, ru, true, pro);
				#ifdef GL
				if (part_idx_to_traj_idx(&ru, i) != UINT_MAX) {
					pro->traj_reset[part_idx_to_traj_idx(&ru, i)] = 1;
				}
				#endif
			}

		}
		#endif

		#ifdef OUT_BOUNCE
		switch (check_particle_in_sim_area(p_n, ru,
		                                   pro->prev_x[ims],
		                                   pro->prev_y[ims],
		                                   pro->prev_z[ims])) {
		#else
		switch (check_particle_in_sim_area(p_n, ru, pro)) {
		#endif

		/* Here follows a fairly yucky bit of preprocessor footwork. */
		#ifdef GL
		case 1:
			if (part_idx_to_traj_idx(&ru, i) != UINT_MAX) {
				pro->traj_reset[part_idx_to_traj_idx(&ru, i)] = 1;
			}
			#ifdef DEBUG_RUTHERFORD
			if (p_n->species == ELECTRON) {
				printf("%g %g %g\n%g %g %g %g ",
				       old_v[0], old_v[1], old_v[2],
				       rutherford_impact_parameter(&(pro->parts[1]), p_n),
				       p_n->vx, p_n->vy, p_n->vz);
			}
			#endif
		break;
		#else
		case 1:
			#ifdef DEBUG_RUTHERFORD
			if (p_n->species == ELECTRON) {
				printf("%g %g %g\n%g %g %g %g ",
				       old_v[0], old_v[1], old_v[2],
				       rutherford_impact_parameter(&(pro->parts[1]), p_n),
				       p_n->vx, p_n->vy, p_n->vz);
			}
			#endif
		break;
		#endif

		#ifdef SPHERICAL

		case -1:
			fprintf(stderr, "%i: warning: particle %u in outer wall "
			        "travelling inwards.\n", pro->rank, i);
		break;

		#else

		case -1:
			fprintf(stderr, "%i: warning: particle %u in left wall "
			        "travelling rightwards.\n", pro->rank, i);
		break;

		case -2:
			fprintf(stderr, "%i: warning: particle %u in right wall "
			        "travelling leftwards.\n", pro->rank, i);
		break;

		case -3:
			fprintf(stderr, "%i: warning: particle %u in top wall "
			        "travelling downwards.\n", pro->rank, i);
		break;

		case -4:
			fprintf(stderr, "%i: warning: particle %u in bottom wall "
			        "travelling upwards.\n", pro->rank, i);
		break;

		case -5:
			fprintf(stderr, "%i: warning: particle %u in front wall "
			        "travelling towards back wall.\n", pro->rank, i);
		break;

		case -6:
			fprintf(stderr, "%i: warning: particle %u in back wall "
			        "travelling towards front wall.\n", pro->rank, i);
		break;

		#endif

		default:
		break;  /* the particle is within the simulation area */

		}
	}

	if (compute_epe) {

		/* Reset the computed total electrostatic potential energy for
		   this process' electrons and ions. */
		pro->epe[ELECTRON] = 0.0;
		pro->epe[ION] = 0.0;

		/* Keep track of the net electrostatic potential energy of this
		   process' electrons and the ions (so skip the dust grain, if
		   necessary, by ensuring `i` is always >= 1. */
		i = p_idxs[0];
		if (!i) {
			i++;
		}
		for (; i <= p_idxs[1]; i++) {
			p_n = &(pro->new_parts[i]);
			#ifdef DUST_FIELD_ONLY
			pro->epe[p_n->species] += QE * CHARGE[p_n->species]
			                          * phi_particle_feels_from_grain(i, pro,
			                                                          &ru);
			#else
			pro->epe[p_n->species] += QE * CHARGE[p_n->species]
			                          * phi_particle_feels(i, pro, &ru, root);
			#endif
		}

		/* Add up the net potential energy experienced by the electrons
		   and the ions across the worker processes.
		   (FIXME: only the master process really needs this information,
		   at least for now. Though this is simpler code.) */
		MPI_Allreduce(MPI_IN_PLACE, &(pro->epe[1]), 2,
		              MPI_LONG_DOUBLE, MPI_SUM, pro->work_comm);

	}

	#ifndef DISABLE_GRAIN
	/* Every process tells the master process how much net mass, charge &
	   momentum the grain gained this time step, as well as the number of
	   electrons absorbed; the master process adds each of these up,
	   updates its `pro->new_parts` and `pro->total_absorbed_e` accordingly,
	   and recomputes the dust grain's electrostatic potential energy
	   if necessary. */
	MPI_Reduce(&absorbed_m, &net_abs_m, 1, MPI_LONG_DOUBLE,
	           MPI_SUM, 0, pro->work_comm);
	MPI_Reduce(&absorbed_q, &net_abs_q, 1, MPI_INT,
	           MPI_SUM, 0, pro->work_comm);
	MPI_Reduce(&absorbed_e, &net_abs_e, 1, MPI_UNSIGNED,
	           MPI_SUM, 0, pro->work_comm);
	MPI_Reduce(&dpx, &net_dpx, 1, MPI_LONG_DOUBLE,
	           MPI_SUM, 0, pro->work_comm);
	MPI_Reduce(&dpy, &net_dpy, 1, MPI_LONG_DOUBLE,
	           MPI_SUM, 0, pro->work_comm);
	MPI_Reduce(&dpz, &net_dpz, 1, MPI_LONG_DOUBLE,
	           MPI_SUM, 0, pro->work_comm);
	if (!pro->rank) {
		/* Note implicit assumption that it's OK to assume changes to
		   `pro->new_parts[0]` here automatically get propagated through to
		   `pro->parts` by the `memcpy` call below. Also, angular momentum
		   isn't checked here; the dust grain can't start spinning even
		   though really it would be free to rotate in reality. */
		pro->dust_grain_m += net_abs_m;
		pro->new_parts[0].vx += net_dpx / pro->dust_grain_m;
		pro->new_parts[0].vy += net_dpy / pro->dust_grain_m;
		pro->new_parts[0].vz += net_dpz / pro->dust_grain_m;
		pro->dust_grain_q += net_abs_q;
		pro->tot_absor_e += net_abs_e;
		#ifdef POIS_POS_COLLECTION
		/* If the plasma's settling period is complete, add a
		   Poisson-distributed number of elementary positive charges to
		   the grain, with the rate set according to OML theory using the
		   cached ion temperature. That cached temperature is only updated
		   when the system's macroscopic statistics are recalculated, but
		   that shouldn't matter at equilibrium (especially as the ion
		   temperature changes relatively slowly anyway).
		   Warning: this code implicitly assumes there's only one species
		   of positive ion, with a charge of +e. */
		if ((time_steps * ru.dt) > ru.settle) {
			kb_ti = KB * pro->Ti;
			phi_normed = CHARGE[ION] * QE * QE * pro->dust_grain_q;
			phi_normed /= 4.0 * PI * EPSILON0 * ru.a * kb_ti;
			lambda = pro->pos_rate_const * sqrt(kb_ti) * (1.0 - phi_normed);
			pro->dust_grain_q += prv(&(pro->ps), lambda);
		}
		#endif
		if (compute_epe) {
			#ifdef DUST_FIELD_ONLY
			pro->epe[DUST] = 0.0;
			#else
			if (!pro->dust_grain_q) {
				pro->epe[DUST] = 0.0;
			} else {
				pro->epe[DUST] = pro->dust_grain_q * QE
				                 * phi_particle_feels(0, pro, &ru, root);
			}
			#endif
		}
	}
	#endif

	if (pro->rank) {
		/* Worker processes send their stepped particle results back
		   to the master process. */
		MPI_Send(&(pro->new_parts[start]), p_idxs[1] - start + 1,
		         MPIParticle, 0, 0, pro->work_comm);
	} else {
		/* Master process overwrites the start of its old particle state
		   array with the results of its own particle stepping, then
		   overwrites the rest with the other processes' results.
		   (Note implicit assumption here that the master process handles the
		    particles at the start of the particle arrays `pro->parts` &
		    `pro->new_parts`.) */
		memcpy(pro->parts, pro->new_parts,
		       sizeof(Particle) * (p_idxs[1] - p_idxs[0] + 1));
		for (i = 1; i < pro->num_workers; i++) {
/* FIXME
			if (MPI_Recv(&(pro->parts[pro->pits[i][0]]), pro->pits[i][1] - pro->pits[i][0] + 1, MPIParticle, i, 0, pro->work_comm, &recv_status) != MPI_SUCCESS) {
*/
			if (MPI_Recv(&(pro->parts[pro->starts[i]]), pro->pits[i][1] - pro->starts[i] + 1, MPIParticle, i, 0, pro->work_comm, &recv_status) != MPI_SUCCESS) {
				fprintf(stderr, "0: failed to receive particle data from process %i.\n", i);
				/* FIXME: could maybe do something with `recv_status` here? */
			}
		}
	}

	#ifdef DEBUG_KEPLER
	/* If solving a two-body stable Keplerish orbit, periodically record the
	   electron and ion's positions & velocities. */
	if (!(time_steps % 1000000)) {
		printf("%Lg \"i\" %.8g %.8g %.8g %.8g %g %g\n", time_steps * ru.dt,
		       pro->parts[1].x, pro->parts[1].y, pro->parts[1].z,
		       pro->parts[1].vx, pro->parts[1].vy, pro->parts[1].vz);
		printf("%Lg \"e\" %g %g %g %g %g %g\n", time_steps * ru.dt,
		       pro->parts[2].x, pro->parts[2].y, pro->parts[2].z,
		       pro->parts[2].vx, pro->parts[2].vy, pro->parts[2].vz);
	}
	#endif

	#ifndef DISABLE_GRAIN
	/* Master process broadcasts the dust grain's newly-calculated
	   charge & mass to the workers. */
	MPI_Bcast(&(pro->dust_grain_q), 1, MPI_INT, 0, pro->work_comm);
	MPI_Bcast(&(pro->dust_grain_m), 1, MPI_LONG_DOUBLE, 0, pro->work_comm);
	#endif

	#ifdef GL
	/* Each worker's tracked whether each of their particles'
	   trajectories needs resetting (because a particle was reinjected).
	   However, only the master process can use this information to
	   reset the stored trajectories. So collect this information from
	   all of the workers, allowing the master process to reset the
	   necessary trajectories. */
	if ((((unsigned int) pro->rank) < ru.N) && ru.track) {
		/* (FIXME: as above, only the master process needs this
		   information, but again this is simpler code.) */
		MPI_Allreduce(MPI_IN_PLACE, pro->traj_reset, ru.track,
					  MPI_UNSIGNED_CHAR, MPI_LOR, pro->work_comm);
	}
	if (!pro->rank) {
		for (i = 0; i < ru.track; i++) {
			if (pro->traj_reset[i]) {
				reset_trajectory(pro, &ru, i);
			}
		}
	}
	#endif
}

unsigned int set_up_workers_for_assignments(Proc* pro, Run ru)
{
	unsigned int i;
	unsigned int* p;
    MPI_Group work_grp;
    MPI_Group world_grp;
    int worker_range[1][3] = { { 0, 0, 1 } };

	/* If there are more processes on offer to do work than there are
	   particles to simulate(!), then the program only uses as many
	   processes as there are particles to simulate, assigning one particle
	   to each process. Set the number of processes doing work, allowing for
	   this possibility. */
	if (((unsigned int) pro->num_procs) > ru.N) {
		pro->num_workers = ru.N;  /* only `ru.N` processes do any work */
	} else {
		pro->num_workers = pro->num_procs;  /* every process works */
	}

	/* Now that the number of worker processes is decided, set up the
	   workers' MPI communicator. */
	worker_range[0][1] = pro->num_workers - 1;
	MPI_Comm_group(MPI_COMM_WORLD, &world_grp);
	MPI_Group_range_incl(world_grp, 1, worker_range, &work_grp);
	MPI_Comm_create(MPI_COMM_WORLD, work_grp, &(pro->work_comm));

	/* The surplus processes don't need to do anything more now that
	   they've calculated the number of worker processes. */
	if (((unsigned int) pro->rank) >= ru.N) {
		pro->pits = NULL;
		return 0;
	}

	/* Allocate a little memory to store a list of array indices giving the
	   range of particles each process is responsible for stepping. */
	if ((p = malloc(2 * sizeof(unsigned int) * pro->num_workers)) == NULL) {
		fprintf(stderr, "%i: can't allocate %lu bytes of memory for particle array indices.\n", pro->rank, pro->num_workers * 2 * sizeof(unsigned int));
		return 1;
	}
	if ((pro->pits = malloc(sizeof(unsigned int*) * pro->num_workers)) == NULL) {
		fprintf(stderr, "%i: can't allocate %lu bytes of memory for particle array index pointers.\n", pro->rank, pro->num_workers * 2 * sizeof(unsigned int*));
		free(p);
		return 1;
	}

	/* Work out which subset of the particles in the particle array each
	   process will step forward, and make a note in this process'
	   `pro->pits` variable. Basically, the total number of particles is
	   divided by the number of available processes, then rounded down; each
	   process gets this many particles. The remainder are assigned to the
	   last process, so it takes up any slack. */
	for (i = 0; i < pro->num_workers; i++) {
		pro->pits[i] = &(p[i * 2]);
		pro->pits[i][0] = i * (ru.N / pro->num_workers);
		if (i == (pro->num_workers - 1)) {
			pro->pits[i][1] = ru.N - 1;
		} else {
			pro->pits[i][1] = (i + 1) * (ru.N / pro->num_workers) - 1;
		}
	}

	return 0;
}

#ifdef RECORD_PHI_DATA
void write_phi_data_to_file(const char* file_path, const Proc* pro, Run ru, unsigned long int time_step, Tree* root)
{
	double azi;  /* azimuthal angle */
	unsigned int azi_idx;
	FILE* fp = fopen(file_path, "ab");
	double r;
	unsigned int r_idx;
	#ifdef SPHERICAL
	double r_step = ru.r * ru.mpp / CELLS_R;
	#else
	double r_step = 0.5 * ru.L * ru.mpp / CELLS_R;
	#endif
	const double x_grain = pro->parts[0].x;
	const double y_grain = pro->parts[0].y;
	const double z_grain = pro->parts[0].z;
	double zen;  /* zenith angle */
	unsigned int zen_idx;

	if (fp == NULL) {
		fprintf(stderr, "Can't open %s to write electric potential data.\n",
		        file_path);
		return;
	}

	for (r_idx = 1; r_idx <= CELLS_R; r_idx++) {
		for (azi_idx = 0; azi_idx < CELLS_AZI; azi_idx++) {
			for (zen_idx = 0; zen_idx <= CELLS_ZEN; zen_idx++) {
				r = r_step * r_idx;
				azi = 2.0 * PI * azi_idx / CELLS_AZI;
				zen = PI * zen_idx / CELLS_ZEN;
				fprintf(fp, "%Lg %g %g %g %g\n",
				        time_step * ru.dt, r, azi, zen,
				        phi_at(x_grain + (r * cos(azi) * sin(zen)),
				               y_grain + (r * sin(azi) * sin(zen)),
				               z_grain + (r * cos(zen)),
				               pro, &ru, root));
			}
		}
	}

	fclose(fp);
}
#endif

#ifdef RECORD_V_R_DATA
void write_v_r_data_to_file(const char* file_path, const Particle* parts, Run ru, unsigned long int time_step)
{
	double direc[3];
	FILE* fp = fopen(file_path, "ab");
	unsigned int i;
	unsigned int n_e[CELLS_R];
	unsigned int n_i[CELLS_R];
	double r;
	unsigned int r_idx;
	double r_step = ru.r * ru.mpp / CELLS_R;
	double v_r;
	double v_r_e[CELLS_R];
	double v_r_i[CELLS_R];
	const double x_grain = parts[0].x;
	const double y_grain = parts[0].y;
	const double z_grain = parts[0].z;

	/* Reset the particle count bins and mean radial velocity bins to zero
	   (maybe not the most upstanding way to do it, but it's easy). */
	memset(n_e, 0, sizeof(unsigned int) * CELLS_R);
	memset(n_i, 0, sizeof(unsigned int) * CELLS_R);
	memset(v_r_e, 0, sizeof(double) * CELLS_R);
	memset(v_r_i, 0, sizeof(double) * CELLS_R);

	/* Skipping the dust grain, iterate over each particle, work out which
	   spatial bin it's in, then increment the ion/electron count for that
	   bin as appropriate. */
	for (i = 1; i < ru.N; i++) {

		/* Calculate the direction vector from the grain's centre to the
		   particle and use it to compute the distance `r` between them. */
		direc[0] = parts[i].x - x_grain;
		direc[1] = parts[i].y - y_grain;
		direc[2] = parts[i].z - z_grain;
		r = pythag(direc[0], direc[1], direc[2]);
		if (r >= (r_step * CELLS_R)) {
			/* Don't count particles outside the full radius of the counting
			   zone (which should simply be the half of the simulation
			   area's smallest dimension). */
			continue;
		}
		r_idx = r / r_step;

		/* Normalize the direction vector to unity length and take its
		   dot product with the particle's velocity to obtain the particle's
		   radial velocity, adding this radial velocity to the running
		   total (the average being taken a bit later). */
		direc[0] /= r;
		direc[1] /= r;
		direc[2] /= r;
		v_r = direc[0] * parts[i].vx + direc[1] * parts[i].vy
		      + direc[2] * parts[i].vz;

		if (parts[i].species == ELECTRON) {
			v_r_e[r_idx] += v_r;
			n_e[r_idx]++;
		} else {
			v_r_i[r_idx] += v_r;
			n_i[r_idx]++;
		}
	}

	if (fp == NULL) {
		fprintf(stderr, "Can't open %s to write radial velocity data.\n",
		        file_path);
		return;
	}

	for (r_idx = 0; r_idx < CELLS_R; r_idx++) {
		/* Record this cell's average electron & ion radial velocity. */
		fprintf(fp, "%Lg %g ", time_step * ru.dt, r_step * (r_idx + 1));
		if (n_e[r_idx]) {
			fprintf(fp, "%g ", v_r_e[r_idx] / n_e[r_idx]);
		} else {
			fputs("NA ", fp);
		}
		fprintf(fp, "%u ", n_e[r_idx]);
		if (n_i[r_idx]) {
			fprintf(fp, "%g ", v_r_i[r_idx] / n_i[r_idx]);
		} else {
			fputs("NA ", fp);
		}
		fprintf(fp, "%u\n", n_i[r_idx]);
	}

	fclose(fp);
}
#endif

#ifdef RECORD_V_EACH_DATA
void write_v_each_data_to_file(const char* file_path, const Particle* parts, Run ru, unsigned long int time_step)
{
	FILE* fp = fopen(file_path, "ab");
	unsigned int i;

	if (fp == NULL) {
		fprintf(stderr, "Can't open %s to write each particle's velocity.\n",
		        file_path);
		return;
	}

	/* Skipping the dust grain, iterate over each particle, writing its
	   velocity to the given file. */
	for (i = 1; i < ru.N; i++) {
		fprintf(fp, "%Lg %i %g %g %g\n",
		        time_step * ru.dt, parts[i].species,
		        parts[i].vx, parts[i].vy, parts[i].vz);
	}

	fclose(fp);
}
#endif

#ifdef RECORD_N_DATA
void write_n_data_to_file(const char* file_path, Run ru, Proc* pro, unsigned long int time_step)
{
	long double azi;  /* azimuthal angle */
	unsigned int azi_idx;
	FILE* fp = fopen(file_path, "ab");
	unsigned int i;
	long double n_e[CELLS_R][CELLS_AZI][CELLS_ZEN];  /* FIXME: use heap? */
	long double n_i[CELLS_R][CELLS_AZI][CELLS_ZEN];  /* FIXME: use heap? */
	const Particle* parts = pro->parts;
	double r;
	unsigned int r_idx;
	#ifdef SPHERICAL
	double r_step = ru.r * ru.mpp / CELLS_R;
	#else
	double r_step = 0.5 * ru.L * ru.mpp / CELLS_R;
	#endif
	const double x_grain = pro->parts[0].x;
	const double y_grain = pro->parts[0].y;
	long double zen;  /* zenith angle */
	unsigned int zen_idx;
	const double z_grain = pro->parts[0].z;

	/* Reset the particle count bins to zero (in a way that's probably not
	   wholly by the book, but is expedient). */
	memset(n_e, 0, sizeof(long double) * CELLS_R * CELLS_AZI * CELLS_ZEN);
	memset(n_i, 0, sizeof(long double) * CELLS_R * CELLS_AZI * CELLS_ZEN);

	/* Skipping the dust grain, iterate over each particle, work out which
	   spatial bin it's in, then increment the ion/electron count for that
	   bin as appropriate. */
	for (i = 1; i < ru.N; i++) {
		r = pythag(parts[i].x - x_grain, parts[i].y - y_grain,
		           parts[i].z - z_grain);
		zen = acosl((parts[i].z - z_grain) / r);
		azi = atan2l(parts[i].y - y_grain, parts[i].x - x_grain);
		if (azi < 0.0) {
			azi += 2 * PI;
		}
		if (r >= (r_step * CELLS_R)) {
			/* Don't count particles outside the full radius of the counting
			   zone (which should simply be half of the simulation area's
			   smallest dimension). */
			continue;
		}
		r_idx = r / r_step;
		azi_idx = azi / (2.0 * PI / CELLS_AZI);
		zen_idx = zen / (PI / CELLS_ZEN);
		if (parts[i].species == ELECTRON) {
			n_e[r_idx][azi_idx][zen_idx]++;
		} else {
			n_i[r_idx][azi_idx][zen_idx]++;
		}
	}

	if (fp == NULL) {
		fprintf(stderr, "Can't open %s to write particle density data.\n",
		        file_path);
		return;
	}

	for (r_idx = 0; r_idx < CELLS_R; r_idx++) {
		for (azi_idx = 0; azi_idx < CELLS_AZI; azi_idx++) {
			for (zen_idx = 0; zen_idx < CELLS_ZEN; zen_idx++) {
#if 0
				/* Compute this cell's outer radius and its volume. */
				r = r_step * (r_idx + 1);
				volume = (powl(r, 3.0) - powl(r_step * r_idx, 3.0)) / 3.0;
				volume *= (2.0 * PI / CELLS_AZI);
				volume *= cosl((PI / CELLS_ZEN) * zen_idx)
				          - cosl((PI / CELLS_ZEN) * (zen_idx + 1));
#endif
				/* Record this cell's location and particle counts. */
				fprintf(fp, "%Lg %g %g %g %Lg %Lg\n",
				        time_step * ru.dt,
				        r_step * (r_idx + 1),
				        2.0 * PI * (azi_idx + 1) / CELLS_AZI,
				        PI * (zen_idx + 1) / CELLS_ZEN,
				        n_e[r_idx][azi_idx][zen_idx],
				        n_i[r_idx][azi_idx][zen_idx]);
			}
		}
	}

	fclose(fp);
}
#endif

#ifdef RECORD_MACRO_DATA

/* Take an array `cms` of four deviations from the mean raised to the powers
   1 to 4, and divide the 2nd, 3rd & 4th deviations as necessary to make
   them the variance, skewness, and excess kurtosis. */
void divide_and_adjust_central_moments(long double* cms, unsigned int N)
{
	cms[1] /= N;
	cms[2] /= N * powl(cms[1], 1.5);
	cms[3] /= N * cms[1] * cms[1];
	cms[3] -= 3.0;  /* make the kurtosis the excess kurtosis */
}

/* Write macroscopic information about the whole simulated system. */
void write_macroscopic_data_to_file(Run ru, Proc* pro, unsigned long int time_step)
{
	long double* am_mean;
	long double am_e_mean[3] = { 0.0, 0.0, 0.0 };
	long double am_i_mean[3] = { 0.0, 0.0, 0.0 };
	double cent;  /* coordinate of simulation domain centre */
	long double debye;
	long double K = 0.0;
	long double Ks;  /* particle KE, shifted by a sample mean */
	long double lm_e_mean[3] = { 0.0, 0.0, 0.0 };
	long double lm_i_mean[3] = { 0.0, 0.0, 0.0 };

	/* Moment-based temperatures of the electrons & ions. */
	double momt_temp_e = 0.0;
	double momt_temp_i = 0.0;

	/* `mvsk_K_e` contains the mean, variance, skewness, and (excess)
	   kurtosis of the electron kinetic energy distribution. */
	long double mvsk_K_e[4] = { 0.0, 0.0, 0.0, 0.0 };

	/* `mvsk_s_e` contains the mean, variance, skewness, and (excess)
	   kurtosis of the electron speed distribution. */
	long double mvsk_s_e[4] = { 0.0, 0.0, 0.0, 0.0 };

	/* `mvsk_s_e` contains the mean, variance, skewness, and (excess)
	   kurtosis of the ion speed distribution. */
	long double mvsk_s_i[4] = { 0.0, 0.0, 0.0, 0.0 };

	unsigned int N_e = 0;
	unsigned int N_i = 0;
	Particle* p;
	unsigned int part_no;
	long double sum_K_ions = 0.0;
	long double sum_sq_K_elec = 0.0;
	long double sum_sq_K_ions = 0.0;
	long double sum_sq_x_elec = 0.0;
	long double sum_x_elec = 0.0;
	long double v;   /* particle speed */
	long double vs;  /* particle speed, shifted by a sample mean */
	long double* v_mean;
	long double v_e_mean[3] = { 0.0, 0.0, 0.0 };
	long double v_i_mean[3] = { 0.0, 0.0, 0.0 };
	#ifdef SPHERICAL
	const long double volume = (4.0 / 3.0) * PI * powl(ru.r * ru.mpp, 3.0);
	#else
	const long double volume = powl(ru.L * ru.mpp, 3.0);
	#endif

	/* If there are hardly any particles in the simulation, don't bother
	   computing macroscopic quantities; there's almost certainly no point
	   and attempting to calculate things like the skewness & kurtosis of
	   the kinetic energy distribution won't work for very small N.
	   Also, the macroscopic data file won't be opened for writing! */
	if (ru.N < 6) {
		if (time_step == 1) {
			e_puts("Particle number very small, "
			       "forgoing writing macroscopic summary data.");
		}
		return;
	}

	/* It may not have been possible to open the macroscopic data file for
	   writing, so check that it's open. */
	if (pro->fp_macro == NULL) {
		fprintf(stderr, "Macroscopic summary data file not open for "
		        "writing (time step %lu).\n", time_step);
		return;
	}

	/* Count the number of electrons & ions and compute their mean
	   velocity, speed, kinetic energy, angular momentum, and linear mom'm.
	   (Note that the particle number index starts at 1 to exclude the
	   dust grain, which has particle number index 0.) */
	#ifdef SPHERICAL
	cent = ru.r * ru.mpp;
	#else
	cent = ru.L * ru.mpp / 2.0;
	#endif
	for (part_no = 1; part_no < ru.N; part_no++) {
		p = &(pro->parts[part_no]);
		v = pythag(p->vx, p->vy, p->vz);
		#ifndef FPE
		/* Manually check for NaN velocity if floating point exceptions
		   aren't being unmasked/detected. */
		if (isnan(v)) {
			/* Occasionally, if the simulation runs with too long a time
			   step and too many particles, things can go haywire. The first
			   indication of this is normally when NaNs start appearing in
			   the macroscopic summary output file. May as well explicitly
			   check for NaN velocities. */
			fprintf(stderr, "Particle %u (%g, %g, %g) has velocity "
			        "(%g, %g, %g) at time step %lu.\n",
			        part_no, p->x, p->y, p->z, p->vx, p->vy, p->vz,
			        time_step);
		}
		#endif
		if (p->species == ELECTRON) {
			mvsk_K_e[0] += 0.5 * MASS[ELECTRON] * v * v;
			mvsk_s_e[0] += v;
			v_mean = v_e_mean;
			am_mean = am_e_mean;
			N_e++;
		} else {
			mvsk_s_i[0] += v;
			v_mean = v_i_mean;
			am_mean = am_i_mean;
			N_i++;
		}
		v_mean[0] += p->vx;
		v_mean[1] += p->vy;
		v_mean[2] += p->vz;
		am_mean[0] += ((p->y - cent) * p->vz) - ((p->z - cent) * p->vy);
		am_mean[1] -= ((p->x - cent) * p->vz) - ((p->z - cent) * p->vx);
		am_mean[2] += ((p->x - cent) * p->vy) - ((p->y - cent) * p->vx);
	}
	mvsk_K_e[0] /= N_e;
	mvsk_s_e[0] /= N_e;
	mvsk_s_i[0] /= N_i;
	v_e_mean[0] /= N_e;
	v_e_mean[1] /= N_e;
	v_e_mean[2] /= N_e;
	v_i_mean[0] /= N_i;
	v_i_mean[1] /= N_i;
	v_i_mean[2] /= N_i;
	am_e_mean[0] *= MASS[ELECTRON];
	am_e_mean[1] *= MASS[ELECTRON];
	am_e_mean[2] *= MASS[ELECTRON];
	am_i_mean[0] *= MASS[ION];
	am_i_mean[1] *= MASS[ION];
	am_i_mean[2] *= MASS[ION];
	lm_e_mean[0] = MASS[ELECTRON] * v_e_mean[0];
	lm_e_mean[1] = MASS[ELECTRON] * v_e_mean[1];
	lm_e_mean[2] = MASS[ELECTRON] * v_e_mean[2];
	lm_i_mean[0] = MASS[ION] * v_i_mean[0];
	lm_i_mean[1] = MASS[ION] * v_i_mean[1];
	lm_i_mean[2] = MASS[ION] * v_i_mean[2];

	/* Compute total kinetic energy (KE), total square of KE, total
	   electrostatic potential energy (PE), and total square of PE for all
	   particles.
	   FIXME?: the 1-pass calculations of means & standard deviations
	   are a bit numerically unstable.
	   Oh, also, compute the variance, skewness & excess kurtosis of 
	   electron speed (this is pass 2 of a 2-pass algorithm), plus the
	   electron & ion "temperatures" based on the second moment of their
	   velocity distributions. */
	for (part_no = 1; part_no < ru.N; part_no++) {
		p = &(pro->parts[part_no]);
		v = pythag(p->vx, p->vy, p->vz);
		K = 0.5 * MASS[p->species] * v * v;
		if (p->species == ELECTRON) {
			sum_sq_K_elec += K * K;
			sum_x_elec += p->x;
			sum_sq_x_elec += p->x * p->x;
			Ks = K - mvsk_K_e[0];
			mvsk_K_e[1] += Ks * Ks;
			mvsk_K_e[2] += Ks * Ks * Ks;
			mvsk_K_e[3] += Ks * Ks * Ks * Ks;
			vs = v - mvsk_s_e[0];
			mvsk_s_e[1] += vs * vs;
			mvsk_s_e[2] += vs * vs * vs;
			mvsk_s_e[3] += vs * vs * vs * vs;
			momt_temp_e += (p->vx - v_e_mean[0]) * (p->vx - v_e_mean[0]);
			momt_temp_e += (p->vy - v_e_mean[1]) * (p->vy - v_e_mean[1]);
			momt_temp_e += (p->vz - v_e_mean[2]) * (p->vz - v_e_mean[2]);
		} else {
			sum_K_ions += K;
			sum_sq_K_ions += K * K;
			vs = v - mvsk_s_i[0];
			mvsk_s_i[1] += vs * vs;
			mvsk_s_i[2] += vs * vs * vs;
			mvsk_s_i[3] += vs * vs * vs * vs;
			momt_temp_i += (p->vx - v_i_mean[0]) * (p->vx - v_i_mean[0]);
			momt_temp_i += (p->vy - v_i_mean[1]) * (p->vy - v_i_mean[1]);
			momt_temp_i += (p->vz - v_i_mean[2]) * (p->vz - v_i_mean[2]);
		}
	}
	divide_and_adjust_central_moments(mvsk_K_e, N_e);
	divide_and_adjust_central_moments(mvsk_s_e, N_e);
	divide_and_adjust_central_moments(mvsk_s_i, N_i);
	momt_temp_e *= (MASS[ELECTRON] / (3.0 * KB)) / N_e;
	momt_temp_i *= (MASS[ION] / (3.0 * KB)) / N_i;

	#ifdef POIS_POS_COLLECTION
	/* Poisson-distributed ion collection requires an estimate of the ion
	   temperature, so cache this function's calculated ion temperature. */
	pro->Ti = momt_temp_i;
	#endif

	/* Output the time step number, mean KEs & PEs, mean (& s.d. of)
	   electron x-coordinate, the dust grain charge, the (cold ion)
	   Debye length, Debye number, moments of the electron & ion speed
	   distributions, T_e & T_i, and the average components of the
	   electrons' & ions' angular momenta. 
	   (Note that the dust grain's KE is computed using `hypot` instead
	    of `pythag` as underflow could be a threat actually worth
	    worrying about here.) */
	debye = sqrtl(EPSILON0 * KB * momt_temp_e /
	              ((N_e / volume) * QE * QE));
	fprintf(pro->fp_macro, "%.7Lg %Lg ", time_step * ru.dt, mvsk_K_e[0]);
	if (isnt_weird(pro->epe[ELECTRON])) {
		/* FIXME: for some reason SIGFPE gets raised here if the simulation
		   hasn't run for long enough to generate good values for the
		   energies in pro->epe. Why doesn't the if-statement check prevent
		   this?
		   FIXME FIXME: if-statement check might work now...? */
		fprintf(pro->fp_macro, "%Lg ", pro->epe[ELECTRON] / N_e);
	} else {
		fputs("NA ", pro->fp_macro);
	}
	fprintf(pro->fp_macro, "%Lg ", sum_K_ions / N_i);
	if (isnt_weird(pro->epe[ION])) {
		fprintf(pro->fp_macro, "%Lg ", pro->epe[ION] / N_i);
	} else {
		fputs("NA ", pro->fp_macro);
	}
	fprintf(pro->fp_macro, "%Lg ", 
	        0.5 * pro->dust_grain_m * hypot(pro->parts[0].vx,
	                                        hypot(pro->parts[0].vy,
	                                              pro->parts[0].vz)));
	if (isnt_weird(pro->epe[DUST])) {
		fprintf(pro->fp_macro, "%Lg ", pro->epe[DUST]);
	} else {
		fputs("NA ", pro->fp_macro);
	}
	fprintf(pro->fp_macro, "%.7Lg %.7Lg %i %Lg %Lg"
	        " %Lg %Lg %Lg %Lg %Lg %Lg %Lg %Lg %g %g %Lg %Lg %Lg"
	        " %Lg %Lg %Lg %Lg %Lg %Lg %Lg %Lg %Lg %u\n",
	        sum_x_elec / N_e,
	        sqrtl(sum_sq_x_elec / N_e
			      - (sum_x_elec / N_e) * sum_x_elec / N_e),
	        pro->dust_grain_q,
	        debye,
	        (4.0 * PI / 3.0) * debye * debye * debye * N_e / volume,
	        mvsk_s_e[0], mvsk_s_e[1], mvsk_s_e[2], mvsk_s_e[3],
	        mvsk_s_i[0], mvsk_s_i[1], mvsk_s_i[2], mvsk_s_i[3],
	        momt_temp_e, momt_temp_i,
	        lm_e_mean[0], lm_e_mean[1], lm_e_mean[2],
	        lm_i_mean[0], lm_i_mean[1], lm_i_mean[2],
	        am_e_mean[0], am_e_mean[1], am_e_mean[2],
	        am_i_mean[0], am_i_mean[1], am_i_mean[2],
	        pro->tot_absor_e);

	/* Try to ensure that the data just written are actually written! */
	if (fflush(pro->fp_macro)) {
		e_puts("Can't force write of macroscopic summary data to disk.");
	}
}

#endif

void write_run_info_to_file(const char* file_path, Run ru)
{
	char date_time[50];
	FILE* fp = fopen(file_path, "ab");
	bool record_macro_data;
	time_t t_t = time(NULL);
	struct tm* t_tm;

	#ifdef RECORD_MACRO_DATA
	record_macro_data = true;
	#else
	record_macro_data = false;
	#endif

	if (fp == NULL) {
		fprintf(stderr,
		        "Can't open %s to write information about this run.\n",
		        file_path);
		return;
	}

	if ((t_tm = localtime(&t_t)) == NULL) {
		fprintf(stderr, "Can't convert Unix time %llu to time structure.\n",
		        (unsigned long long int) t_t);
		return;
	}

	if (!strftime(date_time, sizeof(date_time), "%Y-%m-%d %H:%M:%S", t_tm)) {
		e_puts("Can't represent date & time as string of text.");
		return;
	}

	fprintf(fp,
	        "%s\t%lu\t%Lg\t%g\t%g\t%g\t%g\t%g\t"
	        "%s\t%s\t%u\t"
	        "%g\t%u\t%g\n",
	        date_time, ru.N, ru.dt, ru.soft, ru.mpp,
	        ru.Ti, ru.Te, ru.x_drift,
	        #ifdef SPHERICAL
	        __DATE__, __TIME__, ru.r,
	        #else
	        __DATE__, __TIME__, ru.L,
	        #endif
	        ru.a, (unsigned int) record_macro_data, OAP);

	fclose(fp);
}

double early_oml_charge_time(Run ru, double n_total)
{
	double nu = sqrt(MASS[ION] * ru.Te / (MASS[ELECTRON] * ru.Ti));
	double tim;

	/* This is the naive implicit OML time constant. */
	tim = EPSILON0 * sqrt(2.0 * PI * MASS[ION] * KB
	                      * ru.Te * ru.Te / ru.Ti);
	tim /= (QE * QE * ru.a * n_total / 2.0);

	/* Incorporate the correction factor to give the complete time
	   constant for OML with small t. */
	tim *= log((1.0 - nu) / (pow(nu, 1.0 - exp(-1.0)) - nu));

	return tim;
}

/* Compute the plasma's total particle density (i.e. the density of both
   ions and electrons, not either).
   When computing the particle density `n` the dust grain doesn't count,
   which is why `n` is proportional to (ru.N - 1), not ru.N. */
long double total_number_density(Run ru)
{
	#ifdef SPHERICAL
	return (ru.N - 1) / ((4.0 / 3.0) * PI * powl(ru.r * ru.mpp, 3.0));
	#else
	return (ru.N - 1) / powl(ru.L * ru.mpp, 3.0);
	#endif
}

void write_run_description_to_file(const char* file_path, Run ru)
{
	double ccp;             /* estimated Coulomb coupling parameter */
	double charge_osc_rat;  /* ratio of charging time to elec. osc. time */
	double coulomb_log;     /* estimated Coulomb logarithm */
	long double debye_len;  /* Debye length (m) */
	double debye_number;    /* mean no. of electrons in Debye sphere */
	long double e_colll;    /* estimated electron collision length */
	long double edct;       /* estimated electron dust grain crossing time */
	long double eop;        /* electron plasma oscillation period (s) */
#if 0
	long double ett;        /* estimated electron traversal time */
#endif
	FILE* fp = fopen(file_path, "ab");
	long double htcd;  /* head-on thermal collision distance */
	long double n;     /* mean particle density */
	long double size;  /* characteristic plasma size/smallest dim'n (m) */
	long double wsr;   /* Wigner-Seitz radius */

	if (fp == NULL) {
		fprintf(stderr,
		        "Can't open %s to write this run's initial description.\n",
		        file_path);
		return;
	}

	/* List the runtime settings, marking those the user's set to
	   non-default values with an asterisk.
	   I've never been a fan of the ternary operator but it saves
	   quite a bit of code here! */
	fprintf(fp,
	        "%c total particle count:          %lu\n"
	        "%c time step (s):                 %Lg\n"
	        "%c nominal elec. temperature (K): %g\n"
	        "%c nominal ion temperature (K):   %g\n"
	        "%c softening parameter (m):       %g\n"
	        "%c drift velocity (m*s^-1):       %g\n"
	        "%c pixel size (m):                %g\n"
	        #ifdef SPHERICAL
	        "%c simulation radius (px.):       %u\n"
	        #else
	        "%c simulation length (px.):       %u\n"
	        #endif
	        "%c tracked particles:             %hu\n"
	        "%c number of time steps:          %lu\n"
	        "%c grain radius (m):              %g\n"
		#if defined(PROLATE) || defined(OBLATE)
		"%c aspect ratio:                  %g\n"
		#endif
	        "%c settling time (s):             %g\n\n",
	        (ru.N != DEFAULT_N)             ? '*' : ' ', ru.N,
	        (ru.dt != DEFAULT_DT)           ? '*' : ' ', ru.dt,
	        (ru.Te != DEFAULT_T_E)          ? '*' : ' ', ru.Te,
	        (ru.Ti != DEFAULT_T_I)          ? '*' : ' ', ru.Ti,
	        (ru.soft != DEFAULT_SOFT)       ? '*' : ' ', ru.soft,
	        (ru.x_drift != DEFAULT_X_DRIFT) ? '*' : ' ', ru.x_drift,
	        (ru.mpp != DEFAULT_MPP)         ? '*' : ' ', ru.mpp,
	        #ifdef SPHERICAL
	        (ru.r != DEFAULT_R)             ? '*' : ' ', ru.r,
	        #else
	        (ru.L != DEFAULT_L)             ? '*' : ' ', ru.L,
	        #endif
	        #ifdef GL
	        (ru.track != DEFAULT_TRACK)     ? '*' : ' ', ru.track,
	        #else
	        ' ', 0,
	        #endif
	        (ru.run_for != DEFAULT_RUN_FOR) ? '*' : ' ', ru.run_for,
	        (ru.a != DEFAULT_A)             ? '*' : ' ', ru.a,
		#if defined(PROLATE) || defined(OBLATE)
		(ru.A != DEFAULT_ASPECT_RATIO)  ? '*' : ' ', ru.A,
		#endif
	        (ru.settle != DEFAULT_SETTLE)   ? '*' : ' ', ru.settle);

	/* List estimates of diagnostic plasma parameters.
	   The Coulomb coupling parameter expression comes from U. Schumacher's
	   2005 book chapter "Basics of Plasma Physics". Coulomb logarithm
	   expression's from Susanne Pfalzner's introductory ICF book.
	   (Yep, some more ternary operators here.) */
	n = total_number_density(ru);
	#ifdef SPHERICAL
	size = 2.0 * ru.r * ru.mpp;
	#else
	size = ru.L * ru.mpp;
	#endif
	debye_len = sqrtl(EPSILON0 * KB * ru.Te / ((n / 2.0) * QE * QE));
	debye_number = (4.0 / 3.0) * PI * powl(debye_len, 3.0) * n / 2.0;
	wsr = powl(4.0 * PI * n / 3.0, -1.0 / 3.0);
	if (fabs(CHARGE[ELECTRON]) > fabs(CHARGE[ION])) {
		ccp = QE * CHARGE[ELECTRON] * QE * CHARGE[ELECTRON];
	} else {
		ccp = QE * CHARGE[ION] * QE * CHARGE[ION];
	}
	ccp /= (4.0 * PI * EPSILON0 * wsr * KB * ru.Te);
	coulomb_log = log(9.0 * debye_number);
	e_colll = 16.0 * PI * (n / 2.0) * powl(debye_len, 4.0) / coulomb_log;
	edct = ru.a / sqrtl(3.0 * KB * ru.Te / MASS[ELECTRON]);
	eop = 2.0 * PI * sqrtl(ME * EPSILON0 / (n * QE * QE / 2.0));
	charge_osc_rat = early_oml_charge_time(ru, n) / eop;
	htcd = QE * QE / (12 * PI * EPSILON0 * KB);
	if (ru.Te > ru.Ti) {
		htcd /= ru.Te;
	} else {
		htcd /= ru.Ti;
	}
	fprintf(fp,
	        "  particle density (m^-3):    %Lg\n"
	        "  Debye length (m):           %Lg\n"
	        "  grain radius/DL:            %Lg\n"
	        "  elec. oscil. period (s):    %Lg\n"
	        "  early OML chg. time (s):    %g\n"
	        "  Wigner-Seitz radius (m):    %Lg\n"
	        "  normalized drift speed:     %Lg\n"
	        "  ion traversal time (s):     %Lg\n"
	        "  Coulomb logarithm:          %g\n"
	        "  electron coll. len. (m):    %Lg\n"
	        "%c electron Knudsen number:    %Lg\n"
	        "%c Debye number:               %g\n"
	        "%c Coulomb coupling param.:    %g\n"
	        "%c softening parameter/WSR:    %Lg\n"
	        "%c t. step/elec. cross time:   %Lg\n"
	        "%c DL/characteristic size:     %Lg\n"
	        "%c OML time/elec. osc. period: %g\n"
	        "%c soft. par./head coll. dist: %Lg\n"
	        "%c soft. par./grain radius:    %g\n\n",
	        n, debye_len, ru.a / debye_len, eop,
	        early_oml_charge_time(ru, n), wsr,
	        ru.x_drift / sqrtl(KB * ru.Ti / MASS[ION]),
	        size / sqrtl(3.0 * KB * ru.Ti / MASS[ION]),
	        coulomb_log, e_colll,
	        ((e_colll / size) > 0.5)   ? '!' : ' ', e_colll / size,
	        (debye_number < 2.0)       ? '!' : ' ', debye_number,
	        (ccp > 1.0)                ? '!' : ' ', ccp,
	        ((ru.soft / wsr) > 1.0)    ? '!' : ' ', ru.soft / wsr,
	        ((ru.dt / edct) > 0.5)     ? '!' : ' ', ru.dt / edct,
	        ((debye_len / size) > 0.1) ? '!' : ' ', debye_len / size,
	        (charge_osc_rat < 100.0)   ? '!' : ' ', charge_osc_rat,
	        ((ru.soft / htcd) > 0.5)   ? '!' : ' ', ru.soft / htcd,
	        ((ru.soft / ru.a) > 0.1)   ? '!' : ' ', ru.soft / ru.a);

	/* List more esoteric compile-time simulation parameters. */
	fprintf(fp,
	        "opening angle parameter: %g\n"
	        "ion/electron mass ratio: %Lg\n"
	        "external B field (T):    (%g, %g, %g)\n"
	        "spatial bin counts:      (%i, %i, %i)\n"
	        "max. tree walk length:   %i\n"
	        "grain mass (kg):         %Lg\n\n",
	        OAP, MASS[ION] / MASS[ELECTRON],
	        EXTERNAL_B_X, EXTERNAL_B_Y, EXTERNAL_B_Z,
	        CELLS_R, CELLS_AZI, CELLS_ZEN, MAX_PATH_LEN, MASS[DUST]);

	/* Finally, list the compile-time settings related to binary options,
	   which directory the output files are written to, and the compilation
	   date & time. */

	#ifdef DISABLE_GRAIN
	fputs("grain:                disabled\n", fp);
	#else
	fputs("grain:                enabled\n", fp);
	#endif

	#ifdef GL
	fputs("OpenGL interface:     yes\n", fp);
	fprintf(fp, "trajectory length:    %i data points\n", TRAJ_LEN);
	#else
	fputs("OpenGL interface:     no\n", fp);
	fputs("trajectory length:    N/A\n", fp);
	#endif

	#ifdef OUT_BOUNCE
	fprintf(fp, "containment:          bouncing (c.o.e. %g)\n", OUT_BOUNCE);
	#else
	#ifdef SCEPTIC_REINJECTION
	fputs("containment:          reinjection (SCEPTIC-style)\n", fp);
	#else
	fputs("containment:          reinjection (old-style)\n", fp);
	#endif
	#endif

	#ifdef FPE
	fputs("FP exceptions:        yes\n", fp);
	#else
	fputs("FP exceptions:        no\n", fp);
	#endif

	#ifdef DIPOLES
	fputs("multipole expansion:  dipole\n", fp);
	#else
	fputs("multipole expansion:  monopole\n", fp);
	#endif

	#ifdef SPHERICAL
	fputs("simulation shape:     spherical\n", fp);
	#else
	fputs("simulation shape:     cubic\n", fp);
	#endif

	#ifdef EULER_RICHARDSON
	fputs("numerical integrator: Euler-Richardson\n", fp);
	#else
	#ifdef BORIS
	fputs("numerical integrator: Boris\n", fp);
	#else
	fputs("numerical integrator: velocity Verlet\n", fp);
	#endif
	#endif

	#ifdef POIS_POS_COLLECTION
	fputs("pos. Pois. collec'n:  yes\n", fp);
	#else
	fputs("pos. Pois. collec'n:  no\n", fp);
	#endif

	#ifdef DEBUG_LANGMUIR
	fputs("Langmuir test wave:   yes\n", fp);
	#else
	fputs("Langmuir test wave:   no\n", fp);
	#endif

	#ifdef DUST_FIELD_ONLY
	fputs("P-P interactions:     no\n", fp);
	#else
	fputs("P-P interactions:     yes\n", fp);
	#endif

	#ifdef PROLATE
	fputs("Grain shape:          prolate spheroid\n\n", fp);
	#elif defined(OBLATE)
	fputs("Grain shape:          oblate spheroid\n\n", fp);
	#else
	fputs("Grain shape:          sphere\n\n", fp);
	#endif

	fprintf(fp,
	        "output file dir.:     %s\n"
	        "compiled on:          %s %s\n",
	        OUTPUT_DIR, __DATE__, __TIME__);

	fclose(fp);
}

/* Write complete state information for the simulation's particles. */
void write_sim_state_to_file(unsigned int N, const Proc* pro, const char* path)
{
	double chi_b;  /* normalized electric edge potential */
	FILE* fp;
	unsigned char good_to_go = false;  /* `bool` sadly infeasible here */
	MPI_Datatype MPIPRNG;
	PRNG other_ps;
	int ran;

	#ifdef SCEPTIC_REINJECTION
	chi_b = pro->chi_b;
	#else
	chi_b = 9999.99;
	#endif

	MPIPRNG = prng_mpi_type();

	if (pro->rank) {

		if (pro->rank >= (int) pro->num_workers) {
			/* This process isn't simulating anything, so it doesn't have
			   any PRNG state to send to the master process. (Such processes
			   should already have excused themselves by exiting, so this
			   bit of code should never run, but just in case....) */
			MPI_Type_free(&MPIPRNG);
			return;
		}

		/* Get a signal from the master (rank 0) process to decide
		   whether to send this worker process' master state to
		   the master process. `good_to_go` would ideally be a `bool`
		   rather than an `unsigned char` but my MPI installation
		   doesn't define MPI_C_BOOL! */
		MPI_Bcast(&good_to_go, 1, MPI_UNSIGNED_CHAR, 0, pro->work_comm);
		if (good_to_go) {
			MPI_Send((PRNG*) &(pro->ps), 1, MPIPRNG, 0, 1, pro->work_comm);
		}

	} else {

		/* This process is the master process, and is the only process
		   that actually writes the simulation state to a file. */

		good_to_go = (fp = fopen(path, "wb")) != NULL;
		MPI_Bcast(&good_to_go, 1, MPI_UNSIGNED_CHAR, 0, pro->work_comm);

		if (!good_to_go) {
			fprintf(stderr,
			        "Failed to open simulation state file %s for writing.\n",
		    	    path);
			MPI_Type_free(&MPIPRNG);
			return;
		}

		fprintf(fp, "C %u %i %Lg %lu %Lg %Lg %Lg %u %g\n",
		        pro->num_workers, pro->dust_grain_q, pro->dust_grain_m,
		        pro->t_steps - 1, pro->epe[0], pro->epe[1], pro->epe[2],
		        pro->tot_absor_e, chi_b);

		if (fwrite(&(pro->ps), sizeof(PRNG), 1, fp) != 1) {
			fprintf(stderr, "Failed to write proc. 0's PRNG state to %s.\n",
			        path);
		}
		for (ran = 1; ran < (int) pro->num_workers; ran++) {
			MPI_Recv(&other_ps, 1, MPIPRNG, ran, 1,
			         pro->work_comm, MPI_STATUS_IGNORE);
			if (fwrite(&(other_ps), sizeof(PRNG), 1, fp) != 1) {
				fprintf(stderr,
				        "Failed to write proc. %i's PRNG state to %s.\n",
				        ran, path);
			}
		}
		if (fwrite(pro->parts, sizeof(Particle), N, fp) != N) {
			fprintf(stderr, "Failed to write particle states to %s.\n", path);
		}

		if (fclose(fp)) {
			fprintf(stderr,
			        "Failed to close simulation state output file %s.\n",
			        path);
		}

	}

	MPI_Type_free(&MPIPRNG);
}

/* Load the simulation state from a file. */
unsigned int read_sim_state_from_file(unsigned int N, Proc* pro, const char* path)
{
	unsigned int alleged_num_workers;
	double chi_b;
	int first_char;
	FILE* fp = fopen(path, "rb");
	int ran;
	unsigned int successes = 1;
	PRNG throwaway;

	if (fp == NULL) {
		fprintf(stderr,
		        "%i: failed to open simulation state file %s to read.\n",
		        pro->rank, path);
		successes = 0;
	} else {
		first_char = getc(fp);
		if (first_char == 'C') {
			/* This purports to be a format C state file. */
			if ((fscanf(fp, " %u %i %Lg %lu %Lg %Lg %Lg %u %lg\n",
			            &alleged_num_workers, &(pro->dust_grain_q),
			            &(pro->dust_grain_m), &(pro->t_steps),
			            &(pro->epe[0]), &(pro->epe[1]), &(pro->epe[2]),
			            &(pro->tot_absor_e), &(chi_b))) != 9) {
				fprintf(stderr,
				        "%i: failed to read 9 assorted data values "
				        "from %s.\n",
				        pro->rank, path);
				successes = 0;
			}
			#ifdef SCEPTIC_REINJECTION
			pro->chi_b = chi_b;
			#endif
		} else if (first_char == 'B') {
			/* This purports to be a format B state file. */
			if ((fscanf(fp, " %u %i %lu %Lg %Lg %Lg %u\n",
			            &alleged_num_workers, &(pro->dust_grain_q),
			            &(pro->t_steps),
			            &(pro->epe[0]), &(pro->epe[1]), &(pro->epe[2]),
			            &(pro->tot_absor_e))) != 7) {
				fprintf(stderr,
				        "%i: failed to read 7 assorted data values "
				        "from %s.\n",
				        pro->rank, path);
				successes = 0;
			}
		} else if ((first_char >= '0') && (first_char <= '9')) {
			/* This should be a format A state file. */
			rewind(fp);
			if ((fscanf(fp, "%u %i %lu %Lg %Lg %Lg\n", &alleged_num_workers,
			            &(pro->dust_grain_q), &(pro->t_steps), &(pro->epe[0]),
			            &(pro->epe[1]), &(pro->epe[2]))) != 6) {
				fprintf(stderr,
				        "%i: failed to read 6 assorted data values "
				        "from %s.\n",
				        pro->rank, path);
				successes = 0;
			}
		} else {
			/* This isn't actually a state file. */
			fprintf(stderr, "%i: alleged state file %s isn't.\n",
			        pro->rank, path);
			successes = 0;
		}
		if (successes && (alleged_num_workers != pro->num_workers)) {
			fprintf(stderr,
		    	    "%i: conflicting worker process counts (%u vs. %u).\n",
			        pro->rank, alleged_num_workers, pro->num_workers);
			successes = 0;
		}
	}

	if (successes) {  /* if everything went OK so far... */

		for (ran = 0; ran < (int) pro->num_workers; ran++) {
			if (fread(&throwaway, sizeof(PRNG), 1, fp) != 1) {
				fprintf(stderr,
						"%i: failed to read proc. %i's PRNG state from %s.\n",
						pro->rank, ran, path);
				successes = 0;
				break;
			}
			if (ran == pro->rank) {
				/* The PRNG state just read is the one for this process, so
				   stash it in `pro->ps` instead of just throwing it away.
				   Then check whether that state seems invalid, and if so,
				   sanitize it so it doesn't trigger something nasty like a
				   segmentation fault. */
				copy_prng(&(pro->ps), &throwaway);
				if (sanitize_prng(&(pro->ps))) {
					fprintf(stderr, "%i: read invalid PRNG state index, "
							"attempting to sanitize it.\n", ran);
				}
			}
		}
		if (fread(pro->parts, sizeof(Particle), N, fp) != N) {
			fprintf(stderr, "%i: failed to read particle states from %s.\n",
					pro->rank, path);
			successes = 0;
		}

		if (fclose(fp)) {
			fprintf(stderr, "%i: failed to close simulation state file %s.\n",
					pro->rank, path);
			successes = 0;
		}

	}

	/* Add up the number of worker processes which successfully read a
	   simulation state from the requested file. Return 0 if all
	   succeeded, otherwise return 1. */
	MPI_Allreduce(MPI_IN_PLACE, &successes, 1, MPI_UNSIGNED,
	              MPI_SUM, pro->work_comm);
	return (successes != pro->num_workers);
}

#ifdef GL
void run(GLUquadric* qr, Proc* pro, Run ru)
#else
void run(Proc* pro, Run ru)
#endif
{
	#ifdef GL

	/* The GL(U) viewpoint location (where the camera is, effectively). */
	#ifdef SPHERICAL
	GLdouble cam[3] = { 1.5 * ru.r, 0.0, PI / 2.0 };
	#else
	GLdouble cam[3] = { 1.4 * ru.L, -PI / 2.0, PI };
	#endif

	/* Pointers to the 3 components of a tracked particle's last recorded
	   point in its trajectory. */
	float* l_t[3];

	/* Particle array index corresponding to a given trajectory index. */
	unsigned int p_idx;

	/* The "up" vector for `gluLookAt`. */
	GLdouble up[3];

	#endif

	unsigned int idx;
	MPI_Datatype MPIParticle;
	int MPIParticle_block_lens[] = { 1, 1, 1, 1, 1, 1, 1 };
	MPI_Aint MPIParticle_disps[7];  /* these are calculated below */
	MPI_Datatype MPIParticle_subtypes[] = {  /* must match Particle defn.! */
		MPI_DOUBLE,
		MPI_DOUBLE,
		MPI_DOUBLE,
		MPI_DOUBLE,
		MPI_DOUBLE,
		MPI_DOUBLE,
		MPI_INT
	};
	unsigned char master_proc_orders_a_halt = 0;
	unsigned int num_charges;  /* no. p'cles this proc's responsible for */
	bool record_macro_this_step = false;
	unsigned int start_time = time(NULL);

	/* The 1st process writes information about this run to files. */
	if (!pro->rank) {
		write_run_info_to_file(OUTPUT_DIR "pot-run-init-info.txt", ru);
		write_run_description_to_file(OUTPUT_DIR "pot-run-descrip.txt", ru);
	}

	if (((unsigned int) pro->rank) < ru.N) {

		/* Allocate memory for each particle's mid-time step acceleration. */
		num_charges = pro->pits[pro->rank][1] - pro->pits[pro->rank][0] + 1;
		pro->a_x1 = malloc(3 * sizeof(long double) * num_charges);
		if (pro->a_x1 == NULL) {
			fprintf(stderr, "%i: can't allocate %.1f MiB for mid-time step particle accelerations.\n", pro->rank, 3 * sizeof(long double) * num_charges / 1048576.0);
			abort();
		}
		pro->a_y1 = &(pro->a_x1[num_charges]);
		pro->a_z1 = &(pro->a_x1[2 * num_charges]);
	

		/* Allocate memory for particles' positions at the end of the
		   previous time step. */
		pro->prev_x = malloc(3 * sizeof(double) * num_charges);
		if (pro->prev_x == NULL) {
			fprintf(stderr, "%i: can't allocate %.1f MiB for particles' previous positions.\n", pro->rank, 3 * sizeof(double) * num_charges / 1048576.0);
			abort();
		}
		pro->prev_y = &(pro->prev_x[num_charges]);
		pro->prev_z = &(pro->prev_x[2 * num_charges]);

		/* Allocate memory for 9N nodes in the particle tree.
		   The multiplier of 9 is somewhat arbitrary, but if it's too
		   large the code needlessly allocates too much memory and if it's
		   too small the code wastes time having to reallocate the tree
		   pointer with more memory later. */
		pro->tr_alloc = 9 * ru.N;
		if ((pro->root = malloc(sizeof(Tree) * pro->tr_alloc)) == NULL) {
			fprintf(stderr,
			        "%i: can't allocate %lu bytes for particle tree.\n",
			        pro->rank, sizeof(Tree) * pro->tr_alloc);
			abort();
		}

	}

	/* Decide how many time step steps to use for each species, i.e. the
	   multiple of the time step count on which to actually advance each
	   particle species. Heavier species needn't to be stepped as often,
	   and that should be taken advantage of to reduce CPU work. */
	pro->tss[ELECTRON] = 1;
	pro->tss[ION] = (int) sqrt(ru.Te * MASS[ION] / (ru.Ti * MASS[ELECTRON]));
	pro->tss[DUST] = pro->tss[ION] * (int) sqrt(0.1 * MASS[DUST] / MASS[ION]);

	/* Make an MPI derived datatype for transferring a single Particle's
	   state across processes. */
	MPI_Get_address(&(pro->parts[0].x), &(MPIParticle_disps[0]));
	MPI_Get_address(&(pro->parts[0].y), &(MPIParticle_disps[1]));
	MPI_Get_address(&(pro->parts[0].z), &(MPIParticle_disps[2]));
	MPI_Get_address(&(pro->parts[0].vx), &(MPIParticle_disps[3]));
	MPI_Get_address(&(pro->parts[0].vy), &(MPIParticle_disps[4]));
	MPI_Get_address(&(pro->parts[0].vz), &(MPIParticle_disps[5]));
	MPI_Get_address(&(pro->parts[0].species), &(MPIParticle_disps[6]));
	for (idx = 1; idx < 7; idx++) {
		MPIParticle_disps[idx] -= MPIParticle_disps[0];
	}
	MPIParticle_disps[0] = 0;
	MPI_Type_create_struct(7, MPIParticle_block_lens, MPIParticle_disps,
	                       MPIParticle_subtypes, &MPIParticle);
	MPI_Type_commit(&MPIParticle);

	#ifdef GL
	/* If the user wants particle trajectory tracks, allocate the memory to
	   store the trajectories and reset all of their coordinates to zero. */
	if ((!pro->rank) && ru.track) {

		pro->traj_x = malloc(3 * sizeof(float*) * ru.track);
		if (pro->traj_x == NULL) {
			fprintf(stderr, "Can't allocate %lu bytes for pointers to "
			        "particle trajectories.\n",
			        3 * sizeof(float*) * ru.track);
			abort();
		}
		pro->traj_y = &(pro->traj_x[ru.track]);
		pro->traj_z = &(pro->traj_x[2 * ru.track]);

		pro->traj_x[0] = malloc(3 * sizeof(float) * ru.track * TRAJ_LEN);
		if (pro->traj_x[0] == NULL) {
			fprintf(stderr, "Can't allocate %lu bytes for particle "
			        "trajectories.\n",
			        sizeof(float) * ru.track * TRAJ_LEN);
			abort();
		}
		for (idx = 0; idx < ru.track; idx++) {
			pro->traj_x[idx] = &(pro->traj_x[0][idx * TRAJ_LEN]);
			pro->traj_y[idx] = &(pro->traj_x[0][(ru.track + idx) * TRAJ_LEN]);
			pro->traj_z[idx] = &(pro->traj_x[0][((2 * ru.track) + idx) * TRAJ_LEN]);
		}

		for (idx = 0; idx < ru.track; idx++) {
			reset_trajectory(pro, &ru, idx);
		}
	}
	if ((((unsigned int) pro->rank) < ru.N) && ru.track) {
		if ((pro->traj_reset = malloc(ru.track)) == NULL) {
			fprintf(stderr, "Can't allocate %hu bytes for particle "
			        "trajectory reset statuses.\n", ru.track);
			abort();
		}
	}
	#endif

	/* During & after numeric integration to solve particle motion, the
	   1st process has to receive updated particle state data from the
	   other processes, collect it, and then re-transmit it among the
	   other processes. To do so it needs to know the number of particles
	   each process has updated during a time step, which it does based on
	   starting indices (for worker processes' particle arrays) stored in
	   `pro->starts`. So: allocate the memory for that variable. */
	if (!pro->rank) {
		if (((unsigned int) pro->num_procs) >= ru.N) {
			pro->starts = malloc(sizeof(unsigned int) * ru.N);
			if (pro->starts == NULL) {
				fprintf(stderr, "Can't allocate %lu bytes for updated "
				        "particle indices.\n", sizeof(unsigned int) * ru.N);
				abort();
			}
		} else {
			pro->starts = malloc(sizeof(unsigned int) * pro->num_procs);
			if (pro->starts == NULL) {
				fprintf(stderr, "Can't allocate %lu bytes for updated "
				        "particle indices.\n",
				        sizeof(unsigned int) * pro->num_procs);
				abort();
			}
		}
	}

	/* If the collection of positive ions is to be simulated as a stand-
	   alone Poisson process, an initial estimate of the ion temperature
	   is needed to determine the Poisson process' rate. Use the ion
	   reinjection temperature as that initial ion temperature estimate.
	   Then compute and store the constant part of the Poisson process'
	   rate. */
	#ifdef POIS_POS_COLLECTION
	pro->Ti = ru.Ti;
	pro->pos_rate_const = 4.0 * PI * ru.a * ru.a * CHARGE[ION];
	pro->pos_rate_const /= sqrt(2.0 * PI * MASS[ION]);
	pro->pos_rate_const *= ru.dt * total_number_density(ru) / 2.0;
	#endif

	/* Main loop. */
	while (1) {

		pro->t_steps++;
		if ((ru.run_for) && (pro->t_steps > ru.run_for)) {
			/* Enough time steps have been simulated, stop running. */
			break;
		}

		if (ru.tim_lim > 0.0) {
			if (time(NULL) < start_time) {
				fprintf(stderr,
				        "%i: system time seems to've gone backwards, "
				        "wall-time limit unenforceable!\n",
				        pro->rank);
			} else if ((time(NULL) - start_time) > (60 * ru.tim_lim)) {
				/* The wall-time quota's been used up, stop running. */
				break;
			}
		}

		/* Master process (i.e. the 1st process) broadcasts its flag that
		   tells the other processes whether to halt or not. If the halting
		   flag is set, each process stops simulating once the flag is
		   sent & received. */
		if (MPI_Bcast(&master_proc_orders_a_halt, 1, MPI_UNSIGNED_CHAR, 0, pro->work_comm) != MPI_SUCCESS) {
			e_puts("0: failed to tell other processes whether to continue simulating.");
		}
		if (master_proc_orders_a_halt) {
			/* The master process called a halt to the simulation. */
			break;
		}

		/* Processes that are no longer needed have nothing more to do
		   on this time step; they can just return to the loop's start. */
		if (((unsigned int) pro->rank) >= ru.N) {
			continue;
		}

		/* Master process broadcasts the particles' state to the other
		   workers. */
		if (MPI_Bcast(pro->parts, ru.N, MPIParticle, 0, pro->work_comm) != MPI_SUCCESS) {
			fprintf(stderr,
			        "%i: failed to send/receive particle state data.\n",
			        pro->rank);
			break;
		}

		#ifndef DUST_FIELD_ONLY
		/* Build the particle tree using the latest particle states. */
		if (build_tree(pro, ru) == NULL) {
			fprintf(stderr, "%i:%lu: tree building failed.\n", pro->rank,
			        pro->t_steps);
			break;
		}
		#endif

		#ifdef GL
		if ((!pro->rank) && !(pro->t_steps % 4)) {
			/* Refresh the particle display (and tracked particle
			   trajectories, where necessary) every 4 time steps. */
			for (idx = 0; idx < ru.track; idx++) {
				/* Set convenience pointers to the 3 components of this
				   particle's last recorded trajectory point. */
				l_t[0] = &(pro->traj_x[idx][TRAJ_LEN - 1]);
				l_t[1] = &(pro->traj_y[idx][TRAJ_LEN - 1]);
				l_t[2] = &(pro->traj_z[idx][TRAJ_LEN - 1]);
				p_idx = (idx * ((ru.N - 1) / ru.track)) + 1;
				if (pythag(*(l_t[0]) - (pro->parts[p_idx].x / ru.mpp),
				           *(l_t[1]) - (pro->parts[p_idx].y / ru.mpp),
				           *(l_t[2]) - (pro->parts[p_idx].z / ru.mpp)) > 0.4) {
					/* The particle's moved by more than 0.4 pixels from
					   the last recorded point in its trajectory, so it's
					   worth updating its trajectory. (If the particle
					   hasn't moved that much, this block of code isn't
					   run, so the trajectory isn't updated; this saves
					   on memory by making better use of the `TRAJ_LEN`
					   points that can be stored.) */
					memmove(&(pro->traj_x[idx][0]), &(pro->traj_x[idx][1]),
					        sizeof(float) * (TRAJ_LEN - 1));
					memmove(&(pro->traj_y[idx][0]), &(pro->traj_y[idx][1]),
					        sizeof(float) * (TRAJ_LEN - 1));
					memmove(&(pro->traj_z[idx][0]), &(pro->traj_z[idx][1]),
					        sizeof(float) * (TRAJ_LEN - 1));
					*(l_t[0]) = pro->parts[p_idx].x / ru.mpp;
					*(l_t[1]) = pro->parts[p_idx].y / ru.mpp;
					*(l_t[2]) = pro->parts[p_idx].z / ru.mpp;
				}
			}
			blank_window();
			draw_particles(qr, pro, ru);
			glfwSwapBuffers();
		}
		#endif

		/* If the user requested that the full simulation state be
		   recorded, write the full state to a file periodically. */
		if ((ru.sw_path != NULL) && !(pro->t_steps % 50000)) {
			write_sim_state_to_file(ru.N, pro, ru.sw_path);
		}

		#ifndef RECORD_MACRO_DATA
		/* Decide whether to record macroscopic information about the
		   simulated system on this time step. */
		record_macro_this_step = false;
		#else
		record_macro_this_step = (pro->t_steps <= 200);
		if ((pro->t_steps * ru.dt) > (2.0 * ru.settle)) {
			record_macro_this_step |= !(pro->t_steps % 2000);
		} else {
			record_macro_this_step |= !(pro->t_steps % 200);
		}
		#ifdef DEBUG_LANGMUIR
		/* It's useful to record macroscopic summary data more often when
		   running a Langmuir oscillation test. */
		record_macro_this_step |= !(pro->t_steps % 20);
		#endif
		#endif

		/* Step this process' particles through one time step using
		   the tree of particles just built.
		   Also, if macroscopic summary data are going to be recorded
		   below for this time step, tell `move_particles` to recalculate
		   the total electrostatic potential energy in the system. */
		#ifdef DUST_FIELD_ONLY
		move_particles(pro, ru, MPIParticle,
		               record_macro_this_step, pro->t_steps);
		#else
		move_particles(pro, ru, pro->root, MPIParticle,
		               record_macro_this_step, pro->t_steps);
		#endif

		#ifdef RECORD_PHI_DATA
		if ((!pro->rank) && !(pro->t_steps % 20000)) {
			/* Record the current electric potential profile. */
			write_phi_data_to_file(RECORD_PHI_DATA, pro, ru,
			                       pro->t_steps, pro->root);
		}
		#endif

		#ifdef RECORD_V_R_DATA
		if ((!pro->rank) && !(pro->t_steps % 1000)) {
			write_v_r_data_to_file(RECORD_V_R_DATA, pro->parts, ru,
			                       pro->t_steps);
		}
		#endif

		#ifdef RECORD_V_EACH_DATA
		if (!pro->rank) {
			if ((!(pro->t_steps % 20000)) || ((pro->t_steps <= 50000) && !(pro->t_steps % 5000))) {
				write_v_each_data_to_file(RECORD_V_EACH_DATA, pro->parts,
				                          ru, pro->t_steps);
			}
		}
		#endif

		#ifdef RECORD_N_DATA
		if ((!pro->rank) && !(pro->t_steps % 20000)) {
			write_n_data_to_file(RECORD_N_DATA, ru, pro, pro->t_steps);
		}
		#endif

		#ifdef RECORD_MACRO_DATA
		if ((!pro->rank) && record_macro_this_step) {
			/* Record whole-system information about the simulation. */
			write_macroscopic_data_to_file(ru, pro, pro->t_steps);
		}
		#endif

		/* If the GUI's enabled, check for any keypresses. If the user
		   presses Esc, the 1st process detects this and sets a flag
		   telling the other processes to stop. */
		#ifdef GL
		if (!pro->rank) {
			if (glfwGetKey(GLFW_KEY_ESC) || !glfwGetWindowParam(GLFW_OPENED)) {
				master_proc_orders_a_halt = 1;
			}
			if (glfwGetKey('W') && (cam[0] > 3.0)) {
				cam[0] -= 3.0;
			}
			#ifdef SPHERICAL
			if (glfwGetKey('S') && (cam[0] < (3.0 * ru.r))) {
			#else
			if (glfwGetKey('S') && (cam[0] < (1.5 * ru.L))) {
			#endif
				cam[0] += 3.0;
			}
			if (glfwGetKey('A')) {
				cam[1] += PI / 400.0;
			}
			if (glfwGetKey('D')) {
				cam[1] -= PI / 400.0;
			}
			if (glfwGetKey('Q')) {
				cam[2] += PI / 400.0;
			}
			if (glfwGetKey('E')) {
				cam[2] -= PI / 400.0;
			}
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			/* Orient the camera so that it's placed according to the user's
			   preferred spherical coordinates, looking at the dust grain
			   with its upward face pointing in the theta direction.
			   (This "up" vector's magnitude is unimportant.) */
			up[0] = cam[0] * cos(cam[1]) * cos(cam[2]);
			up[1] = cam[0] * sin(cam[1]) * cos(cam[2]);
			up[2] = cam[0] * -sin(cam[2]);
			#ifdef FPE
			disable_fpes();
			#endif
			#ifdef SPHERICAL
			gluLookAt(ru.r + cam[0] * cos(cam[1]) * sin(cam[2]),
			          ru.r + cam[0] * sin(cam[1]) * sin(cam[2]),
			          ru.r + cam[0] * cos(cam[2]),
			          ru.r, ru.r, ru.r,
			          up[0], up[1], up[2]);
			#else
			gluLookAt((ru.L / 2.0) + cam[0] * cos(cam[1]) * sin(cam[2]),
			          (ru.L / 2.0) + cam[0] * sin(cam[1]) * sin(cam[2]),
			          (ru.L / 2.0) + cam[0] * cos(cam[2]),
			          ru.L / 2.0, ru.L / 2.0, ru.L / 2.0,
			          up[0], up[1], up[2]);
			#endif
			#ifdef FPE
			enable_fpes();
			#endif
		}
		#endif
	}

	#ifdef GL
	/* If the user requested particle trajectory tracking, free the memory
	   allocated above for that. */
	if ((((unsigned int) pro->rank) < ru.N) && ru.track) {
		free(pro->traj_reset);
		if (!pro->rank) {
			free(pro->traj_x[0]);
			free(pro->traj_x);
		}
	}
	#endif

	/* Free the memory for the MPI particle type set up earlier. */ 
	MPI_Type_free(&MPIParticle);

	/* Done running; if this process was working on stepping particles
	   it now frees the memory for its list of particle array indices,
	   particles' mid-time step accelerations, particle tree, and particles'
	   positions on the previous time step. */
    if (((unsigned int) pro->rank) < ru.N) {
		free(pro->root);
		free(pro->a_x1);
		free(pro->pits[0]);
		free(pro->pits);
		free(pro->prev_x);
	}

	/* Free the memory the 1st process allocated to store the start of the
	   part of each process' particle array being transmitted during a
	   given time step. */
	if (!pro->rank) {
		free(pro->starts);
	}

	/* Record the final microscopic simulation state, if requested. */
	if (ru.sw_path != NULL) {
		write_sim_state_to_file(ru.N, pro, ru.sw_path);
	}

	#ifdef RECORD_N_DATA
	if (!pro->rank) {
		/* Record the final spatial particle distribution. */
		/* FIXME: need to check pot didn't already write this information
		   above when the time step was a round number. Otherwise, when the
		   final results get aggregated by another program, double
		   counting may occur! */
		write_n_data_to_file(RECORD_N_DATA, ru, pro, pro->t_steps - 1);
	}
	#endif

	#ifdef RECORD_MACRO_DATA
	if ((!pro->rank) && !record_macro_this_step) {
		/* Record final whole-system information about the simulation. */
		/* FIXME: see FIXME just above re: spatial particle distribution. */
		/* FIXME: if simulation ends unexpectedly, the electrostatic pot'l
		   energy totals in `pro.epe` will be outdated. */
		write_macroscopic_data_to_file(ru, pro, pro->t_steps - 1);
	}
	#endif
}

unsigned int set_run_params_from_args(unsigned int argc, const char* argv[], Run* ru)
{
	unsigned int n;

	/* Set reasonable defaults (for my computer, anyway). */
	ru->N = DEFAULT_N;
	#ifdef SPHERICAL
	ru->r = DEFAULT_R;
	#else
	ru->L = DEFAULT_L;
	#endif
	ru->dt = DEFAULT_DT;
	ru->mpp = DEFAULT_MPP;
	ru->run_for = DEFAULT_RUN_FOR;  /* if zero, run indefinitely */
	ru->soft = DEFAULT_SOFT;
	ru->Te = DEFAULT_T_E;
	ru->Ti = DEFAULT_T_I;
	ru->x_drift = DEFAULT_X_DRIFT;
	ru->a = DEFAULT_A;
	#ifdef GL
	ru->track = DEFAULT_TRACK;
	#endif
	ru->settle = DEFAULT_SETTLE;
	ru->sw_path = NULL;
	ru->sr_path = NULL;
	ru->tim_lim = DEFAULT_TIM_LIM;  /* if zero, no runtime limit */
	#if defined(PROLATE) || defined(OBLATE)
	ru->A = DEFAULT_ASPECT_RATIO;
	#endif

	for (n = 1; n < argc; n++) {
		if (argv[n][0] == '-') {
			switch (argv[n][1]) {

			/* FIXME: may be best to check `errno` properly after
			   these `strtoul` calls at some point. */

			case 'a':
				if (++n >= argc) {
					return 1;
				}
				ru->a = strtod(argv[n], NULL);
				if (ru->a <= 0.0) {
					e_puts("Grain radius must be positive.");
					return 2;
				} else if (ru->a == HUGE_VAL) {
					e_puts("Given grain radius is too huge!");
					return 2;
				}
			break;

			#ifdef PROLATE
			/* I'm not sure what `HUGE_VAL` or `HUGE_VALL` are, so will
			   omit any comparison to these elusive numbers - JH */
			case 'A':
				if (++n >= argc) {
					return 1;
				}
				ru->A = strtod(argv[n], NULL);
				if (ru->A <= 1.0) {
					e_puts("Prolate aspect ratio must be greater than 1.");
					return 2;
				}
			break;
			#elif defined(OBLATE)
			case 'A':
				if (++n >= argc) {
					return 1;
				}
				ru->A = strtod(argv[n], NULL);
				if (ru->A >= 1.0) {
					e_puts("Oblate aspect ratio must be less than 1.");
					return 2;
				}
			break;
			#else
			case 'A':
				e_puts("Aspect ratio is not a parameter for a sphere.");
				return 2;
			break;
			#endif

			case 'd':
				if (++n >= argc) {
					return 1;
				}
				ru->dt = strtold(argv[n], NULL);
				if (ru->dt <= 0.0) {
					e_puts("Time step must be positive.");
					return 2;
				} else if (ru->dt == HUGE_VALL) {
					e_puts("Given time step is too huge!");
					return 2;
				}
			break;

			case 'E':
				if (++n >= argc) {
					return 1;
				}
				ru->Te = strtold(argv[n], NULL);
				if (ru->Te < 0.0) {
					e_puts("Electron temperature must be non-negative.");
					return 2;
				} else if (ru->Te == HUGE_VALL) {
					e_puts("Given electron temperature is too huge!");
					return 2;
				}
			break;

			case 'e':
				if (++n >= argc) {
					return 1;
				}
				ru->settle = strtold(argv[n], NULL);
				if (ru->settle < 0.0) {
					e_puts("Settling time must be non-negative.");
					return 2;
				} else if (ru->settle == HUGE_VALL) {
					e_puts("Given settling time is too huge!");
					return 2;
				}
			break;

			case 'f':
				if (++n >= argc) {
					return 1;
				}
				ru->sw_path = malloc(strlen(argv[n]) + 1);
				if (ru->sw_path == NULL) {
					e_puts("Can't allocate memory for -f's argument.");
					return 2;
				}
				strcpy(ru->sw_path, argv[n]);
			break;

			case 'g':
				if (++n >= argc) {
					return 1;
				}
				/* Read a path from which to read the initial simulation
				   state. Notice that the file path isn't copied into a
				   separate string; `sr_path` just points directly into
				   `argv`. This should be OK since this path's only going
				   to be used once or twice by the `main` function (and
				   possibly `output_args_for_run`) shortly after this
				   function returns. */
				ru->sr_path = (char*) argv[n];
			break;

			case 'I':
				if (++n >= argc) {
					return 1;
				}
				ru->Ti = strtold(argv[n], NULL);
				if (ru->Ti < 0.0) {
					e_puts("Ion temperature must be non-negative.");
					return 2;
				} else if (ru->Ti == HUGE_VALL) {
					e_puts("Given ion temperature is too huge!");
					return 2;
				}
			break;

			case 'i':
				if (++n >= argc) {
					return 1;
				}
				if ((ru->run_for = strtoul(argv[n], NULL, 10)) == 0) {
					e_puts("Number of time steps to run for must be a valid nonzero integer.");
					return 2;
				}
			break;

			#ifndef SPHERICAL
			case 'L':
				if (++n >= argc) {
					return 1;
				}
				if ((ru->L = strtoul(argv[n], NULL, 10)) == 0) {
					e_puts("Simulation length in pixels must be a valid nonzero integer.");
					return 2;
				}
			break;
			#endif

			case 'l':
			return 3;  /* tell main() to list default `ru` values */

			case 'm':
				if (++n >= argc) {
					return 1;
				}
				ru->mpp = strtod(argv[n], NULL);
				if (ru->mpp <= 0.0) {
					e_puts("Pixel size must be positive.");
					return 2;
				} else if (ru->mpp == HUGE_VAL) {
					e_puts("Given pixel size is too huge!");
					return 2;
				}
			break;

			case 'N':
				if (++n >= argc) {
					return 1;
				}
				if ((ru->N = strtoul(argv[n], NULL, 10)) == 0) {
					e_puts("Number of particles must be a valid nonzero integer.");
					return 2;
				}
				#ifdef DEBUG_FIELDS
				if (ru->N <= DEBUG_FIELDS) {
					fprintf(stderr, "Particle number must be at least %u.\n",
					        DEBUG_FIELDS + 1);
					return 2;
				}
				#endif
			break;

			#ifdef SPHERICAL
			case 'r':
				if (++n >= argc) {
					return 1;
				}
				if ((ru->r = strtoul(argv[n], NULL, 10)) == 0) {
					e_puts("Simulation radius in pixels must be a valid nonzero integer.");
					return 2;
				}
			break;
			#endif

			case 's':
				if (++n >= argc) {
					return 1;
				}
				ru->soft = strtod(argv[n], NULL);
				if (ru->soft < 0.0) {
					e_puts("Softening parameter must be non-negative.");
					return 2;
				} else if (ru->soft == HUGE_VAL) {
					e_puts("Given softening parameter is too huge!");
					return 2;
				}
			break;

			#ifdef GL
			case 't':
				if (++n >= argc) {
					return 1;
				}
				if (sscanf(argv[n], "%hu", &(ru->track)) != 1) {
					e_puts("Can't read number of trajectories to track.");
					return 2;
				}
			break;
			#endif

			case 'w':
				if (++n >= argc) {
					return 1;
				}
				ru->tim_lim = strtod(argv[n], NULL);
				if (ru->tim_lim < 0.0) {
					e_puts("Wall-time limit must be non-negative.");
					return 2;
				} else if (ru->tim_lim == HUGE_VAL) {
					e_puts("Wall-time limit is too huge!");
					return 2;
				}
			break;

			case 'x':
				if (++n >= argc) {
					return 1;
				}
				ru->x_drift = strtod(argv[n], NULL);
				if (ru->x_drift == HUGE_VAL) {
					e_puts("Given drift velocity is too huge!");
					return 2;
				}
			break;

			default:
			return 1;

			}
		} else {
			/* Argument given isn't a flag or a flag's argument,
			   that's not right! */
			return 1;
		}
	}

	#ifdef GL
	if (ru->track >= ru->N) {
		fprintf(stderr, "%hu particle tracks requested for %lu particles, "
		        "recording only %lu tracks.\n", ru->track, ru->N, ru->N - 1);
		ru->track = ru->N - 1;
	}
	#endif

	return 0;
}

void usage(const char* prog_name)
{
	fprintf(stderr, "Usage: %s [OPTIONS]\n", prog_name);
	e_puts("  -a [SIZE]     grain radius in metres");
	e_puts("  -d [NUMBER]   time step in seconds");
	e_puts("  -E [NUMBER]   initial electron temperature in Kelvin");
	e_puts("  -e [NUMBER]   initial settling/equilibration time");
	e_puts("  -f [PATH]     write final simulation state to [PATH]");
	e_puts("  -g [PATH]     read initial simulation state from [PATH]");
	e_puts("  -I [NUMBER]   initial ion temperature in Kelvin");
	e_puts("  -i [INTEGER]  number of time steps to run for");
	#ifndef SPHERICAL
	e_puts("  -L [INTEGER]  simulation length in pixels");
	#endif
	e_puts("  -l            show runtime settings as command-line arguments");
	e_puts("  -m [SIZE]     pixel length in metres");
	e_puts("  -N [INTEGER]  number of particles to simulate");
	#ifdef SPHERICAL
	e_puts("  -r [INTEGER]  simulation radius in pixels");
	#endif
	e_puts("  -s [NUMBER]   softening parameter in metres");
	#ifdef GL
	e_puts("  -t [NUMBER]   number of particle trajectories to display");
	#endif
	e_puts("  -w [NUMBER]   limit on wall-time spent running in minutes");
	e_puts("  -x [NUMBER]   drift velocity in x-direction");
	#ifdef PROLATE
	e_puts("  -A [NUMBER]   aspect ratio of prolate spheroid");
	#endif
	#ifdef OBLATE
	e_puts("  -A [NUMBER]   aspect ratio of oblate spheroid");
	#endif
}

void output_args_for_run(const char* prog_name, Run ru)
{
	printf("Settings: %s -N %lu", prog_name, ru.N);

	#ifdef SPHERICAL
	printf(" -r %u", ru.r);
	#else
	printf(" -L %u", ru.L);
	#endif

	printf(" -d %Lg -m %g -s %g -E %g -I %g -x %g -a %g -e %g",
	       ru.dt, ru.mpp, ru.soft, ru.Te, ru.Ti, ru.x_drift, ru.a, ru.settle);

	#if defined(PROLATE) || defined(OBLATE)
	printf(" -A %g", ru.A);
	#endif

	if (ru.run_for) {
		printf(" -i %lu", ru.run_for);
	}

	#ifdef GL
	printf(" -t %hu", ru.track);
	#endif

	if (ru.sw_path) {
		printf(" -f %s", ru.sw_path);
	}
	if (ru.sr_path) {
		printf(" -g %s", ru.sr_path);
	}

	if (ru.tim_lim) {
		printf(" -w %g", ru.tim_lim);
	}

	putchar('\n');
}

#ifdef DEBUG_LANGMUIR
/* Initiate a Langmuir wave (electron plasma oscillation) by shifting some
   of the electrons in the x-direction. */
void set_up_langmuir_wave(Proc* pro, Run ru)
{
	unsigned int p_idx;

	#ifdef SPHERICAL
	#error "for Langmuir wave tests, use a cubical domain"
	#else
	for (p_idx = 1; p_idx < ru.N; p_idx++) {
		if (pro->parts[p_idx].species != ELECTRON) {
			continue;
		} else if (pro->parts[p_idx].x < (0.425 * ru.mpp * ru.L)) {
			continue;
		} else if (pro->parts[p_idx].x > (0.575 * ru.mpp * ru.L)) {
			continue;
		}
		pro->parts[p_idx].x += 0.15 * ru.mpp * ru.L;
	}
	#endif
}
#endif

/* As part of cleanly terminating a process, free now-unneeded memory from
   `ru` & `pro`, and close the MPI environment. */
void clean_up(Run* ru, Proc* pro)
{
	int init_flag;

	#ifdef RECORD_MACRO_DATA
	if ((!pro->rank) && (ru->N >= 6) && (pro->fp_macro != NULL)) {
		/* The master process closes the macroscopic summary data file.
		   (Not much point checking whether it succeeds, frankly). */
		fclose(pro->fp_macro);
		pro->fp_macro = NULL;
	}
	#endif

	if (pro->new_parts != NULL) {
		free(pro->new_parts);
		pro->new_parts = NULL;
	}
	if (pro->parts != NULL) {
		free(pro->parts);
		pro->parts = NULL;
	}

	if (ru->sw_path != NULL) {
		free(ru->sw_path);
		ru->sw_path = NULL;
	}

	MPI_Initialized(&init_flag);
	if (init_flag) {
		MPI_Finalize();
	}
}

int main(int argc, char* argv[])
{
	Proc pro;
	#ifdef GL
	GLUquadric* qr = NULL;
	#endif
	Run ru;
	FILE* urandom_fp;
	uint32_t urandom_seed[16];
	#ifdef GL
	char window_title[75];
	#endif

	/* Initialize MPI and fetch this process' MPI ID. */
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		e_puts("Failed to initialize MPI environment.");
		return EXIT_FAILURE;
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &(pro.rank));

	/* Fetch the number of processes running, store it, and check that
	   this process' rank is less than the number of running processes. */
	MPI_Comm_size(MPI_COMM_WORLD, &(pro.num_procs));
	if (pro.num_procs < (pro.rank + 1)) {
		fprintf(stderr, "%i: MPI environment reports %i processes running!\n",
		        pro.rank, pro.num_procs);
		abort();
	}

	/* Before parsing the program's arguments, see whether there are any. */
	if ((argc < 2) && !pro.rank) {
		/* No arguments given. The user can do this, but let them know
		   in case this was accidental. */
		e_puts("Run without arguments; attempting to use default settings.");
	}

	/* Set the particle array data pointers to NULL before `clean_up`
	   might be called, so that if `clean_up` is called it refrains from
	   trying to free these pointers. */
	pro.parts = NULL;
	pro.new_parts = NULL;
	pro.fp_macro = NULL;

	/* Try to configure program settings based on the program's arguments.
	   (As this is quick, it makes sense for each process to do this.) */
	/* FIXME: if something goes wrong during the function call, an error
	   message gets displayed once for every single process! */
	switch (set_run_params_from_args(argc, (const char**) argv, &ru)) {

	case 1:  /* remind user of usage syntax, then exit */
		if (!pro.rank) {
			usage(argv[0]);
		}
		clean_up(&ru, &pro);
	return EXIT_FAILURE;

	case 2:  /* specific error message already shown, just exit */
		clean_up(&ru, &pro);
	return EXIT_FAILURE;

	case 3:  /* user would like to see default settings */
		if (!pro.rank) {
			output_args_for_run(argv[0], ru);
		}
		clean_up(&ru, &pro);
	return EXIT_SUCCESS;

	default:  /* everything OK, carry on */
	break;

	}

	/* Check that the dust grain isn't bigger than the simulation itself! */
	#ifdef SPHERICAL
		#ifdef PROLATE
		if (ru.a * pow(ru.A,2.0/3.0) >= (ru.r * ru.mpp)) {
			fprintf(stderr,
			        "Grain half-length (%g m) must be less than "
			        "simulation region radius (%g m).\n",
			        ru.a * pow(ru.A,2.0/3.0), ru.r * ru.mpp);
		#elif defined(OBLATE)
		if (ru.a * pow(ru.A,-1.0/3.0) >= (ru.r * ru.mpp)) {
			fprintf(stderr,
			        "Grain radius (%g m) must be less than "
			        "simulation region radius (%g m).\n",
			        ru.a * pow(ru.A,-1.0/3.0), ru.r * ru.mpp);
		#else
		if (ru.a >= (ru.r * ru.mpp)) {
			fprintf(stderr,
			        "Grain radius (%g m) must be less than "
			        "simulation region radius (%g m).\n",
			        ru.a, ru.r * ru.mpp);
		#endif
	#else
		#ifdef PROLATE
		if (ru.a * pow(ru.A,2.0/3.0) >= (ru.L * ru.mpp / 2.0)) {
			fprintf(stderr,
			        "Grain half-length (%g m) must be less than "
			        "simulation region half-size (%g m).\n",
			        ru.a * pow(ru.A,2.0/3.0), ru.L * ru.mpp / 2.0);
		#elif defined(OBLATE)
		if (ru.a * pow(ru.A,-1.0/3.0) >= (ru.L * ru.mpp / 2.0)) {
			fprintf(stderr,
			        "Grain half-length (%g m) must be less than "
			        "simulation region half-size (%g m).\n",
			        ru.a * pow(ru.A,-1.0/3.0), ru.L * ru.mpp / 2.0);
		#else
		if (ru.a >= (ru.L * ru.mpp / 2.0)) {
			fprintf(stderr,
			        "Grain radius (%g m) must be less than "
			        "simulation region half-size (%g m).\n",
			        ru.a, ru.L * ru.mpp / 2.0);
		#endif
	#endif
		clean_up(&ru, &pro);
		return EXIT_FAILURE;
	}

	#ifdef FPE
	/* Disable masking of most floating point exceptions.
	   This is not ANSI C! */
	enable_fpes();
	#endif

	/* Work out which subset of particles in the array each process is
	   responsible for stepping, and set up the MPI communicator for
	   inter-worker communication. */
	if (set_up_workers_for_assignments(&pro, ru)) {
		abort();
	}

	/* If the user's running this program over more processes than
	   there are particles, the surplus processes excuse themselves 
	   and tell the user that they're not needed.
	   Notice that this is only done once `pro.work_comm` is set by
	   the `set_up_workers_for_assignments` call above; setting up
	   the `pro.work_comm` communicator wouldn't work otherwise, and
	   the program would freeze. */
	if (((unsigned int) pro.rank) >= ru.N) {
		fprintf(stderr, "%i: process unneeded as there are only "
		        "%lu particles, ending.\n", pro.rank, ru.N);
		clean_up(&ru, &pro);
		return EXIT_SUCCESS;
	}

	/* Allocate memory to store particles' states in `pro.parts`, as well
	   as temporary storage for time-stepped particles in `pro.new_parts`. */
	if ((pro.parts = malloc(sizeof(Particle) * ru.N)) == NULL) {
		fprintf(stderr,
		        "%i: can't allocate %.1f MiB for %lu particles' data.\n",
		        pro.rank, sizeof(Particle) * ru.N / 1048576.0, ru.N);
		clean_up(&ru, &pro);
		return EXIT_FAILURE;
	}
	if ((pro.new_parts = malloc(sizeof(Particle) * ru.N)) == NULL) {
		fprintf(stderr, "%i: can't allocate %.1f MiB for "
		        "%lu updated particle states.\n",
		        pro.rank, sizeof(Particle) * ru.N / 1048576.0, ru.N);
		clean_up(&ru, &pro);
		return EXIT_FAILURE;
	}

	/* If the user wants this program to take its initial simulation state
	   from a file, read that file and load that state. If not, generate a
	   seed for the PRNG, initialize the PRNG, and set the number of elapsed
	   time steps to zero. (A little later the master process will also need
	   to generate a random initial state for the particles, but not yet.) */
	if (ru.sr_path != NULL) {
		/* Try to load an initial simulation state (particles and PRNG)
		   from a simulation state file. */
		if (read_sim_state_from_file(ru.N, &pro, ru.sr_path)) {
			clean_up(&ru, &pro);
			return EXIT_FAILURE;
		}
		/* Has the user requested the simulation run for a particular
		   number of time steps? And if so, has the simulation just loaded
		   already run for that long? If both conditions are true, there's
		   no point running the simulation any further, so stop. */
		if (ru.run_for && (pro.t_steps >= ru.run_for)) {
			printf("%i: simulation's already run for %lu steps, exiting.\n",
			       pro.rank, ru.run_for);
			clean_up(&ru, &pro);
			return EXIT_SUCCESS;
		}
	} else {
		/* Try to get a high-quality PRNG seed from the file/device at the
		   path `RANDOMNESS_SOURCE` (/dev/urandom by default). If that's
		   not possible, fall back on using the MPI process ID as a (totally
		   inadequate) seed. */
		if (((urandom_fp = fopen(RANDOMNESS_SOURCE, "rb")) != NULL)
		    && (fread(urandom_seed, sizeof(uint32_t),
			          16, urandom_fp) == 16)) {
			init_prng(&(pro.ps), urandom_seed);
			fclose(urandom_fp);  /* little point checking the return value */
		} else {
			fprintf(stderr, "%i: can't read PRNG seed from %s,"
			                " PRNG may behave inadequately.\n",
			        pro.rank, RANDOMNESS_SOURCE);
			urandom_seed[0] = pro.rank;
			init_prng(&(pro.ps), urandom_seed);
		}
		/* This is the beginning of a new simulation, so set the number of
		   elapsed time steps and absorbed electrons to zero, initialize
		   the grain's charge & mass, and (if necessary) set an initial
		   estimate of the electric potential at the simulation edge. */
		pro.t_steps = 0;
		pro.tot_absor_e = 0;
		pro.dust_grain_q = CHARGE[DUST];
		pro.dust_grain_m = MASS[DUST];
		#ifdef SCEPTIC_REINJECTION
		pro.chi_b = 0.0;
		#endif
	}

	#ifdef DEBUG_RUTHERFORD
	/* Write column headings for the Rutherford scattering test results. */
	if (!pro.rank) {
		puts("B VX1 VY1 VZ1 VX2 VY2 VZ2");
	}
	#endif

	if (!pro.rank) {

		if (ru.sr_path == NULL) {
			/* Initialize the particles with random initial states. */
			init_particles(&pro, ru);
		}

		#ifdef DEBUG_LANGMUIR
		/* Set up a test Langmuir wave. */
		set_up_langmuir_wave(&pro, ru);
		#endif

		#ifdef RECORD_MACRO_DATA
		/* If enough particles are being simulated to warrant doing so,
		   open the macroscopic summary data output file and, if that file's
		   newly created, write its header. */
		if (ru.N >= 6) {
			if ((pro.fp_macro = fopen(RECORD_MACRO_DATA, "ab")) == NULL) {
				fprintf(stderr,
				        "Can't open %s to write macroscopic data.\n",
				        RECORD_MACRO_DATA);
			}
			if ((pro.fp_macro != NULL) && !ftell(pro.fp_macro)) {
				fputs("TIME KE PE KI PI KG PG XE XESD DGQ DL DN "
				      "SE1 SE2 SE3 SE4 SI1 SI2 SI3 SI4 TE TI LMEX LMEY LMEZ "
				      "LMIX LMIY LMIZ AMEX AMEY AMEZ AMIX AMIY AMIZ TOTE\n",
				      pro.fp_macro);
			}
		}
		#endif

		/* If the GUI's enabled, initialize GLFW and open a window. */
		#ifdef GL
		snprintf(window_title, 75,
		         "Plasma octree simulation (N = %lu)", ru.N);
		/* FIXME: do these have return values that need checking? */
		if (glfwInit() != GL_TRUE) {
			e_puts("Failed to initialize GLFW.");
			clean_up(&ru, &pro);
			return EXIT_FAILURE;
		}
		#ifdef SPHERICAL
		glfwOpenWindow(ru.r, ru.r, 0, 0, 0, 0, 0, 0, GLFW_WINDOW);
		#else
		glfwOpenWindow(ru.L / 2.0, ru.L / 2.0, 0, 0, 0, 0, 0, 0, GLFW_WINDOW);
		#endif
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);  /* stops hole appearing in alpha'd spheres */
		glfwSetWindowTitle(window_title);
		if (!(qr = gluNewQuadric())) {
			e_puts("Failed to initialize quadric for displaying particles.");
			clean_up(&ru, &pro);
			return EXIT_FAILURE;
		}
		#ifdef SPHERICAL
		glViewport(0, 0, ru.r, ru.r);
		#else
		glViewport(0, 0, ru.L / 2.0, ru.L / 2.0);
		#endif
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
/*
		glOrtho(0, 1000, 0, 1000, -1, 1);
*/
		/* Apparently `gluPerspective` is deprecated in OpenGL 3.1, so at
		   some point it might be wise to replace it with, well, whatever
		   ought to replace it. But for the time being, use it to set up
		   the perspective projection matrix (and I should look up what
		   that is precisely). */
		#ifdef SPHERICAL
		gluPerspective(65.0, 1.0, 0.5, 6.0 * ru.r);
		#else
		gluPerspective(65.0, 1.0, 0.5, 3.0 * ru.L);
		#endif

		/* Prepare for `run` to point the OpenGL camera. */
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		#endif
	}

	/* Run the main loop. */
	#ifdef GL
	run(qr, &pro, ru);
	#else
	run(&pro, ru);
	#endif

	/* If the GUI's enabled, the master process shuts down GLFW and
	   deletes the now-unneeded quadric. */
	#ifdef GL
	if (!pro.rank) {
		gluDeleteQuadric(qr);
		glfwTerminate();
	}
	#endif

	/* Finish by freeing now-unneeded memory from `ru` & `pro`, and
	   ending the MPI environment. */
	clean_up(&ru, &pro);

	return EXIT_SUCCESS;
}
