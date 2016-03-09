#include <math.h>

#ifdef SPHERICAL
void velocity_in_cartesian_coords(double* v, double phi, double theta, double* vx, double* vy, double* vz);
#endif

#ifdef SCEPTIC_REINJECTION
void rotate_vec_about_xy_dir(const double* d, double s, double c, double *x, double* y, double* z, double offset);
#endif

#ifdef BORIS
void rotate_vec_about_dir(const double* d, double s, double c, double *x, double* y, double* z, double offset);
#endif
