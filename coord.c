#include "coord.h"

#ifdef SPHERICAL
/* Convert a velocity in spherical coordinates to a velocity in Cartesian
   coordinates, at a position given in spherical coordinates. */
void velocity_in_cartesian_coords(double* v, double phi, double theta, double* vx, double* vy, double* vz)
{
	*vx = (v[0] * cos(phi) * sin(theta))
	      - (v[1] * sin(phi)) + (v[2] * cos(phi) * cos(theta));
	*vy = (v[0] * sin(phi) * sin(theta))
	      + (v[1] * cos(phi)) + (v[2] * sin(phi) * cos(theta));
	*vz = (v[0] * cos(theta)) - (v[2] * sin(theta));
}
#endif

#ifdef SCEPTIC_REINJECTION
/* Rotate the position vector (`x`, `y`, `z`) about the direction vector `d`
   (with no z-component!) through an angle with sine `s` and cosine `c`.
   Cf. <http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/>,
   end of section 5.2. */
void rotate_vec_about_xy_dir(const double* d, double s, double c, double *x, double* y, double* z, double offset)
{
	const double omc = 1.0 - c;
	double temp[3];

	temp[0] = d[0] * ((d[0] * *x) + (d[1] * *y)) * omc;
	temp[0] += (*x * c) + (d[1] * *z * s);
	temp[1] = d[1] * ((d[0] * *x) + (d[1] * *y)) * omc;
	temp[1] += (*y * c) - (d[0] * *z * s);
	temp[2] = *z * c;
	temp[2] += ((d[0] * *y) - (d[1] * *x)) * s;
	*x = offset + temp[0];
	*y = offset + temp[1];
	*z = offset + temp[2];
}
#endif

#ifdef BORIS
/* Rotate the position vector (`x`, `y`, `z`) about an arbitrary
   direction vector `d` through an angle with sine `s` and cosine `c`.
   This is used only by the Boris particle-stepping algorithm, to
   compute Larmor orbits for magnetic field sub-steps.
   Cf. <http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/>,
   end of section 5.2. */
void rotate_vec_about_dir(const double* d, double s, double c, double *x, double* y, double* z, double offset)
{
	const double omc = 1.0 - c;
	double temp[3];

	temp[0] = d[0] * ((d[0] * *x) + (d[1] * *y) + (d[2] * *z)) * omc;
	temp[0] += (*x * c) + (((d[1] * *z) - (d[2] * *y)) * s);
	temp[1] = d[1] * ((d[0] * *x) + (d[1] * *y) + (d[2] * *z)) * omc;
	temp[1] += (*y * c) + (((d[2] * *x) - (d[0] * *z)) * s);
	temp[2] = d[2] * ((d[0] * *x) + (d[1] * *y) + (d[2] * *z)) * omc;
	temp[2] += (*z * c) + (((d[0] * *y) - (d[1] * *x)) * s);
	*x = offset + temp[0];
	*y = offset + temp[1];
	*z = offset + temp[2];
}
#endif
