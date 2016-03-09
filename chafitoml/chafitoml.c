/* chafitoml: fit an OML charging curve to pot output */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* for `memmove` */

#define EPSILON0 8.8541878176e-12
#define ME 9.1093829e-31
#define PI 3.14159265359
#define QE 1.60217657e-19
#define KB 1.380649e-23

/* `fscanf` specification of the columns/fields formerly used in pot's
   macroscopic output file. */
#define OLD_FIELDS \
	"%lg %*g %*g %*g %*g %*g %*s %*g %*g %i" \
	" %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %lg %lg" \
	" %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g"

/* `fscanf` specification of the columns/fields now used in pot's
   macroscopic output file. */
#define NEW_FIELDS (OLD_FIELDS " %*u")

typedef struct Data {
	unsigned int N;
	double* t;
	int* q;
	double* T_e;
	double* T_i;
	bool any_charging;
} Data;

void init_Data(Data* d)
{
	d->N = 0;
	d->t = NULL;
	d->q = NULL;
	d->any_charging = false;
}

void free_Data_memory(Data* d)
{
	d->N = 0;
	if (d->t != NULL) {
		free(d->t);
	}
	if (d->q != NULL) {
		free(d->q);
	}
}

void seek_to_second_line_of(FILE* f)
{
	rewind(f);
	while ((!feof(f)) & !ferror(f)) {
		if (getc(f) == '\n') {
			return;
		}
	}
}

void read_macroscopic_data_file(Data* d, const char* prog_name, const char* path, double start_time)
{
	char buf[1024];
	size_t bytes_read;
	FILE* f = fopen(path, "rb");
	char* correct_field_format = NULL;
	unsigned int i;
	unsigned int init_N = d->N;
	unsigned int lines = 0;

	if (f == NULL) {
		fprintf(stderr, "%s: can't open %s\n", prog_name, path);
		return;
	}

	/* Quickly scan through the file, counting newlines. As there should be
	   one data point per line, this allows pre-estimation of the maximum
	   number of data points to read (and hence memory to allocate). */
	while ((!feof(f)) && !ferror(f)) {
		bytes_read = fread(buf, 1, 1024, f);
		for (i = 0; i < bytes_read; i++) {
			if (buf[i] == '\n') {
				d->N++;
			}
		}
	}

	/* Check there're actual data lines to be read in this file.
	   (Also, Linux will blithely open any directories passed to the
	   program as argument by mistake, and this seems to catch that.) */
	if ((d->N - init_N) < 2) {
		fprintf(stderr, "%s: %s has no data\n", prog_name, path);
		fclose(f);
		d->N = init_N;
		return;
	}
	d->N--;  /* first line's just a header, remember */

	/* Allocate enough memory for the data. */
	if ((d->t = realloc(d->t, 3 * sizeof(double) * d->N)) == NULL) {
		fprintf(stderr, "%s: can't allocate %lu KiB of memory for data\n",
		        prog_name, 3 * sizeof(double) * d->N / 1024);
		fclose(f);
		d->N = init_N;
		return;
	}
	if ((d->q = realloc(d->q, sizeof(int) * d->N)) == NULL) {
		fprintf(stderr, "%s: can't allocate %lu KiB of memory for data\n",
		        prog_name, sizeof(int) * d->N / 1024);
		fclose(f);
		d->N = init_N;
		return;
	}

	/* Go back to the file's second line (where the data begin) and use the
	   first line of data to discern the file's format (the newer format has
	   an extra field) by counting the number of spaces in it (i.e. the
	   number of fields minus one). */
	seek_to_second_line_of(f);
	i = 0;
	while ((!feof(f)) && !ferror(f)) {
		buf[0] = getc(f);
		if (buf[0] == '\n') {
			break;
		}
		i += (buf[0] == ' ');
	}
	switch (i) {
	case 33: correct_field_format = OLD_FIELDS; break;
	case 34: correct_field_format = NEW_FIELDS; break;
	default:
		fprintf(stderr, "%s: %s doesn't look like a valid data file\n",
		        prog_name, path);
		fclose(f);
		d->N = init_N;
	return;
	}

	/* Return to the start of the data. */
	seek_to_second_line_of(f);

	/* Point the temperature data pointers to the right subsection of
	   the main data array. Allow for the fact that temperature data
	   from previous data files may have been read into the main data
	   array already; move previously read data to their new home. */
	d->T_e = &(d->t[d->N]);
	d->T_i = &(d->t[2 * d->N]);
	if (init_N) {
		memmove(d->T_i, &(d->t[2 * init_N]), sizeof(double) * init_N);
		memmove(d->T_e, &(d->t[init_N]), sizeof(double) * init_N);
	}

	/* Read in the data. */
	i = init_N;
	while (((lines + init_N) < d->N) && (!feof(f)) && !ferror(f)) {
		if (fscanf(f, correct_field_format,
		           &(d->t[i]), &(d->q[i]), &(d->T_e[i]), &(d->T_i[i])) != 4) {
			fprintf(stderr, "%s: %s:%u: invalid line of data\n",
			        prog_name, path, lines + 1);
			break;
		}
		lines++;
		if (d->t[i] < start_time) {
			/* This datum is too early to be included.
			   And it should be zero anyway. */
			if (d->q[i]) {
				fprintf(stderr, "%s: %s:%u: charge is %i, should be 0\n",
				        prog_name, path, lines + 1, d->q[i]);
			}
			continue;
		}
		if (d->q[i]) {
			d->any_charging = true;  /* some charging occurred */
		}
		i++;
	}
	d->N = i;  /* update with actual number of data points included */

	fclose(f);
}

double uint_pow(double x, unsigned int y)
{
	double result = 1.0;

	while (y--) {
		result *= x;
	}

	return result;
}

unsigned int factorial(unsigned int x)
{
	unsigned int result = 1;

	while (x > 1) {
		result *= x;
		x--;
	}

	return result;
}

double prob_of_delta_q(int delta, double e_rate, double i_rate)
{
	unsigned int e_min;
	unsigned int inc = 0;
	double p;  /* probability of one particular way of getting `delta` */
	double tot_p = 0;  /* total estimated probability */

	/* Precomputed terms in the product of the electron Poisson PMF and
	   the ion Poisson PMF. */
	double e_rate_pow_fac_e_min;
	double exp_e_rate_i_rate = exp(-e_rate - i_rate);
	double i_rate_pow_fac_i_min;

	if (delta < 0) {
		e_min = -delta;
	} else {
		e_min = 0;
	}

	if (e_min < 13) {
		e_rate_pow_fac_e_min = uint_pow(e_rate, e_min)
		                       / (double) factorial(e_min);
	} else {
		/* The result of the factorial would be too huge for `factorial` to
		   successfully return, so roll out the heavy artillery. */
		e_rate_pow_fac_e_min = uint_pow(e_rate, e_min)
		                       / tgamma(1.0 + e_min);
	}
	i_rate_pow_fac_i_min = uint_pow(i_rate, delta + e_min)
	                       / (double) factorial(delta + e_min);
	do {
		p = exp_e_rate_i_rate * e_rate_pow_fac_e_min
		    * i_rate_pow_fac_i_min;
		tot_p += p;
		inc++;

		/* Update the combined power-factorial terms
		   for the next go-around. */
		e_rate_pow_fac_e_min *= e_rate / (e_min + inc);
		i_rate_pow_fac_i_min *= i_rate / (delta + e_min + inc);
	} while (p > 1e-10);

	return tot_p;
}

double log_like(const Data* d, double a, double n, double mime, int Z)
{
	double chi_3 = 4 * PI * a * a * n * sqrt(KB / (2 * PI * ME));
	double chi_4 = 4 * PI * EPSILON0 * KB * a / (QE * QE);
	unsigned int i;
	double ll;
	double rate[2];  /* expec'd rates of electron (0) & ion (1) collec. */
	double time_step;

	if (d->N < 2) {
		return 0.0;  /* Not enough data! */
	}

	/* Determine the log likelihood of the first datum given, assuming the
	   charge was zero at all preceding times. Note the presumption that
	   the time step equals the time between the first datum and the second
	   (since there's no way of taking the time between the imaginary
	   zeroth datum and the first). */
	rate[0] = chi_3 * sqrt(d->T_e[0]) * (d->t[1] - d->t[0]);
	rate[1] = chi_3 * sqrt(d->T_i[0] / mime) * (d->t[1] - d->t[0]);
	ll = log(prob_of_delta_q(d->q[0], rate[0], rate[1]));

	for (i = 1; i < d->N; i++) {
		time_step = d->t[i] - d->t[i-1];
		#ifdef DEBUG
		if (time_step > 1e-6) {
			fprintf(stderr, "time step %u is %g, worryingly long\n",
			        i, time_step);
		}
		#endif
		rate[0] = chi_3 * sqrt(d->T_e[i-1]);
		rate[0] *= exp(d->q[i-1] / (chi_4 * d->T_e[i-1]));
		rate[0] *= time_step;
		rate[1] = chi_3 * sqrt(d->T_i[i-1] / mime);
		rate[1] *= 1.0 - (Z * d->q[i-1] / (chi_4 * d->T_i[i-1]));
		rate[1] *= time_step;
		ll += log(prob_of_delta_q(d->q[i] - d->q[i-1], rate[0], rate[1]));
	}

	return ll;
}

void adjust_params_using_gradients(const double* grad, double* param, double* rss, char* points)
{
	unsigned int i;
	#ifdef DEBUG
	const char* param_name[3] = { "a", "n_e", "m_i / m_e" };
	#endif

	for (i = 1; i < 3; i++) {
		if ((grad[i] < 0.0) && (points[i] >= 0)) {
			points[i]++;
			if (points[i] == 5) {
				rss[i] *= 2.0;
				points[i] = 0;
				#ifdef DEBUG
				printf("doubling `rss[%u]`\n", i);
				#endif
			}
		} else if ((grad[i] > 0.0) && (points[i] <= 0)) {
			points[i]--;
			if (points[i] == -5) {
				rss[i] /= 2.0;
				points[i] = 0;
				#ifdef DEBUG
				printf("halving `rss[%u]`\n", i);
				#endif
			}
		} else {
			points[i] = 0;
		}
	}
	for (i = 0; i < 3; i++) {
		param[i] -= rss[i] * grad[i];
		while (param[i] < 0.0) {
			#ifdef DEBUG
			printf("quartering `rss[%u]`, it was so big %s < 0\n",
			       i, param_name[i]);
			#endif
			param[i] += rss[i] * grad[i];
			rss[i] /= 4.0;
			param[i] -= rss[i] * grad[i];
		}
	}
}

unsigned int fit_by_grad_desc(const char* prog_name, const Data* d, double* param, double* initial_nll, double* final_nll)
{
	double dee[3] = { param[0] / 1e4, param[1] / 1e4, param[2] / 1e4 };
	double delta_nll;
	double grad[3];
	double nll[5];
	char points[3] = { 0, 0, 0 };
/*
	double rss[3] = { param[0] / 1e9, param[1] * 1e12, param[2] * 100 };
*/
	double rss[3] = { param[0] / 1e9, param[1] * 1e12, param[2] * 50 };
	unsigned int iterations = 0;
	#ifdef DEBUG
	FILE* debug_out;
	#endif

	if (param[0] < 0.0) {
		fprintf(stderr, "%s: initial a %g < 0\n", prog_name, param[0]);
		return 0;
	}
	if (param[1] < 0.0) {
		fprintf(stderr, "%s: initial n_e %g < 0\n", prog_name, param[1]);
		return 0;
	}
	if (param[2] < 0.0) {
		fprintf(stderr, "%s initial m_i / m_e %g < 0\n", prog_name, param[2]);
		return 0;
	}

	nll[0] = -log_like(d, param[0], param[1], param[2], 1);
	*initial_nll = nll[0];

	if (!isfinite(nll[0])) {
		fprintf(stderr, "%s: initial parameter values give invalid NLL\n",
		        prog_name);
		fprintf(stderr, "%s: a = %g, n_e = %g, m_i / m_e = %g\n",
		        prog_name, param[0], param[1], param[2]);
		return 0;
	}

	#ifdef DEBUG
	printf("delta a:       %g\n", dee[0]);
	printf("delta n_e:     %g\n", dee[1]);
	printf("delta m_i/m_e: %g\n", dee[2]);
	puts("a           n_e         m_i/m_e old NLL "
	     "delta NLL  d(NLL)/da  d(NLL)/dn    d(NLL)/d(MIME) "
	     "RSS a    RSS n_e  RSS MIME");
	printf("%#.6g %#.6g %#.6g %#.6g\n",
	       param[0], param[1], param[2], nll[0]);
	if ((debug_out = fopen("debug.dat", "wb")) == NULL) {
		fprintf(stderr,
		        "%s: can't open debug.dat to write debugging output\n",
		        prog_name);
	} else {
		fprintf(debug_out, "%g %g %g %g\n",
		        param[0], param[1], param[2], nll[0]);
	}
	#endif

	do {
		/* Take finite differences to estimate the NLL function's
		   gradient at the current parameter estimates. */
		nll[1] = -log_like(d, param[0] + dee[0], param[1], param[2], 1);
		nll[2] = -log_like(d, param[0], param[1] + dee[1], param[2], 1);
		nll[3] = -log_like(d, param[0], param[1], param[2] + dee[2], 1);
		grad[0] = (nll[1] - nll[0]) / dee[0];
		grad[1] = (nll[2] - nll[0]) / dee[1];
		grad[2] = (nll[3] - nll[0]) / dee[2];

		/* Adjust the parameter estimates based on the gradient estimate.
		   Compute the new estimates' NLL. */
		adjust_params_using_gradients(grad, param, rss, points);
		nll[4] = -log_like(d, param[0], param[1], param[2], 1);
		delta_nll = nll[4] - nll[0];

		#ifdef DEBUG
		printf("%#.6g %#.6g %#.6g %#.6g %# 9.4f"
		       "  %# .4g %# .6g %# .8f    %#.3g %#.3g %#.3g\n",
		       param[0], param[1], param[2], nll[4], delta_nll,
		       grad[0], grad[1], grad[2], rss[0], rss[1], rss[2]);
		if (debug_out != NULL) {
			fprintf(debug_out, "%g %g %g %g\n",
			        param[0], param[1], param[2], nll[4]);
		}
		#endif

		if (delta_nll < 0.0) {
			/* This gradient descent step was successful! */
			rss[0] *= 1.5;
			rss[1] *= 1.5;
			rss[2] *= 1.5;
			nll[0] = nll[4];  /* already computed the new values' NLL */
		} else {
			/* This gradient descent step failed. */
			param[0] += rss[0] * grad[0];
			param[1] += rss[1] * grad[1];
			param[2] += rss[2] * grad[2];
			rss[0] /= 2.0;
			rss[1] /= 2.0;
			rss[2] /= 2.0;
		}
		iterations++;
	} while ((delta_nll > 0.0) || (delta_nll < -8e-5));

	*final_nll = nll[0];

	#ifdef DEBUG
	if (debug_out != NULL) {
		fclose(debug_out);
	}
	#endif

	return iterations;
}

void ses_from_log_like(const Data* d, double* param, double final_nll, double* ses)
{
	double cent_diff;
	double dee[3] = { param[0] / 1e4, param[1] / 1e4, param[2] / 1e4 };
	unsigned int i;
	double nll[2];

	for (i = 0; i < 3; i++) {
		/* Approximate the second derivative of the negative log likelihood
		   function using a central difference. */
		param[i] -= dee[i];
		nll[0] = -log_like(d, param[0], param[1], param[2], 1);
		param[i] += (2.0 * dee[i]);
		nll[1] = -log_like(d, param[0], param[1], param[2], 1);
		cent_diff = nll[1] - (2.0 * final_nll) + nll[0];
		ses[i] = 1.0 / (sqrt(cent_diff) / dee[i]);
		param[i] -= dee[i];
	}
}

#ifdef ETA
void write_eta(FILE* out, const Data* d, const double* param, double actual_a, double start_time)
{
	double dt;
	unsigned int i;
	const double mu = sqrt(param[2]);
	const double mult = QE * QE * param[0] * param[1] / (EPSILON0 * sqrt(2 * PI * ME * KB));
	double theta;
	double upsilon;
	double y = 0.0;

	for (i = 0; i < d->N; i++) {
		if (d->t[i] < start_time) {
			continue;
		}
		theta = d->T_i[i] / d->T_e[i];
		upsilon = mult / sqrt(d->T_e[i]);
		if (!i) {
			dt = d->t[0] - start_time;
		} else {
			dt = d->t[i] - d->t[i-1];
		}
		y += dt * upsilon * (sqrt(theta) * (1.0 - (y / theta)) / mu - exp(y));
		fprintf(out, "%g %g\n", d->t[i], y * param[0] / actual_a);
	}
}
#endif

int main(int argc, char* argv[])
{
	Data d;
	#ifdef ETA
	FILE* eta_out;
	#endif
	double final_nll;
	int i;
	double initial_nll;
	unsigned int iterations;
	double new_param[3];
	double param[3];  /* a, n_e, and m_i / m_e respectively */
	double prev_nll;
	double prev_param[3];
	double ses[3];

	if (argc < 5) {
		fprintf(stderr,
		        "Usage: %s [A] [N_E] [START TIME] [FILES]\n", argv[0]);
		return EXIT_FAILURE;
	}

	init_Data(&d);

	/* I could use `strtod` or `sscanf` instead of `atof` in this function
	   to be really fastidious but it scarcely seems necessary. */
	for (i = 4; i < argc; i++) {
		read_macroscopic_data_file(&d, argv[0], argv[i], atof(argv[3]));
		if (!d.N) {
			return EXIT_FAILURE;
		}
	}
	if (!(d.any_charging)) {
		fprintf(stderr, "%s: no grain charging occurred\n", argv[0]);
	}

	param[0] = atof(argv[1]);
	param[1] = atof(argv[2]);
	param[2] = 1836.15;

	prev_param[0] = param[0];
	prev_param[1] = param[1];
	prev_param[2] = param[2];

	iterations = fit_by_grad_desc(argv[0], &d, param,
	                              &initial_nll, &final_nll);
	printf("Initial NLL: %g (for %u data points)\n", initial_nll, d.N);
	if (!isfinite(initial_nll)) {
		fprintf(stderr, "%s: bad initial NLL %g, failing\n",
		        argv[0], initial_nll);
		return EXIT_FAILURE;
	}

	while (1) {
		new_param[0] = param[0] + (0.8 * (param[0] - prev_param[0]));
		if (new_param[0] < 0.0) {
			new_param[0] = param[0];
		}
		new_param[1] = param[1] + (0.8 * (param[1] - prev_param[1]));
		if (new_param[1] < 0.0) {
			new_param[1] = param[1];
		}
		new_param[2] = param[2] + (0.8 * (param[2] - prev_param[2]));
		if (new_param[2] < 0.0) {
			new_param[2] = param[2];
		}
		prev_nll = final_nll;
		prev_param[0] = param[0];
		prev_param[1] = param[1];
		prev_param[2] = param[2];
		iterations += fit_by_grad_desc(argv[0], &d, new_param,
		                               &initial_nll, &final_nll);
		if ((final_nll - prev_nll) > -0.01) {
			final_nll = prev_nll;
			break;
		}
		param[0] = new_param[0];
		param[1] = new_param[1];
		param[2] = new_param[2];
	}

	ses_from_log_like(&d, param, final_nll, ses);

	puts("  a           n_e         m_i / m_e   NLL      iterations");
	printf("  %11g %11g %11g %8g %u\n",
	       param[0], param[1], param[2], final_nll, iterations);
	printf("\xc2\xb1 %11.3g %11.3g %11.3g\n", ses[0], ses[1], ses[2]);

	#ifdef ETA
	if ((eta_out = fopen(ETA, "wb")) == NULL) {
		fprintf(stderr, "%s: can't open %s to write eta\n", argv[0], ETA);
		free_Data_memory(&d);
		return EXIT_FAILURE;
	}
	write_eta(eta_out, &d, param, atof(argv[1]), atof(argv[3]));
	fclose(eta_out);
	#endif

	free_Data_memory(&d);

	return EXIT_SUCCESS;
}
