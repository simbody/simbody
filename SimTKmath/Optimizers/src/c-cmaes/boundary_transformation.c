#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <float.h> /* DBL_MAX */
#include <stdio.h>
#include "boundary_transformation.h"

static int do_assertions = 1;
static unsigned long _index(boundary_transformation_t *t, unsigned long i);
static void _FatalError(const char *s);
static double default_lower[1];
static double default_upper[1];

void boundary_transformation_init(boundary_transformation_t *t,
		double const *lower_bounds, double const *upper_bounds, unsigned long len_of_bounds)
{
	unsigned i;
	double const *ub, *lb;
	default_lower[0] = DBL_MAX / -1e2;
	default_upper[0] = DBL_MAX / 1e2;
	t->lower_bounds = lower_bounds;
	t->upper_bounds = upper_bounds;
	t->len_of_bounds = len_of_bounds;

	if (lower_bounds == NULL && len_of_bounds <= 1) /* convenience default */
		t->lower_bounds = default_lower;
	if (upper_bounds == NULL && len_of_bounds <= 1)
		t->upper_bounds = default_upper;
	if (len_of_bounds == 0) {
		t->lower_bounds = default_lower;
		t->upper_bounds = default_upper;
		t->len_of_bounds = 1;
	}

	if (t->lower_bounds == NULL || t->upper_bounds == NULL)
		_FatalError("init: input upper_bounds or lower_bounds was NULL and len_of_bounds > 1");

	/* compute boundaries in pre-image space, al and au */
	t->al = calloc(t->len_of_bounds, sizeof(double));
	t->au = calloc(t->len_of_bounds, sizeof(double));
	if (!t->al || !t->au)
		_FatalError(" in _init(): could not allocate memory");

	lb = t->lower_bounds;
	ub = t->upper_bounds;
	for(i = 0; i < t->len_of_bounds; ++i) {
		if (lb[i] == ub[i])
			_FatalError("in _init: lower and upper bounds must be different in all variables");
		/* between lb+al and ub-au transformation is the identity */
		t->al[i] = fmin((ub[i] - lb[i]) / 2., (1. + fabs(lb[i])) / 20.);
		t->au[i] = fmin((ub[i] - lb[i]) / 2., (1. + fabs(ub[i])) / 20.);
	}
}
void boundary_transformation_exit(boundary_transformation_t *t)
{
	if(t->al)
		free(t->al);
	if(t->au)
		free(t->au);
}
void boundary_transformation(boundary_transformation_t *t,
		double const *x, double *y, unsigned long len)
{
	double lb, ub, al, au;
	unsigned long i;
	boundary_transformation_shift_into_feasible_preimage(t, x, y, len);
	for(i = 0; i < len; ++i) {
		lb = t->lower_bounds[_index(t, i)];
		ub = t->upper_bounds[_index(t, i)];
		al = t->al[_index(t, i)];
		au = t->au[_index(t, i)];
        if (y[i] < lb + al)
            y[i] = lb + (y[i] - (lb - al)) * (y[i] - (lb - al)) / 4. / al;
        else if (y[i] > ub - au)
            y[i] = ub - (y[i] - (ub + au)) * (y[i] - (ub + au)) / 4. / au;
	}
}
void boundary_transformation_shift_into_feasible_preimage(
			boundary_transformation_t *t, double const *x, double *y, unsigned long len)
{
	double lb, ub, al, au, r, xlow, xup;
	unsigned long i;

	for(i = 0; i < len; ++i) {
		lb = t->lower_bounds[_index(t, i)];
		ub = t->upper_bounds[_index(t, i)];
		al = t->al[_index(t, i)];
		au = t->au[_index(t, i)];
		xlow = lb - 2 * al - (ub - lb) / 2.0;
		xup = ub + 2 * au + (ub - lb) / 2.0;
        r = 2 * (ub - lb + al + au); /* == xup - xlow == period of the transformation */

        y[i] = x[i];

		if (y[i] < xlow) { /* shift up */
			y[i] += r * (1 + (int)((xlow - y[i]) / r));
		}
		if (y[i] > xup) { /* shift down */
			y[i] -= r * (1 + (int)((y[i] - xup) / r));
			/* printf(" \n%f\n", fmod(y[i] - ub - au, r)); */
		}
		if (y[i] < lb - al) /* mirror */
			y[i] += 2 * (lb - al - y[i]);
        if (y[i] > ub + au)
        	y[i] -= 2 * (y[i] - ub - au);

		if ((y[i] < lb - al - 1e-15) || (y[i] > ub + au + 1e-15)) {
			printf("BUG in boundary_transformation_shift_into_feasible_preimage: lb=%f, ub=%f, al=%f au=%f, y=%f\n",
					lb, ub, al, au, y[i]);
			_FatalError("BUG");
		}
	}
}
void boundary_transformation_inverse(boundary_transformation_t *t,
		double const *x, double *y, unsigned long len)
{
	double lb, ub, al, au;
	unsigned long i;

	for (i = 0; i < len; ++i) {
		lb = t->lower_bounds[_index(t, i)];
		ub = t->upper_bounds[_index(t, i)];
		al = t->al[_index(t, i)];
		au = t->au[_index(t, i)];
		y[i] = x[i];
		if (y[i] < lb + al)
			y[i] = (lb - al) + 2 * sqrt(al * (y[i] - lb));
		else if (y[i] > ub - au)
			y[i] = (ub + au) - 2 * sqrt(au * (ub - y[i]));
	}
	if (11 < 3 || do_assertions) {
		double *z = calloc(len, sizeof(double));
		for (i = 0; i < len; ++i)
			z[i] = y[i];
		boundary_transformation(t, z, y, len);
		for (i = 0; i < len; ++i)
			if (fabs(y[i] - x[i]) > 1e-14)
				printf("  difference for index %ld should be zero, is %f ", i, y[i] - x[i]);
		for (i = 0; i < len; ++i)
			y[i] = z[i];
		free(z);
	}
}
/*
static void _manage_len(boundary_transformation_t *t, unsigned int len)
{
	if (t->len_of_return_value < len) {
		if (t->return_value)
			free(t->return_value);
		t->return_value = calloc(len, sizeof(double));
		t->len_of_return_value = len;
	}
	if (!t->return_value)
		_FatalError("could not allocate memory");
}
*/
static unsigned long _index(boundary_transformation_t *t, unsigned long i)
{
	return i < t->len_of_bounds ? i : t->len_of_bounds - 1;
}
static void _FatalError(const char *s)
{
	printf("Fatal error in boundary_transformation: %s\n", s);
	exit(1);
}

