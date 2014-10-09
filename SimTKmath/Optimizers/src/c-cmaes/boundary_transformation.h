/*
 * INTERFACE: the relevant functions are
 *
 * * boundary_transformation_init(this, l_bound, u_bound, len)
 * * boundary_transformation_exit(this)
 * * boundary_transformation(this, x, y, len)
 *
 * implements a smooth mapping double *x -> double *y that guarantees
 * elements of y to be within specified boundaries. The mapping is piecewise
 * either linear or quadratic and can achieve values arbitrarily close to and
 * on the boundaries. The middle of the domain l_bound + (u_bound-l_bound) / 2.0
 * always maps to itself. Typically, 90% of feasible values (those not close 
 * to the boundaries) are mapped to themselves, preventing any numerical subtleties. 
 * Specifically, al, au > 0 are internally chosen offsets. The mapping
 * [l_bound - al, u_bound + au] <-> [l_bound, u_bound] is monotonous, bijective
 * and invertible. It is the identity within [l_bound + al, u_bound - au] and 
 * quadratic for [l_bound - 3*al, l_bound + al] (with l_bound - al -> l_bound)
 * and for [u_bound - au, u_bound + 3*au] (with u_bound + au -> u_bound).
 *
 * The method is robust against very small/large boundary values, say
 * -1e99 and/or 1e99, to emulated unbounded variables. In this case values
 * between -1e98 and 1e98 are never changed, i.e. mapped to itself.
 *
 */

typedef struct {
	double const *lower_bounds; /* array of size len_of_bounds */
	double const *upper_bounds; /* array of size len_of_bounds */
	unsigned long len_of_bounds; /* in case, last value is recycled */
	double *al; /* "add"-on to lower boundary preimage, same length as bounds */
	double *au; /* add-on to upper boundary preimage, same length as bounds */
} boundary_transformation_t;

/* set lower and upper bounds, the values lower_bounds[len_of_bounds - 1] and
 * upper_bounds[len_of_bounds - 1] are recycled for any element >= len_of_bounds.
 * If len_of_bounds == 0, no bounds are assumed. If len_of_bounds == 1, the
 * zero pointer is allowed for lower_bounds or upper_bounds and indicates no
 * respective bounds. "no bounds" is "emulated" using the very small/large value
 * of DBL_MAX / -1e2 and DBL_MAX / 1e2, respectively. */
void boundary_transformation_init(boundary_transformation_t *,
		double const *lower_bounds, double const *upper_bounds, unsigned long len_of_bounds);

/* release memory */
void boundary_transformation_exit(boundary_transformation_t *);

/* on return, y is guaranteed to have all values within the boundaries.
 * The caller inputs x and is responsible for having allocated y in that
 * y[len-1] = x[len-1] is a valid operation.  x==y is valid input, but
 * will fail together with cmaes when x is an element of the population
 * returned by cmaes_SamplePopulation (these elements are of type
 * double const * for a reason).
 * */
void boundary_transformation(boundary_transformation_t *,
		double const *x, double *y, unsigned long len); /* new value into y */

/* after
 *   boundary_transformation(b,x,y,l) ;
 * the two consecutive calls
 *   boundary_transformation_inverse(b,y,x,l) ; boundary_transformation(b,x,y,l) ;
 * have no effect on y anymore (but they might change x!).
 * */
void boundary_transformation_inverse(boundary_transformation_t *t,
		double const *y, double *x, unsigned long len); /* new value into x */

/* used by function boundary_transformation. After applying the shift,
 *   boundary_transformation_shift_into_feasible_preimage(b,x,x,l)
 * the two consecutive calls
 *   boundary_transformation(b,x,y,l) ; boundary_transformation_inverse(b,y,x,l) ;
 * have no effect on x anymore */
void boundary_transformation_shift_into_feasible_preimage(boundary_transformation_t *t,
		double const *x, double *x_shifted, unsigned long len); /* new value into x_shifted */


