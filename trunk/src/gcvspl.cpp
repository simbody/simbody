/* gcvspl.f -- translated by f2c (version of 16 February 1991  0:35:15).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "SimTKcommon.h"

int gcvspl_(const SimTK_Real *, const SimTK_Real *, int *, const SimTK_Real *, const SimTK_Real *, int *, int *,
	 		int *, int *, SimTK_Real *, SimTK_Real *, int *, SimTK_Real *, int *);

SimTK_Real splder_(int *, int *, int *, SimTK_Real *, const SimTK_Real *, const SimTK_Real *, int *, SimTK_Real *);

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


/* Table of constant values */

static SimTK_Real c_b6 = 1e-15;

 int gcvspl_(const SimTK_Real *x, const SimTK_Real *y, int *ny, 
	const SimTK_Real *wx, const SimTK_Real *wy, int *m, int *n, int *k, 
	int *md, SimTK_Real *val, SimTK_Real *c, int *nc, SimTK_Real *
	wk, int *ier)
{

    static int m2 = 0;
    static int nm1 = 0;
    static SimTK_Real el = 0.;

    /* System generated locals */
    int y_dim1, y_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    static int nm2m1, nm2p1;
    extern SimTK_Real splc_(int *, int *, int *, const SimTK_Real *, 
	    int *, const SimTK_Real *, const SimTK_Real *, int *, SimTK_Real *, 
	    SimTK_Real *, SimTK_Real *, SimTK_Real *, int *, SimTK_Real *,
	     SimTK_Real *, SimTK_Real *, SimTK_Real *, SimTK_Real *);
    extern /* Subroutine */ int prep_(int *, int *, const SimTK_Real *, 
	    const SimTK_Real *, SimTK_Real *, SimTK_Real *);
    static int i, j;
    static SimTK_Real alpha;
    extern /* Subroutine */ int basis_(int *, int *, const SimTK_Real *, 
	    SimTK_Real *, SimTK_Real *, SimTK_Real *);
    static SimTK_Real r1, r2, r3, r4;
    static int ib;
    static SimTK_Real gf2, gf1, gf3, gf4;
    static int iwe;
    static SimTK_Real err;


    /* Parameter adjustments */
    --wk;
    
    c_dim1 = *nc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    
    --wy;
    --wx;
    
    y_dim1 = *ny;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    
    --x;

    /* Function Body */


    *ier = 0;
    if (abs(*md) > 4 || *md == 0 || abs(*md) == 1 && *val < 0. || abs(*md) == 
	    3 && *val < 0. || abs(*md) == 4 && (*val < 0. || *val > (
	    SimTK_Real) (*n - *m))) {
	*ier = 3;
	return 0;
    }
    if (*md > 0) {
	m2 = *m << 1;
	nm1 = *n - 1;
    } else {
	if (m2 != *m << 1 || nm1 != *n - 1) {
	    *ier = 3;
	    return 0;
	}
    }
    if (*m <= 0 || *n < m2) {
	*ier = 1;
	return 0;
    }
    if (wx[1] <= 0.) {
	*ier = 2;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if (wx[i] <= 0. || x[i - 1] >= x[i]) {
	    *ier = 2;
	}
	if (*ier != 0) {
	    return 0;
	}
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	if (wy[j] <= 0.) {
	    *ier = 2;
	}
	if (*ier != 0) {
	    return 0;
	}
    }


    nm2p1 = *n * (m2 + 1);
    nm2m1 = *n * (m2 - 1);
    ib = nm2p1 + 7;
    iwe = ib + nm2m1;

    if (*md > 0) {
	basis_(m, n, &x[1], &wk[ib], &r1, &wk[7]);
	prep_(m, n, &x[1], &wx[1], &wk[iwe], &el);
	el /= r1;
    }
    if (abs(*md) != 1) {
	goto L20;
    }
    r1 = *val;
    goto L100;


L20:
    if (*md < -1) {
	r1 = wk[4];
    } else {
	r1 = 1. / el;
    }
    r2 = r1 * 2.;
    gf2 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r2, &
	    c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]);
L40:
    gf1 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r1, &
	    c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]);
    if (gf1 > gf2) {
	goto L50;
    }
    if (wk[4] <= 0.) {
	goto L100;
    }
    r2 = r1;
    gf2 = gf1;
    r1 /= 2.;
    goto L40;
L50:
    r3 = r2 * 2.;
L60:
    gf3 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r3, &
	    c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]);
    if (gf3 > gf2) {
	goto L70;
    }
    if (wk[4] >= 999999999999999.88) {
	goto L100;
    }
    r2 = r3;
    gf2 = gf3;
    r3 *= 2.;
    goto L60;
L70:
    r2 = r3;
    gf2 = gf3;
    alpha = (r2 - r1) / 1.618033983;
    r4 = r1 + alpha;
    r3 = r2 - alpha;
    gf3 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r3, &
	    c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]);
    gf4 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r4, &
	    c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]);
L80:
    if (gf3 <= gf4) {
	r2 = r4;
	gf2 = gf4;
	err = (r2 - r1) / (r1 + r2);
	if (err * err + 1. == 1. || err <= 1e-6) {
	    goto L90;
	}
	r4 = r3;
	gf4 = gf3;
	alpha /= 1.618033983;
	r3 = r2 - alpha;
	gf3 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r3, &
		c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]
		);
    } else {
	r1 = r3;
	gf1 = gf3;
	err = (r2 - r1) / (r1 + r2);
	if (err * err + 1. == 1. || err <= 1e-6) {
	    goto L90;
	}
	r3 = r4;
	gf3 = gf4;
	alpha /= 1.618033983;
	r4 = r1 + alpha;
	gf4 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r4, &
		c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]
		);
    }
    goto L80;
L90:
    r1 = (r1 + r2) * .5;


L100:
    gf1 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r1, &
	    c_b6, &c[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7]);

    return 0;
}




/* BASIS.FOR, 1985-06-03 */

int basis_(int *m, int *n, const SimTK_Real *x, SimTK_Real 
	*b, SimTK_Real *bl, SimTK_Real *q)
{
    /* System generated locals */
    int b_dim1, b_offset, q_offset, i__1, i__2, i__3, i__4;
    SimTK_Real d__1;

    /* Local variables */
    static int nmip1, i, j, k, l;
    static SimTK_Real u, v, y;
    static int j1, j2, m2, ir, mm1, mp1;
    static SimTK_Real arg;



    /* Parameter adjustments */
    q_offset = 1 - *m;
    q -= q_offset;
    
    b_dim1 = *m - 1 - (1 - *m) + 1;
    b_offset = 1 - *m + b_dim1;
    b -= b_offset;
    
    --x;

    if (*m == 1) {
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    b[i * b_dim1] = 1.;
	}
	*bl = 1.;
	return 0;
    }

    mm1 = *m - 1;
    mp1 = *m + 1;
    m2 = *m << 1;
    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *m;
	for (j = -mm1; j <= i__2; ++j) {
	    q[j] = 0.;
	}
	q[mm1] = 1.;
	if (l != 1 && l != *n) {
	    q[mm1] = 1. / (x[l + 1] - x[l - 1]);
	}
	arg = x[l];
	i__2 = m2;
	for (i = 3; i <= i__2; ++i) {
	    ir = mp1 - i;
	    v = q[ir];
	    if (l < i) {
		i__3 = i;
		for (j = l + 1; j <= i__3; ++j) {
		    u = v;
		    v = q[ir + 1];
		    q[ir] = u + (x[j] - arg) * v;
		    ++ir;
		}
	    }
	    i__3 = l - i + 1;
	    j1 = max(i__3,1);
	    i__3 = l - 1, i__4 = *n - i;
	    j2 = min(i__3,i__4);
	    if (j1 <= j2) {
		if (i < m2) {
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			y = x[i + j];
			u = v;
			v = q[ir + 1];
			q[ir] = u + (v - u) * (y - arg) / (y - x[j]);
			++ir;
		    }
		} else {
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			u = v;
			v = q[ir + 1];
			q[ir] = (arg - x[j]) * u + (x[i + j] - arg) * v;
			++ir;
		    }
		}
	    }
	    nmip1 = *n - i + 1;
	    if (nmip1 < l) {
		i__3 = l - 1;
		for (j = nmip1; j <= i__3; ++j) {
		    u = v;
		    v = q[ir + 1];
		    q[ir] = (arg - x[j]) * u + v;
		    ++ir;
		}
	    }
	}
	i__2 = mm1;
	for (j = -mm1; j <= i__2; ++j) {
	    b[j + l * b_dim1] = q[j];
	}
    }

    i__1 = mm1;
    for (i = 1; i <= i__1; ++i) {
	i__2 = mm1;
	for (k = i; k <= i__2; ++k) {
	    b[-k + i * b_dim1] = 0.;
	    b[k + (*n + 1 - i) * b_dim1] = 0.;
	}
    }

    *bl = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = mm1;
	for (k = -mm1; k <= i__2; ++k) {
	    *bl += (d__1 = b[k + i * b_dim1], abs(d__1));
	}
    }
    *bl /= *n;

    return 0;
}



/* PREP.FOR, 1985-07-04 */


int prep_(int *m, int *n, const SimTK_Real *x, const SimTK_Real *
	w, SimTK_Real *we, SimTK_Real *el)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    SimTK_Real d__1;

    /* Local variables */
    static SimTK_Real f;
    static int i, j, k, l;
    static SimTK_Real y, f1;
    static int i1, i2, m2;
    static SimTK_Real ff;
    static int jj, jm, kl, nm, ku;
    static SimTK_Real wi;
    static int n2m, mp1, i2m1, inc, i1p1, m2m1, m2p1;



/* WE(-M:M,N) */
    /* Parameter adjustments */
    --we;
    --w;
    --x;

    /* Function Body */
    m2 = *m << 1;
    mp1 = *m + 1;
    m2m1 = m2 - 1;
    m2p1 = m2 + 1;
    nm = *n - *m;
    f1 = -1.;
    if (*m != 1) {
	i__1 = *m;
	for (i = 2; i <= i__1; ++i) {
	    f1 = -f1 * i;
	}
	i__1 = m2m1;
	for (i = mp1; i <= i__1; ++i) {
	    f1 *= i;
	}
    }

    i1 = 1;
    i2 = *m;
    jm = mp1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	inc = m2p1;
	if (j > nm) {
	    f1 = -f1;
	    f = f1;
	} else {
	    if (j < mp1) {
		inc = 1;
		f = f1;
	    } else {
		f = f1 * (x[j + *m] - x[j - *m]);
	    }
	}
	if (j > mp1) {
	    ++i1;
	}
	if (i2 < *n) {
	    ++i2;
	}
	jj = jm;
	ff = f;
	y = x[i1];
	i1p1 = i1 + 1;
	i__2 = i2;
	for (i = i1p1; i <= i__2; ++i) {
	    ff /= y - x[i];
	}
	we[jj] = ff;
	jj += m2;
	i2m1 = i2 - 1;
	if (i1p1 <= i2m1) {
	    i__2 = i2m1;
	    for (l = i1p1; l <= i__2; ++l) {
		ff = f;
		y = x[l];
		i__3 = l - 1;
		for (i = i1; i <= i__3; ++i) {
		    ff /= y - x[i];
		}
		i__3 = i2;
		for (i = l + 1; i <= i__3; ++i) {
		    ff /= y - x[i];
		}
		we[jj] = ff;
		jj += m2;
	    }
	}
	ff = f;
	y = x[i2];
	i__2 = i2m1;
	for (i = i1; i <= i__2; ++i) {
	    ff /= y - x[i];
	}
	we[jj] = ff;
	jj += m2;
	jm += inc;
    }

    kl = 1;
    n2m = m2p1 * *n + 1;
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	ku = kl + *m - i;
	i__2 = ku;
	for (k = kl; k <= i__2; ++k) {
	    we[k] = 0.;
	    we[n2m - k] = 0.;
	}
	kl += m2p1;
    }

    jj = 0;
    *el = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	wi = w[i];
	i__2 = m2p1;
	for (j = 1; j <= i__2; ++j) {
	    ++jj;
	    we[jj] /= wi;
	    *el += (d__1 = we[jj], abs(d__1));
	}
    }
    *el /= *n;

    return 0;
}



/* SPLC.FOR, 1985-12-12 */

SimTK_Real splc_(int *m, int *n, int *k, const SimTK_Real *y, int *
	ny, const SimTK_Real *wx, const SimTK_Real *wy, int *mode, SimTK_Real *val, 
	SimTK_Real *p, SimTK_Real *eps, SimTK_Real *c, int *nc, 
	SimTK_Real *stat, SimTK_Real *b, SimTK_Real *we, SimTK_Real *el, 
	SimTK_Real *bwe)
{
    /* System generated locals */
    int y_dim1, y_offset, c_dim1, c_offset, b_dim1, b_offset, we_dim1, 
	    we_offset, bwe_dim1, bwe_offset, i__1, i__2, i__3, i__4;
    SimTK_Real ret_val, d__1;

    /* Local variables */
    static int i, j, l;
    extern SimTK_Real trinv_(SimTK_Real *, SimTK_Real *, int *, int *)
	    ;
    static SimTK_Real dp;
    static int km;
    static SimTK_Real dt;
    static int kp;
    extern /* Subroutine */ int bandet_(SimTK_Real *, int *, int *), 
	    bansol_(SimTK_Real *, const SimTK_Real *, int *, SimTK_Real *, 
	    int *, int *, int *, int *);
    static SimTK_Real pel, esn, trn;



/* ***  Check on p-value */

    /* Parameter adjustments */
    bwe_dim1 = *m - (-(*m)) + 1;
    bwe_offset = -(*m) + bwe_dim1;
    bwe -= bwe_offset;
    
    we_dim1 = *m - (-(*m)) + 1;
    we_offset = -(*m) + we_dim1;
    we -= we_offset;
    
    b_dim1 = *m - 1 - (1 - *m) + 1;
    b_offset = 1 - *m + b_dim1;
    b -= b_offset;
    
    --stat;
    
    c_dim1 = *nc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    
    --wy;
    --wx;
    
    y_dim1 = *ny;
    y_offset = y_dim1 + 1;
    y -= y_offset;

    /* Function Body */
    dp = *p;
    stat[4] = *p;
    pel = *p * *el;
    if (pel < *eps) {
	dp = *eps / *el;
	stat[4] = 0.;
    }
    if (pel * *eps > 1.) {
	dp = 1. / (*el * *eps);
	stat[4] = dp;
    }

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *m, i__3 = i - 1;
	km = -min(i__2,i__3);
	i__2 = *m, i__3 = *n - i;
	kp = min(i__2,i__3);
	i__2 = kp;
	for (l = km; l <= i__2; ++l) {
	    if (abs(l) == *m) {
		bwe[l + i * bwe_dim1] = dp * we[l + i * we_dim1];
	    } else {
		bwe[l + i * bwe_dim1] = b[l + i * b_dim1] + dp * we[l + i * 
			we_dim1];
	    }
	}
    }

    bandet_(&bwe[bwe_offset], m, n);
    bansol_(&bwe[bwe_offset], &y[y_offset], ny, &c[c_offset], nc, m, n, k);
    stat[3] = trinv_(&we[we_offset], &bwe[bwe_offset], m, n) * dp;
    trn = stat[3] / *n;

    esn = 0.;
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    dt = -y[i + j * y_dim1];
	    i__3 = *m - 1, i__4 = i - 1;
	    km = -min(i__3,i__4);
	    i__3 = *m - 1, i__4 = *n - i;
	    kp = min(i__3,i__4);
	    i__3 = kp;
	    for (l = km; l <= i__3; ++l) {
		dt += b[l + i * b_dim1] * c[i + l + j * c_dim1];
	    }
	    esn += dt * dt * wx[i] * wy[j];
	}
    }
    esn /= *n * *k;

    stat[6] = esn / trn;
    stat[1] = stat[6] / trn;
    stat[2] = esn;
    if (abs(*mode) != 3) {
	stat[5] = stat[6] - esn;
	if (abs(*mode) == 1) {
	    ret_val = 0.;
	}
	if (abs(*mode) == 2) {
	    ret_val = stat[1];
	}
	if (abs(*mode) == 4) {
	    ret_val = (d__1 = stat[3] - *val, abs(d__1));
	}
    } else {
	stat[5] = esn - *val * (trn * 2. - 1.);
	ret_val = stat[5];
    }

    return ret_val;
}



/* BANDET.FOR, 1985-06-03 */

int bandet_(SimTK_Real *e, int *m, int *n)
{
    /* System generated locals */
    int e_dim1, e_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int i, k, l;
    static SimTK_Real di, dl;
    static int mi, km, lm;
    static SimTK_Real du;



    /* Parameter adjustments */
    
    e_dim1 = *m - (-(*m)) + 1;
    e_offset = -(*m) + e_dim1;
    e -= e_offset;

    /* Function Body */
    if (*m <= 0) {
	return 0;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	di = e[i * e_dim1];
	i__2 = *m, i__3 = i - 1;
	mi = min(i__2,i__3);
	if (mi >= 1) {
	    i__2 = mi;
	    for (k = 1; k <= i__2; ++k) {
		di -= e[-k + i * e_dim1] * e[k + (i - k) * e_dim1];
	    }
	    e[i * e_dim1] = di;
	}
	i__2 = *m, i__3 = *n - i;
	lm = min(i__2,i__3);
	if (lm >= 1) {
	    i__2 = lm;
	    for (l = 1; l <= i__2; ++l) {
		dl = e[-l + (i + l) * e_dim1];
		i__3 = *m - l, i__4 = i - 1;
		km = min(i__3,i__4);
		if (km >= 1) {
		    du = e[l + i * e_dim1];
		    i__3 = km;
		    for (k = 1; k <= i__3; ++k) {
			du -= e[-k + i * e_dim1] * e[l + k + (i - k) * e_dim1]
				;
			dl -= e[-l - k + (l + i) * e_dim1] * e[k + (i - k) * 
				e_dim1];
		    }
		    e[l + i * e_dim1] = du;
		}
		e[-l + (i + l) * e_dim1] = dl / di;
	    }
	}
    }

    return 0;
} 



/* BANSOL.FOR, 1985-12-12 */

int bansol_(SimTK_Real *e, const SimTK_Real *y, int *ny, 
	SimTK_Real *c, int *nc, int *m, int *n, int *k)
{
    /* System generated locals */
    int e_dim1, e_offset, y_dim1, y_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static SimTK_Real d;
    static int i, j, l, mi, nm1;



/* ***  Check on special cases: M=0, M=1, M>1 */

    /* Parameter adjustments */
    
    c_dim1 = *nc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    
    y_dim1 = *ny;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    
    e_dim1 = *m - (-(*m)) + 1;
    e_offset = -(*m) + e_dim1;
    e -= e_offset;

    /* Function Body */
    nm1 = *n - 1;
    if ((i__1 = *m - 1) < 0) {
	goto L10;
    } else if (i__1 == 0) {
	goto L40;
    } else {
	goto L80;
    }

L10:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
	    c[i + j * c_dim1] = y[i + j * y_dim1] / e[i * e_dim1];
	}
    }
    return 0;

L40:
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	c[j * c_dim1 + 1] = y[j * y_dim1 + 1];
	i__2 = *n;
	for (i = 2; i <= i__2; ++i) {
	    c[i + j * c_dim1] = y[i + j * y_dim1] - e[i * e_dim1 - 1] * c[i - 
		    1 + j * c_dim1];
	}
	c[*n + j * c_dim1] /= e[*n * e_dim1];
	for (i = nm1; i >= 1; --i) {
	    c[i + j * c_dim1] = (c[i + j * c_dim1] - e[i * e_dim1 + 1] * c[i 
		    + 1 + j * c_dim1]) / e[i * e_dim1];
	}
    }
    return 0;


L80:
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	c[j * c_dim1 + 1] = y[j * y_dim1 + 1];
	i__2 = *n;
	for (i = 2; i <= i__2; ++i) {
	    i__3 = *m, i__4 = i - 1;
	    mi = min(i__3,i__4);
	    d = y[i + j * y_dim1];
	    i__3 = mi;
	    for (l = 1; l <= i__3; ++l) {
		d -= e[-l + i * e_dim1] * c[i - l + j * c_dim1];
	    }
	    c[i + j * c_dim1] = d;
	}
	c[*n + j * c_dim1] /= e[*n * e_dim1];
	for (i = nm1; i >= 1; --i) {
	    i__2 = *m, i__3 = *n - i;
	    mi = min(i__2,i__3);
	    d = c[i + j * c_dim1];
	    i__2 = mi;
	    for (l = 1; l <= i__2; ++l) {
		d -= e[l + i * e_dim1] * c[i + l + j * c_dim1];
	    }
	    c[i + j * c_dim1] = d / e[i * e_dim1];
	}
    }
    return 0;
}



/* TRINV.FOR, 1985-06-03 */


SimTK_Real trinv_(SimTK_Real *b, SimTK_Real *e, int *m, int *n)
{
    /* System generated locals */
    int b_dim1, b_offset, e_dim1, e_offset, i__1, i__2, i__3;
    SimTK_Real ret_val;

    /* Local variables */
    static int i, j, k;
    static SimTK_Real dd, dl;
    static int mi;
    static SimTK_Real du;
    static int mn, mp;



/* ***  Assess central 2*M+1 bands of E**-1 and store in array E */

    /* Parameter adjustments */
    
    e_dim1 = *m - (-(*m)) + 1;
    e_offset = -(*m) + e_dim1;
    e -= e_offset;
    
    b_dim1 = *m - (-(*m)) + 1;
    b_offset = -(*m) + b_dim1;
    b -= b_offset;

    /* Function Body */
    e[*n * e_dim1] = 1. / e[*n * e_dim1];
    for (i = *n - 1; i >= 1; --i) {
	i__1 = *m, i__2 = *n - i;
	mi = min(i__1,i__2);
	dd = 1. / e[i * e_dim1];
	i__1 = mi;
	for (k = 1; k <= i__1; ++k) {
	    e[k + *n * e_dim1] = e[k + i * e_dim1] * dd;
	    e[-k + e_dim1] = e[-k + (k + i) * e_dim1];
	}
	dd += dd;
	for (j = mi; j >= 1; --j) {
	    du = 0.;
	    dl = 0.;
	    i__1 = mi;
	    for (k = 1; k <= i__1; ++k) {
		du -= e[k + *n * e_dim1] * e[j - k + (i + k) * e_dim1];
		dl -= e[-k + e_dim1] * e[k - j + (i + j) * e_dim1];
	    }
	    e[j + i * e_dim1] = du;
	    e[-j + (j + i) * e_dim1] = dl;
	    dd -= e[j + *n * e_dim1] * dl + e[-j + e_dim1] * du;
	}
	e[i * e_dim1] = dd * .5;
    }

    dd = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *m, i__3 = i - 1;
	mn = -min(i__2,i__3);
	i__2 = *m, i__3 = *n - i;
	mp = min(i__2,i__3);
	i__2 = mp;
	for (k = mn; k <= i__2; ++k) {
	    dd += b[k + i * b_dim1] * e[-k + (k + i) * e_dim1];
	}
    }
    ret_val = dd;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	e[k + *n * e_dim1] = 0.;
	e[-k + e_dim1] = 0.;
    }
    return ret_val;
} 



/* SPLDER.FOR, 1985-06-11 */


SimTK_Real splder_(int *ider, int *m, int *n, SimTK_Real *t, 
	const SimTK_Real *x, const SimTK_Real *c, int *l, SimTK_Real *q)
{
    /* System generated locals */
    int i__1, i__2;
    SimTK_Real ret_val;

    /* Local variables */
    static int lk1i1;
    static SimTK_Real xjki;
    static int i, j, k;
    static SimTK_Real z;
    static int i1, j1, k1, j2, m2, ii, jj, ki, jl, lk, mi, nk, lm, ml, jm,
	     ir, ju;
    extern /* Subroutine */ int search_(int *, const SimTK_Real *, SimTK_Real *,
	     int *);
    static SimTK_Real tt;
    static int lk1, mp1, m2m1, jin, nki, npm, lk1i, nki1;



/* ***  Derivatives of IDER.ge.2*M are alway zero */

    /* Parameter adjustments */
    
    --q;
    --c;
    --x;

    /* Function Body */
    m2 = *m << 1;
    k = m2 - *ider;
    if (k < 1) {
	ret_val = 0.;
	return ret_val;
    }

    search_(n, &x[1], t, l);

    tt = *t;
    mp1 = *m + 1;
    npm = *n + *m;
    m2m1 = m2 - 1;
    k1 = k - 1;
    nk = *n - k;
    lk = *l - k;
    lk1 = lk + 1;
    lm = *l - *m;
    jl = *l + 1;
    ju = *l + m2;
    ii = *n - m2;
    ml = -(*l);
    i__1 = ju;
    for (j = jl; j <= i__1; ++j) {
	if (j >= mp1 && j <= npm) {
	    q[j + ml] = c[j - *m];
	} else {
	    q[j + ml] = 0.;
	}
    }

    if (*ider > 0) {
	jl -= m2;
	ml += m2;
	i__1 = *ider;
	for (i = 1; i <= i__1; ++i) {
	    ++jl;
	    ++ii;
	    j1 = max(1,jl);
	    j2 = min(*l,ii);
	    mi = m2 - i;
	    j = j2 + 1;
	    if (j1 <= j2) {
		i__2 = j2;
		for (jin = j1; jin <= i__2; ++jin) {
		    --j;
		    jm = ml + j;
		    q[jm] = (q[jm] - q[jm - 1]) / (x[j + mi] - x[j]);
		}
	    }
	    if (jl >= 1) {
		goto L6;
	    }
	    i1 = i + 1;
	    j = ml + 1;
	    if (i1 <= ml) {
		i__2 = ml;
		for (jin = i1; jin <= i__2; ++jin) {
		    --j;
		    q[j] = -q[j - 1];
		}
	    }
L6:
	    ;
	}
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    q[j] = q[j + *ider];
	}
    }

    if (k1 >= 1) {
	i__1 = k1;
	for (i = 1; i <= i__1; ++i) {
	    nki = nk + i;
	    ir = k;
	    jj = *l;
	    ki = k - i;
	    nki1 = nki + 1;
	    if (*l >= nki1) {
		i__2 = *l;
		for (j = nki1; j <= i__2; ++j) {
		    q[ir] = q[ir - 1] + (tt - x[jj]) * q[ir];
		    --jj;
		    --ir;
		}
	    }
	    lk1i = lk1 + i;
	    j1 = max(1,lk1i);
	    j2 = min(*l,nki);
	    if (j1 <= j2) {
		i__2 = j2;
		for (j = j1; j <= i__2; ++j) {
		    xjki = x[jj + ki];
		    z = q[ir];
		    q[ir] = z + (xjki - tt) * (q[ir - 1] - z) / (xjki - x[jj])
			    ;
		    --ir;
		    --jj;
		}
	    }
	    if (lk1i <= 0) {
		jj = ki;
		lk1i1 = 1 - lk1i;
		i__2 = lk1i1;
		for (j = 1; j <= i__2; ++j) {
		    q[ir] += (x[jj] - tt) * q[ir - 1];
		    --jj;
		    --ir;
		}
	    }
	}
    }

    z = q[k];
    if (*ider > 0) {
	i__1 = m2m1;
	for (j = k; j <= i__1; ++j) {
	    z *= j;
	}
    }
    ret_val = z;
    return ret_val;
}



/* SEARCH.FOR, 1985-06-03 */


int search_(int *n, const SimTK_Real *x, SimTK_Real *t, 
	int *l)
{
    static int il, iu;

    /* Parameter adjustments */
    
    --x;

    /* Function Body */
    
    if (*t < x[1]) {
	*l = 0;
	return 0;
    }
    if (*t >= x[*n]) {
	*l = *n;
	return 0;
    }
    *l = max(*l,1);
    if (*l >= *n) {
	*l = *n - 1;
    }

    if (*t >= x[*l]) {
	goto L5;
    }
    --(*l);
    if (*t >= x[*l]) {
	return 0;
    }

    il = 1;
L3:
    iu = *l;
L4:
    *l = (il + iu) / 2;
    if (iu - il <= 1) {
	return 0;
    }
    if (*t < x[*l]) {
	goto L3;
    }
    il = *l;
    goto L4;
L5:
    if (*t < x[*l + 1]) {
	return 0;
    }
    ++(*l);
    if (*t < x[*l + 1]) {
	return 0;
    }
    il = *l + 1;
    iu = *n;
    goto L4;

}
