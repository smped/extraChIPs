#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <float.h>
#include "lookup_tables.h"

/* Modified from FMStable by Geoff Robertson */
/* The specific use case here is the harmonic mean p and as such */
/* Tables and if statements which are not required have been incrementally */
/* commented out, deleted and a handful of constants defined to save recalculation */
/* The gain in performance is currently about a 15x speed-up */

/*
 * OI is the order of interpolation to be used
 * HOI is half the order of interpolation to be used
 * OIM1 is OI-1
 */
#define OI 16
#define HOI 8
#define OIM1 15

#define FALSE 0
#define TRUE 1

static const double SQRT_PI = 1.77245385090551602729816748334115; // sqrt(M_PI)
static const double xlowlimit = -44.900242184179755123; // -(1 + log(M_PI_2 * 1.0e30)) * M_2_PI

/* Define all constants & variables for table lookups */
static double previous_alpha = -999.;
static double previous_oneminusalpha = -999.;
static double previous_twominusalpha = -999.;

static double sa2, nu, Calpha_M, midpoint, xi, ximid;
static double Clogd, alphastar, eta;
static double ffound, dfound, logapprox;

/* Boolean variables stored as integers. */
static int distributiontabulated;

/* Interpolation table i gives values for density as di and for distribution */
/* function (or its complement) as fi for nxi values of x and nyi values of y*/

/* First tables are for alpha < 1, alpha*xi > 1/5 in C parametrization */
/* Vy1 is alpha, Vx1 is proportional to 1/(alpha*xi)	*/
#define nx1 70
#define ny1 20
static double f1[nx1], d1[nx1];
static double xdenom1[(nx1 - OIM1) * OI];
static double ydenom1[(ny1 - OIM1) * OI];

/* Second tables are for alpha < .5, alpha*xi < 1/5 and x<1 in C parametrization */
/* Vy2 is alpha, Vx2 is proportional to x**(-1/alpha) */
#define nx2 20
#define ny2 20
static double xdenom2[(nx2 - OIM1) * OI];
static double ydenom2[(ny2 - OIM1) * OI];

/* Third tables are for alpha < .5, x > 1 in C parametrization */
/* Vy3 is alpha, Vx3 is proportional to x**(-1/alpha) */
#define nx3 20
#define ny3 20
static double xdenom3[(nx3 - OIM1) * OI];
static double ydenom3[(ny3 - OIM1) * OI];

/* Fourth tables are for 1.7 <alpha < 2, -1.3 < x (M=S0) < 20 */
/* Use the difference from the alpha=2 (normal) distribution for all larger x,
 scaling so that the asymptotic behaviour as x tends to infinity is handled
 accurately. */
/* Vy4 is alpha, Vx4 is based on x in the S0 parametrization */
#define nx4 100
#define ny4 17
static double xdenom4[(nx4 - OIM1) * OI];
static double ydenom4[(ny4 - OIM1) * OI];
static double f4_alpha2[nx4];
static double d4_alpha2[nx4];

/* Fifth tables are for 1.7 <alpha < 2, x (M=S0) > 20 */
/* Use Zolotarev 2.5.6 */
/* Vy5 is alpha, Vx5 is proportional to y**(-1/alpha) where x=y+eta y**(1-alpha) */
#define nx5 20
#define ny5 17
static double xdenom5[(nx5 - OIM1) * OI];
static double ydenom5[(ny5 - OIM1) * OI];

/* Sixth tables are for 0.5 <alpha < 1.7, x (M=S0) > 5 */
/* Use Zolotarev 2.5.6 */
/* Vy6 is alpha, Vx6 is proportional to y**(-1/alpha) where x=y+eta y**(1-alpha) */
#define nx6 20
#define ny6 40
static double f6[nx6], d6[nx6];
static double xdenom6[(nx6 - OIM1) * OI];
static double ydenom6[(ny6 - OIM1) * OI];

/* Seventh tables are for 0.5 <alpha < 1.7, from xi=2/5 to x=7.3 */
/* 	Table 1 goes from xi=.2/alpha to xi=infinity, so regions overlap. */
/* Vy7 is alpha, Vx7 is x (S0 parametrization) */
#define nx7 60
#define ny7 40
static double f7[nx7], d7[nx7];
static double xdenom7[(nx7 - OIM1) * OI];
static double ydenom7[(ny7 - OIM1) * OI];


/*====================================================================== */
/**
 * Computes Mill's ratio, which is the ratio of the right tail
 * probability to the density for a standard normal distribution.
 *
 * @param x Input value (should be non-negative)
 * @return The computed Mill's ratio value
 */
double millsratio(double x)
{
    /* Data points for Thiele interpolation */
    static const double rhodiff[26] = {
        1.0, 3.94766755064487339, 0.262699373075432297, -4.79530994886257864,
        -0.472387694632788044, -4.64178612251027508, 0.888029797897972282E-1,
        3.34506626993363627, -0.253018663557510120, -3.22271919537673200,
        0.367907414060725404, 0.515135649489030960, -0.782205852853700406,
        -1.66658903030166903, -0.317783299853388699, 3.44266858016113165,
        0.110601063267646032, -113.675951225046565, -0.475674678528123486E-2,
        891.504346794291090, 0.155839071476740257E-3, -52392.5264048076612,
        -0.145796580254323180E-4, 5331522.33187226848, 0.169353718097630593E-6,
        -83341339.7482781260
    };

    static const double xdata[25] = {
        0.0, 1.0, 0.2704, 0.5776, 0.0784, 0.04, 0.0576, 0.0256, 0.1024,
        0.1296, 0.16, 0.1936, 0.2304, 0.0064, 0.3136, 0.36, 0.4096, 0.4624,
        0.5184, 0.0144, 0.64, 0.7056, 0.7744, 0.8464, 0.9216
    };

    /* Check for negative input */
    if (x < 0.0) {
        return 0.0;
    }

    /* Compute Thiele interpolation using iterative approach */
    double t = 1.0 / (x + 1.0);
    double result = rhodiff[25];  /* Start with the last coefficient */

    /* Calculate the continued fraction from the bottom up */
    for (int i = 24; i >= 0; i--) {
        double denominator = result;

        /* Handle the near-zero case */
        if (fabs(denominator) < DBL_EPSILON * 100.0) {
            denominator = DBL_EPSILON * 100.0;
        }

        /* For the last step (i = 0), we don't need to add (t - xdata[i]) */
        if (i == 0) {
            result = rhodiff[i];
        } else {
            result = rhodiff[i] + (t - xdata[i-1]) / denominator;
        }
    }

    return t * result;
}
/*=====================================================================*/

/* Finds the right-hand tail probability for a standard normal distribution. */
double normaltail(double z)
{
    /* Standard normal density at zero: 1/sqrt(2π) */
    static const double STANDARD_NORMAL_DENSITY = 0.39894228040143267793994605993438187;

    double val;

    if (z < 0.0) {
        /* For negative z, use the fact that P(X > z) = 1 - P(X > -z) */
        val = 1.0 - STANDARD_NORMAL_DENSITY * exp(-0.5 * z * z) * millsratio(-z);
    } else {
        /* For positive z, use Mills ratio directly */
        val = STANDARD_NORMAL_DENSITY * exp(-0.5 * z * z) * millsratio(z);
    }

    return val;
}

/*=====================================================================*/

double LogGamma(double x) {

    /* Abramowitz & Stegun equation 6.1.48 */

    static const double c0 = 0.9189385332046727417803297364056177;
    static const double a[7] = {
        1.0 / 12.0, 1.0 / 30.0, 53.0 / 210.0, 195.0 / 371.0,
        22999.0 / 22737.0, 29944523.0 / 19733142.0, 109535241009.0 / 48264275462.0
    };

    double Z, val;

    if (x > 10) {
        val = (x - 0.5) * log(x) + c0 - x;
        val += a[0] / (x + a[1] / (x + a[2] / (x + a[3] / (x + a[4] / (x + a[5] / (x + a[6] / x))))));
    } else {
        Z = x + 9;
        val = (Z - 0.5) * log(Z) + c0 - Z;
        val += a[0] / (Z + a[1] / (Z + a[2] / (Z + a[3] / (Z + a[4] / (Z + a[5] / (Z + a[6] / Z))))));
        val -= log(x * (x + 1) * (x + 2) * (x + 3) * (x + 4) * (x + 5) * (x + 6) * (x + 7) * (x + 8));
    }

    return val;
}

/*====================================================================== */
void calc_recip_denom(int nx, const double x[], double denom[]) {

    /* Calculates reciprocals of denominators for use in later interpolation */

    int i, j, k, offset;
    double product;

    for (i = 0; i < nx - OIM1; i++) {
        offset = i;

        for (j = 0; j < OI; j++) {
            product = 1.0;

            for (k = 0; k < OI; k++) {
                if (k != j) {
                    product *= (x[j + offset] - x[k + offset]);
                }
            }

            denom[i * OI + j] = 1.0 / product;
        }
    }
}
/*========================================================================= */
/**
 * Performs interpolation to find values at point x using provided data points.
 *
 * @param x The point at which to interpolate
 * @param f Pointer to store the interpolated function value
 * @param d Pointer to store the interpolated derivative value
 * @param nxn Number of data points in the arrays
 * @param xn Array of x coordinates
 * @param fn Array of function values
 * @param dn Array of derivative values
 * @param xdenomn Array of precomputed reciprocal denominators
 */
void interpolate(double x, double *f, double *d, int nxn, const double xn[],
                 double fn[], double dn[], double xdenomn[])
{
    /* Find interpolated values for x using vectors xn, fn and dn with */
    /* reciprocals of denominators specified by array xdenomn */
    double difference[OI], product, weight;
    int low, high, mid, offset, start, k;

    /* Check if x is outside the range */
    if (nxn <= 0 || x > xn[nxn - 1]) {
        return;  /* Function returns without setting f and d */
    }

    /* First find the smallest index of a larger value of xn by binary search */
    low = 0;
    high = nxn - 1;

    while (low < high) {
        mid = low + (high - low) / 2;  /* Avoids potential overflow */

        if (xn[mid] >= x) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }

    /* Calculate the optimal starting position for the interpolation window */
    start = (high >= OI) ? ((high - OI + 1 < nxn - OI) ? high - OI + 1 : nxn - OI) : 0;
    offset = start;

    /* Compute differences and their product */
    product = 1.0;
    for (k = 0; k < OI; k++) {
        difference[k] = x - xn[k + offset];
        product *= difference[k];
    }

    /* Handle the case where x exactly matches one of the nodes */
    if (fabs(product) < DBL_EPSILON) {
        for (k = 0; k < OI; k++) {
            if (fabs(x - xn[k + offset]) < DBL_EPSILON) {
                *f = fn[k + offset];
                *d = dn[k + offset];
                return;
            }
        }
    }

    /* Perform the interpolation */
    *f = 0.0;
    *d = 0.0;
    for (k = 0; k < OI; k++) {
        /* Skip division by zero */
        if (fabs(difference[k]) < DBL_EPSILON) continue;

        weight = product * xdenomn[start * OI + k] / difference[k];
        *f += weight * fn[k + offset];
        *d += weight * dn[k + offset];
    }
}
/*=================================================================== */
void interpolate_over_alpha(
        int nx, int nalpha, const double alphalist[],
        double thisalpha, const double tablef[], const double tabled[],
        double thisf[], double thisd[],double denom[]
)
{
    /* To interpolate the tables tablef and tabled over alpha */
    double weight, product, difference[OI];
    int i, j, k, start, offset;

    /* First find the smallest index of a larger value of alpha */
    /* (This could be made faster with a binary chop algorithm) */
    for (j=0; j<nalpha; j++){
        if(alphalist[j] > thisalpha)break;
    }
    start = fmin(fmax(0, j - HOI), nalpha - OI);
    offset = start;
    product = 1;
    for (k=0; k<OI; k++){
        difference[k] = thisalpha - alphalist[k+offset];
        product = product * difference[k];
    }

    /* Minor option: when thisalpha is a tabulated value */
    if(product == 0){
        for (k=0; k<OI; k++){
            if(thisalpha == alphalist[k+offset]){
                for (i=0; i<nx; i++){
                    thisf[i] = tablef[i * nalpha + (k + offset)];
                    thisd[i] = tabled[i * nalpha + (k + offset)];
                }
                break;
            }
        }
    }

    /* Minor option: when thisalpha is a tabulated value */
    if(product == 0){
        for (k = 0; k < OI; k++){
            if(thisalpha == alphalist[k + offset]){
                for (i = 0; i < nx; i++){
                    thisf[i] = tablef[i * nalpha + (k + offset)];
                    thisd[i] = tabled[i * nalpha + (k + offset)];
                }
                break;
            }
        }
    }

    /* Major option: need to interpolate */
    else{
        for (i=0; i<nx; i++){
            thisf[i]=0;
            thisd[i]=0;
        }
        for (k=0; k< OI; k++){
            weight = product * denom[start * OI + k] / difference[k];
            for (i=0; i<nx; i++){
                thisf[i] += weight * tablef[i * nalpha + (k + offset)];
                thisd[i] += weight * tabled[i * nalpha + (k + offset)];
            }
        }
    }


}
/*====================================================================== */

/*========================================================================= */
static void setalpha(double alpha, double oneminusalpha, double twominusalpha)
{
    /* functions setalpha might be made faster by initializing arrays */
    double sinangle;
    int i;

    /* If alpha is 2 or the same as the last time, then do nothing. */
    if ((twominusalpha == 0) | ((alpha == previous_alpha) &
    (oneminusalpha == previous_oneminusalpha) &
    (twominusalpha == previous_twominusalpha) ))return;

    /* On first call, read array margins as vectors and tabulated arrays */
    if (previous_alpha == -999.){

        /* Compute reciprocals of denominators for later interpolation */
        calc_recip_denom(nx1, Vx1, xdenom1);
        calc_recip_denom(ny1, Vy1, ydenom1);
        calc_recip_denom(nx2, Vx2, xdenom2);
        calc_recip_denom(ny2, Vy2, ydenom2);
        calc_recip_denom(nx3, Vx3, xdenom3);
        calc_recip_denom(ny3, Vy3, ydenom3);
        calc_recip_denom(nx4, Vx4, xdenom4);
        calc_recip_denom(ny4, Vy4, ydenom4);
        calc_recip_denom(nx5, Vx5, xdenom5);
        calc_recip_denom(ny5, Vy5, ydenom5);
        calc_recip_denom(nx6, Vx6, xdenom6);
        calc_recip_denom(ny6, Vy6, ydenom6);
        calc_recip_denom(nx7, Vx7, xdenom7);
        calc_recip_denom(ny7, Vy7, ydenom7);

        /* Also calculate Gaussian distribution for tabulated x's for use with table4 */
        for (i=0; i<nx4; i++){
            f4_alpha2[i] = normaltail(Vx4[i] * M_SQRT1_2);
            d4_alpha2[i] = 1 / (2 * SQRT_PI * exp(-Vx4[i] * Vx4[i] * .25));
        }

    }
    /* end of initialization */
    /* ===================================================== */
    /* Store standard numbers which vary with alpha */
    previous_alpha = alpha;
    previous_oneminusalpha = oneminusalpha;
    previous_twominusalpha = twominusalpha;
    distributiontabulated = FALSE;

    /* Case when alpha > .5 */
    alphastar = alpha;
    ximid = .4;
    midpoint = (-log(M_PI_2 * ximid) - 1) * M_2_PI;
    nu = 1;
    eta = 0;

    /* Lower limit where xi=10**30; take density to be zero below here */
    // xlowlimit = -(1 + log(M_PI_2 * 1.E30)) * M_2_PI; // Define as const above

    sa2 = twominusalpha / (2 * alpha);
    // Clogd = log(nu / sqrt(2 * M_PI * alpha)); // Subtract log(sqrt(2 * M_PI * alpha))??
    Clogd = log(nu) - 0.5 * log(2 * M_PI * alpha);
    sinangle = sin(M_PI_2 * twominusalpha);
    Calpha_M = exp(LogGamma(alpha)) * sinangle * M_1_PI;

    interpolate_over_alpha(nx1, ny1, Vy1, alphastar, tablef1, tabled1, f1, d1, ydenom1);
    interpolate_over_alpha(nx6, ny6, Vy6, alpha, tablef6, tabled6, f6, d6, ydenom6);
    interpolate_over_alpha(nx7, ny7, Vy7, alpha, tablef7, tabled7, f7, d7, ydenom7);
}
/*========================================================================= */
/**
 * Computes density, distribution function and complement for a maximally skew
 * stable distribution skewed to the right.
 */
void tailsMSS(int n, double x[], double cF[], double location)
{
    /* Constants */
    static const double NEWTON_TOLERANCE = 1.0e-10;

    double z, y, dy, t;
    int i;

    // Try declaring other variables here
    double F[n];
    double logF[n];
    double logcF[n];

    // Fixed parameters for this specific case
    double alpha = 1.0;
    double oneminusalpha = 0.0;
    double twominusalpha = 1.0;

    /* Input validation */
    if (n <= 0 || !x || !cF) {
        return;
    }

    /* Set up parameters for the current alpha value */
    setalpha(alpha, oneminusalpha, twominusalpha);

    /* Process each point */
    for (i = 0; i < n; i++) {
        /* Scale the input */
        z = (x[i] - location) * M_2_PI;

        /* Case 1: z below lower limit where xi can be calculated */
        if (z < xlowlimit) {
            F[i] = 0.0;
            logF[i] = -DBL_MAX;
            cF[i] = 1.0;
            logcF[i] = 0.0;
        }
        /* Case 2: Low range for x - use table 1 */
        else if (z < midpoint) {
            xi = exp(-1.0 - M_PI_2 * z) * M_2_PI;
            t = 0.2 / (alphastar * xi);

            interpolate(t, &ffound, &dfound, nx1, Vx1, f1, d1, xdenom1);

            /* Calculate and store all required values */
            logF[i] = -0.5 * log(2.0 * M_PI * alpha * xi) - xi + log(ffound);
            F[i] = exp(logF[i]);
            logcF[i] = log1p(-F[i]);
            cF[i] = 1.0 - F[i];
        }
        /* Case 3: Middle range for x - use table 7 */
        else if (z < 7.3) {
            t = (z - midpoint) / (7.3 - midpoint);

            interpolate(t, &ffound, &dfound, nx7, Vx7, f7, d7, xdenom7);

            /* Calculate and store all required values */
            logcF[i] = ffound;
            cF[i] = exp(ffound);
            F[i] = 1.0 - cF[i];
            logF[i] = log1p(-cF[i]);
        }
        /* Case 4: Upper range for x - use table 6 */
        else {
            /* Use Newton-Raphson to solve for y */
            y = z;  /* Initial guess */

            do {
                /* Newton-Raphson iteration */
                dy = (z - y - log(y) * M_2_PI) / (1.0 + 1.0 / (y * M_PI_2));
                y = y + dy;
            } while (fabs(dy) > NEWTON_TOLERANCE * y);

            t = pow((0.2 * y), (-alpha));

            interpolate(t, &ffound, &dfound, nx6, Vx6, f6, d6, xdenom6);

            /* Calculate and store all required values */
            logapprox = log(2.0 * Calpha_M) - alpha * log(y);
            logcF[i] = logapprox + log(ffound);
            cF[i] = exp(logcF[i]);
            F[i] = 1.0 - cF[i];
            logF[i] = log1p(-cF[i]);
        }
    }
}
/*====================================================================== */

SEXP RtailsMSS(SEXP Rlocation, SEXP Rx) {
    // Protect R objects from garbage collection
    PROTECT(Rlocation = coerceVector(Rlocation, REALSXP));
    PROTECT(Rx = coerceVector(Rx, REALSXP));

    // Get C values from R objects
    double location = REAL(Rlocation)[0];
    double *x = REAL(Rx);
    int n = LENGTH(Rx);

    // Create output R object for cF only
    SEXP RcF = PROTECT(allocVector(REALSXP, n));
    double *cF = REAL(RcF);

    // Call the C function
    tailsMSS(n, x, cF, location);

    // Unprotect all protected objects
    UNPROTECT(3);

    // Since we only need cF, we just return that directly
    return RcF;
}
