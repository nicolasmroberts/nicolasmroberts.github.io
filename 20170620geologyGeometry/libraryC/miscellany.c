


/*
Copyright 2016 Joshua R. Davis

Licensed under the Apache License, Version 2.0 (the "License"); you may not 
use this file except in compliance with the License. You may obtain a copy of 
the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software 
distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
License for the specific language governing permissions and limitations under 
the License.
*/



#include <time.h>

#define DEGREE (M_PI / 180.0)

// Input cannot alias output.
// a is index of first entry of left sublist.
// b is index of first entry of right sublist.
// c is index of first entry after right sublist; c - a is length of list.
void doubleMergeSortMerge(double *in, unsigned int a, unsigned int b,
		unsigned int c, double *out) {
	// Merge the two sublists, until one runs out of material.
	int ai = a, bi = b, i = a;
	while (ai < b && bi < c)
		if (in[ai] < in[bi]) {
			out[i] = in[ai];
			ai += 1;
			i += 1;
		}
		else {
			out[i] = in[bi];
			bi += 1;
			i += 1;
		}
	// Throw in the extra entries, if any.
	while (ai < b) {
		out[i] = in[ai];
		ai += 1;
		i += 1;
	}
	while (bi < c) {
		out[i] = in[bi];
		bi += 1;
		i += 1;
	}
}

// Sorts the in-list starting at index a, going until index b - 1 (so b - a is 
// the length of the sub-list to be sorted), storing the result in the 
// corresponding entries of out, using aux for auxiliary storage.
// in, out, and aux cannot safely alias each other.
void doubleMergeSortRecursive(double *in, unsigned int a, unsigned int b,
		double *aux, double *out) {
	unsigned int c;
	if (a + 1 == b)
		out[a] = in[a];
	else {
		c = (a + b) / 2;
		doubleMergeSortRecursive(in, a, c, out, aux);
		doubleMergeSortRecursive(in, c, b, out, aux);
		doubleMergeSortMerge(aux, a, c, b, out);
	}
}

// Sorts list of n doubles, placing result in out, using auxiliary buffer.
// Does not alter in. The arrays in, out, and aux cannot alias each other.
// If aux is NULL, then this function dynamically allocates it.
// Returns 0 if no error, 1 if error (in allocating memory).
int doubleMergeSort(unsigned int n, double *in, double *aux, double *out) {
	double *auxiliary;
	if (aux == NULL) {
		auxiliary = (double *)malloc(n * sizeof(double));
		if (auxiliary == NULL) {
			fprintf(stderr, "error: doubleMergeSort could not malloc\n");
			return 1;
		}
		else {
			doubleMergeSortRecursive(in, 0, n, auxiliary, out);
			free(auxiliary);
			return 0;
		}
	}
	else {
		doubleMergeSortRecursive(in, 0, n, aux, out);
		return 0;
	}
}

// !! doesn't handle memory allocation error well
double doubleMedian(unsigned int n, double *in) {
  if (n <= 100) {
    double aux[100], out[100];
    doubleMergeSort(n, in, aux, out);
    if (n % 2 == 0)
      return out[n / 2];
    else
      return 0.5 * (out[n / 2] + out[n / 2 + 1]);
  } else {
    double *aux, *out;
    aux = (double *)malloc(2 * n * sizeof(double));
    out = &(aux[n]);
    doubleMergeSort(n, in, aux, out);
    if (n % 2 == 0)
      return out[n / 2];
    else
      return 0.5 * (out[n / 2] + out[n / 2 + 1]);
  }
}

// Returns value q from ascending-sorted list that is greater than fraction p of the values in the list.
// Approximate inverse to doublePValue.
double doubleQuantile(unsigned int n, double *sorted, double p) {
    unsigned int i;
    i = (unsigned int)round(n * p);
    if (i >= n)
        return sorted[n];
    else
        return sorted[i];
}

// Returns fraction p such that q is greater than p of the values in the given ascending-sorted list. Approximate inverse to doubleQuantile.
double doublePValue(unsigned int n, double *sorted, double q) {
    unsigned int i = 0;
    while (i < n && sorted[i] < q)
        i += 1;
    return (double)i / n;
}

void printTime(FILE *fp) {
    time_t seconds;
    seconds = time(NULL);
    fprintf(fp, "%ld", seconds);
}

double arcCos(double c) {
	if (c >= 1.0)
		return 0.0;
	else if (c <= -1.0)
		return M_PI;
	else
		return acos(c);
}

double sqRoot(double x) {
	if (x <= 0.0)
		return 0.0;
	else
		return sqrt(x);
}

void cartesianFromSpherical(double rho, double phi, double theta, double *out) {
	double sinphi;
	sinphi = sin(phi);
	out[0] = rho * sinphi * cos(theta);
	out[1] = rho * sinphi * sin(theta);
	out[2] = rho * cos(phi);
}

void sphericalFromCartesian(double *x, double *rho, double *phi, double *theta) {
	*rho = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	if (*rho == 0.0) {
		// x is the zero vector.
		*phi = 0.0;
		*theta = 0.0;
	} else {
		*phi = arcCos(x[2] / *rho);
		double rhoSinPhi = *rho * sin(*phi);
		if (rhoSinPhi == 0.0)
			if (x[2] >= 0.0) {
				// x points toward the north pole.
				*rho = x[2];
				*phi = 0.0;
				*theta = 0.0;
			} else {
				// x points toward the south pole.
				*rho = -x[2];
				*phi = M_PI;
				*theta = 0.0;
			}
		else {
			// This is the typical case.
			*theta = atan2(x[1], x[0]);
		}
	}
}

void sphericalFromTrendPlungeDeg(double trendDeg, double plungeDeg, double *phi, double *theta) {
	*phi = (plungeDeg + 90.0) * DEGREE;
	*theta = (90.0 - trendDeg) * DEGREE;
	while (*theta >= 2.0 * M_PI)
		*theta -= 2.0 * M_PI;
	while (*theta < 0.0)
		*theta += 2.0 * M_PI;
}

void trendPlungeDegFromSpherical(double phi, double theta, double *trendDeg, double *plungeDeg) {
	*trendDeg = (90.0 - theta) / DEGREE;
	while (*trendDeg >= 360.0)
		*trendDeg -= 360.0;
	while (*trendDeg < 0.0)
		*trendDeg += 360.0;
	*plungeDeg = phi / DEGREE - 90.0;
}

void cartesianFromTrendPlungeDeg(double trendDeg, double plungeDeg, double *out) {
	double phi, theta;
	sphericalFromTrendPlungeDeg(trendDeg, plungeDeg, &phi, &theta);
	cartesianFromSpherical(1.0, phi, theta, out);
}

void trendPlungeDegFromCartesian(double *x, double *trendDeg, double *plungeDeg) {
	double rho, phi, theta;
	sphericalFromCartesian(x, &rho, &phi, &theta);
	trendPlungeDegFromSpherical(phi, theta, trendDeg, plungeDeg);
}

void cartesianFromHorizontal(double theta, double z, double *xOut, double *yOut, double *zOut) {
    double r;
    r = sqrt(1.0 - z * z);
  	*xOut = r * cos(theta);
  	*yOut = r * sin(theta);
  	*zOut = z;
}

void horizontalFromCartesian(double x, double y, double z, double *thetaOut, double *zOut) {
    *thetaOut = atan2(y, x);
    *zOut = z;
}

// Modified from https://www.securecoding.cert.org/confluence/display/c/MSC30-C.+Do+not+use+the+rand%28%29+function+for+generating+pseudorandom+numbers
void initializeRandom() {
    srand48(time(NULL));
    //struct timespec ts;
    //if (timespec_get(&ts, TIME_UTC) == 0)
    //    fprintf(stderr, "error: initializeRandom couldn't get time?\n");
    //srandom(ts.tv_nsec ^ ts.tv_sec);
    srandom(time(NULL));
}

// Before calling this function, ensure that initializeRandom has been called.
double doubleUniform(double a, double b) {
    return a + (b - a) * drand48();
}

// Before calling this function, ensure that initializeRandom has been called.
// Modified from https://en.wikipedia.org/wiki/Marsaglia_polar_method.
double doubleNormal(double mean, double stddev) {
	static int hasSpare = 0;
	static double spare;
	if (hasSpare) {
		hasSpare = 0;
		return mean + stddev * spare;
	}
	else {
		hasSpare = 1;
		static double u, v, s, root;
		do {
			u = doubleUniform(-1.0, 1.0);
			v = doubleUniform(-1.0, 1.0);
			s = u * u + v * v;
		} while (s >= 1.0 || s == 0.0);
		root = sqrt(-2.0 * log(s) / s);
		spare = v * root;
		return mean + stddev * u * root;
	}
}

// Before calling this function, ensure that initializeRandom has been called.
void unitUniform(double *out) {
    cartesianFromHorizontal(
      doubleUniform(-M_PI, M_PI),
      doubleUniform(-1.0, 1.0),
      &(out[0]), &(out[1]), &(out[2]));
}

// Helper function for unsignedUniform; do not call directly.
// From http://stackoverflow.com/questions/2509679/how-to-generate-a-random-number-from-within-a-range
// Assumes 0 <= max <= RAND_MAX; returns in the closed interval [0, max]
long random_at_most(long max) {
    unsigned long
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    num_bins = (unsigned long) max + 1,
    num_rand = (unsigned long) RAND_MAX + 1,
    bin_size = num_rand / num_bins,
    defect   = num_rand % num_bins;
    long x;
    do {
        // This POSIX random() function is apparently better than C's rand().
        x = random();
    }
    // This is carefully written not to overflow.
    while (num_rand - defect <= (unsigned long)x);
    // Truncated division is intentional.
    return x/bin_size;
}

// Returns unsigned integer drawn uniformly from the closed interval
// [a, b] = {a, a + 1, ..., b - 1, b}. Assumes a <= b.
// Before calling this function, ensure that initializeRandom has been called.
unsigned int unsignedUniform(unsigned int a, unsigned int b) {
    return a + random_at_most(b - a);
}

// From the integers {0, ..., n - 1}, chooses k at random.
// Assumes 0 <= k <= n.
// choices is an array of k integers, pre-allocated before this function is called.
// See Wikipedia: Reservoir sampling.
void subsetUniform(unsigned int n, unsigned int k, unsigned int *choices) {
  unsigned int i, j;
  // For starters, choose the first k numbers.
  for (i = 0; i < k; k += 1)
    choices[i] = i;
  // Give the others a chance to creep in.
  for (i = k; i < n; i += 1) {
    j = unsignedUniform(0, i - 1);
    if (j < k)
      choices[j] = i;
  }
}

// xs and ys are arrays of length n. xs is sorted in strictly increasing order.
// Regards y as a function of x, and guesses the value at the given x.
double doubleInterpolated(unsigned int n, double *xs, double *ys, double x) {
  // Clamp to the most extreme x-values.
  if (x <= xs[0])
    return xs[0];
  else if (x >= xs[n - 1])
    return xs[n - 1];
  else {
    // Binary search.
    int a, b, c;
    a = 0;
    b = n - 1;
    while (b - a >= 2) {
      c = (b + a) / 2;
      if (x <= xs[c])
        b = c;
      else
        a = c;
    }
    // Linearly interpolate.
    return ys[a] + (ys[b] - ys[a]) * (x - xs[a]) / (xs[b] - xs[a]);
  }
}


