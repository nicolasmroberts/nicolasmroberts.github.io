


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



#include <stdio.h>



/// CONVERSIONS AMONG REPRESENTATIONS ///

void rotationFromAngleAxis(double angle, double *axis, double *out) {
    double U[9];
    U[0 +3* 0] = 0.0;
    U[0 +3* 1] = -axis[2];
    U[0 +3* 2] = axis[1];
    U[1 +3* 0] = -U[0 +3* 1];
    U[1 +3* 1] = 0.0;
    U[1 +3* 2] = -axis[0];
    U[2 +3* 0] = -U[0 +3* 2];
    U[2 +3* 1] = -U[1 +3* 2];
    U[2 +3* 2] = 0.0;
    // Q = I + (sin theta) U + (1 - cos theta) U^2.
    matrixMatrixMultiply(U, U, out);
    matrixScale(1.0 - cos(angle), out, out);
    matrixScale(sin(angle), U, U);
    matrixAdd(U, out, out);
    out[0 +3* 0] += 1.0;
    out[1 +3* 1] += 1.0;
    out[2 +3* 2] += 1.0;
}

void angleAxisFromRotation(double *r, double *angle, double *axis) {
	double cosine, square;
	cosine = 0.5 * (trace(r) - 1.0);
	if (cosine == 1.0) {
		*angle = 0.0;
		vectorSet(0.0, 0.0, -1.0, axis);
	}
	else {
		square = 1.0 + (r[0 +3* 0] - 1.0) / (1.0 - cosine);
		if (square <= 0.0)
			axis[0] = 0.0;
		else if (r[2 +3* 1] < r[1 +3* 2])
			axis[0] = -sqrt(square);
		else
			axis[0] = sqrt(square);
		square = 1.0 + (r[1 +3* 1] - 1.0) / (1.0 - cosine);
		if (square <= 0.0)
			axis[1] = 0.0;
		else if (r[0 +3* 2] < r[2 +3* 0])
			axis[1] = -sqrt(square);
		else
			axis[1] = sqrt(square);
		square = 1.0 + (r[2 +3* 2] - 1.0) / (1.0 - cosine);
		if (square <= 0.0)
			axis[2] = 0.0;
		else if (r[1 +3* 0] < r[0 +3* 1])
			axis[2] = -sqrt(square);
		else
			axis[2] = sqrt(square);
		vectorNormalize(axis, axis);
		*angle = arcCos(cosine);
	}
}

void volumetricFromAngleAxis(double angle, double *axis, double *vol) {
	double rho = pow((angle - sin(angle)) * 3.0 / (4.0 * M_PI * M_PI), 1.0 / 3.0);
	vectorScale(rho, axis, vol);
}

// Solve x - sin(x) == d on [0, pi]. Solution may not be in [0, pi] though.
// A reasonable tolerance is 0.00000001?
double xMinusSinXSolution(double d, double tolerance) {
	double x0, x1;
	// On [0, pi], x^3 / 9 is close to x - sin(x). Use that to choose seed.
	x0 = pow(9.0 * d, 1.0 / 3.0);
	// Now proceed by Newton's method.
	x1 = x0 - (x0 - sin(x0) - d) / (1.0 - cos(x0));
	while (fabs(x1 - x0) > tolerance) {
		x0 = x1;
		x1 = x0 - (x0 - sin(x0) - d) / (1.0 - cos(x0));
	}
	return x1;
}

void angleAxisFromVolumetric(double *vol, double *angle, double *axis) {
	double rho, d;
	rho = sqrt(scalarProduct(vol, vol));
	if (rho == 0.0) {
		*angle = 0.0;
		vectorSet(0.0, 0.0, -1.0, axis);
	} else {
		d = rho * rho * rho * 4.0 * M_PI * M_PI / 3.0;
		*angle = xMinusSinXSolution(d, 0.00000001);
		vectorScale(1.0 / rho, vol, axis);
	}
}

void volumetricFromRotation(double *rot, double *vol) {
	double angle, axis[3];
	angleAxisFromRotation(rot, &angle, axis);
	volumetricFromAngleAxis(angle, axis, vol);
}

void rotationFromVolumetric(double *vol, double *rot) {
	double angle, axis[3];
	angleAxisFromVolumetric(vol, &angle, axis);
	rotationFromAngleAxis(angle, axis, rot);
}

void rotationExp(double *w, double *out) {
	double alpha, wT[9], wTw[9], ww[9];
	// alpha = sqrt(tr(W^T W) / 2).
	matrixTranspose(w, wT);
	matrixMatrixMultiply(wT, w, wTw);
	alpha = sqrt(trace(wTw) * 0.5);
	if (alpha == 0.0)
		matrixIdentity(out);
	else {
		// exp(w) = I + (sin alpha) / alpha W + (1 - cos alpha) / alpha^2 W^2.
		matrixScale(sin(alpha) / alpha, w, out);
		matrixMatrixMultiply(w, w, ww);
		matrixScale((1 - cos(alpha)) / (alpha * alpha), ww, ww);
		matrixAdd(out, ww, out);
		out[0 +3* 0] += 1.0;
		out[1 +3* 1] += 1.0;
		out[2 +3* 2] += 1.0;
	}
}

void rotationLog(double *r, double *out) {
	double alpha, u[3];
	angleAxisFromRotation(r, &alpha, u);
	out[0 +3* 0] = 0.0;
	out[0 +3* 1] = alpha * -u[2];
	out[0 +3* 2] = alpha * u[1];
	out[1 +3* 0] = alpha * u[2];
	out[1 +3* 1] = 0.0;
	out[1 +3* 2] = alpha * -u[0];
	out[2 +3* 0] = alpha * -u[1];
	out[2 +3* 1] = alpha * u[0];
	out[2 +3* 2] = 0.0;
}

void rotationLogOLD(double *r, double *out) {
	double alpha, rT[9];
	alpha = arcCos((trace(r) - 1.0) * 0.5);
	if (alpha == 0.0)
		matrixZero(out);
	else {
		matrixTranspose(r, rT);
		matrixSubtract(r, rT, out);
		matrixScale(alpha / (2.0 * sin(alpha)), out, out);
	}
}

void antiFromVector(double *v, double *out) {
  out[0 +3* 0] = 0.0;
	out[0 +3* 1] = -v[2];
	out[0 +3* 2] = v[1];
	out[1 +3* 0] = v[2];
	out[1 +3* 1] = 0.0;
	out[1 +3* 2] = -v[0];
	out[2 +3* 0] = -v[1];
	out[2 +3* 1] = v[0];
	out[2 +3* 2] = 0.0;
}



/// GEODESIC METHODS ///

double rotationDistance(double *q, double *r) {
	double qT[9], qTr[9];
	matrixTranspose(q, qT);
	matrixMatrixMultiply(qT, r, qTr);
	// distance = arccos((tr Q^T R - 1) / 2).
	return arcCos((trace(qTr) - 1.0) * 0.5);
}

double rotationVariance(unsigned int n, double *rots, double *center) {
  double dist, var = 0.0;
  int i;
  for (i = 0; i < n; i += 1) {
    dist = rotationDistance(&(rots[9 * i]), center);
    var += dist * dist;
  }
  return var / (2.0 * n);
}

void rotationMean(unsigned int n, double *rots, double *seed,
		double epsilon, unsigned int iterations,
		double *mean, double *error, unsigned int *used) {
	double r[9], q[9], qT[9], qTr[9], logqTr[9], w[9], errorSquared;
	int i;
	// Q = seed, W = 0.
	matrixCopy(seed, q);
	matrixZero(w);
	*used = 0;
	do {
		*used += 1;
		// Q = Q exp(W).
		rotationExp(w, r);
		matrixMatrixMultiply(q, r, q);
		// W = (log(Q^T R1) + ... + log(Q^T Rn)) / n.
		matrixTranspose(q, qT);
		matrixMatrixMultiply(qT, rots, qTr);
		rotationLog(qTr, w);
		for (i = 1; i < n; i += 1) {
			matrixMatrixMultiply(qT, &(rots[9 * i]), qTr);
			rotationLog(qTr, logqTr);
			matrixAdd(w, logqTr, w);
		}
		matrixScale(1.0 / n, w, w);
		// errorSquared = |W|^2.
		errorSquared = frobenius(w, w);
	} while (errorSquared >= epsilon * epsilon && *used < iterations);
	*error = sqrt(errorSquared);
  matrixCopy(q, mean);
}

void rotationMeanMany(unsigned int n, double *rots, unsigned int numSeeds,
  	double epsilon, unsigned int iterations,
		double *mean, double *error, unsigned int *used, double *variance) {
  double bestMean[9], newMean[9], newError, newVariance;
  unsigned int newUsed, i, j;
  // All actual rotation variances are less than 5.
  *variance = 5.0;
  for (i = 0; i < numSeeds; i += 1) {
    j = unsignedUniform(0, n - 1);
    rotationMean(n, rots, &(rots[9 * j]), epsilon, iterations, newMean, &newError, &newUsed);
    newVariance = rotationVariance(n, rots, newMean);
    if (newVariance < *variance) {
      matrixCopy(newMean, bestMean);
      *error = newError;
      *used = newUsed;
      *variance = newVariance;
    }
  }
  matrixCopy(bestMean, mean);
}

void tangentVectorFromRotation(double *rot, double *center, double *vec) {
	double m[9], angle, axis[3];
	matrixTranspose(center, m);
	matrixMatrixMultiply(m, rot, m);
	angleAxisFromRotation(m, &angle, axis);
	vectorScale(angle, axis, vec);
}

void rotationCovariance(unsigned int n, double *rots, double *center, double *covar) {
	int i;
	double vec[3];
	matrixZero(covar);
	for (i = 0; i < n; i += 1) {
		tangentVectorFromRotation(&(rots[i * 9]), center, vec);
		covar[0 +3* 0] += vec[0] * vec[0];
		covar[0 +3* 1] += vec[0] * vec[1];
		covar[0 +3* 2] += vec[0] * vec[2];
		covar[1 +3* 1] += vec[1] * vec[1];
		covar[1 +3* 2] += vec[1] * vec[2];
		covar[2 +3* 2] += vec[2] * vec[2];
	}
	covar[1 +3* 0] = covar[0 +3* 1];
	covar[2 +3* 0] = covar[0 +3* 2];
	covar[2 +3* 1] = covar[1 +3* 2];
	matrixScale(1.0 / n, covar, covar);
}

int rotationPrincipalComponents(unsigned int n, double *rots, double *center, double *mags, double *dirs) {
    double covar[9];
    rotationCovariance(n, rots, center, covar);
    if (symmetricEigensystem(3, covar, dirs, mags) != 0)
        return 1;
    else {
        mags[0] = sqRoot(mags[0]);
        mags[1] = sqRoot(mags[1]);
        mags[2] = sqRoot(mags[2]);
        return 0;
    }
}

double rotationMahalanobisNorm(double *rot, double *center, double *covarInv) {
	double v[3], CIv[3];
	tangentVectorFromRotation(rot, center, v);
	vectorMatrixMultiply(covarInv, v, CIv);
	return sqrt(scalarProduct(v, CIv));
}

// Inputs: n, rots, center.
// aux is 3 * n doubles used as space for scratch work.
// Outputs: covarInv (3x3), percentiles (101) of Mahalanobis norm: 0.900, 0.901, 0.902, ..., 0.999, 1.000.
void pvalueRotationMahalanobis(unsigned int n, double *rots, double *center, double *aux,
    double *covarInv, double *percentiles) {
  int i;
  double covar[9], *norms, *sorted;
  norms = &(aux[n]);
  sorted = &(aux[2 * n]);
  rotationCovariance(n, rots, center, covar);
  symmetricInverse(covar, covarInv);
  for (i = 0; i < n; i += 1)
    norms[i] = rotationMahalanobisNorm(&(rots[9 * i]), center, covarInv);
  doubleMergeSort(n, norms, aux, sorted);
  for (i = 0; i < 100; i += 1)
    percentiles[i] = sorted[(int)floor((0.9 + i / 1000.0) * n)];
  percentiles[100] = sorted[n - 1];
}



/// SAMPLING FROM DISTRIBUTIONS ///

// Generates a random rotation matrix.
// Rotation is chosen uniformly from rotations with angle bound epsilon.
// To choose uniformly on all of SO(3), use epsilon = M_PI.
// Before calling this function, ensure that randomInitialize has been called.
// Based on rand_rotation by Jim Arvo (1991), Graphics Gems III.
// I think it's the subgroup algorithm of Diaconis and Shahshahani (1987).
void rotationUniform(double epsilon, double *out) {
	double x[3];
	x[0] = doubleUniform(0.0, epsilon / M_PI);
	x[1] = doubleUniform(0.0, 1.0);
	x[2] = doubleUniform(0.0, epsilon / M_PI);
    double theta = x[0] * 2.0 * M_PI; /* Rotation about the pole (Z).   */
    double phi   = x[1] * 2.0 * M_PI; /* Direction of pole deflection.  */
    double z     = x[2] * 2.0;      /* For magnitude of pole deflection.*/
    /* Compute a vector V used for distributing points over the sphere  */
    /* via the reflection I - V Transpose(V).  This formulation of V    */
    /* will guarantee that if x[1] and x[2] are uniformly distributed,  */
    /* the reflected points will be uniform on the sphere.  Note that V */
    /* has length sqrt(2) to eliminate the 2 in the Householder matrix. */
    double r  = sqrt(z);
    double Vx = sin(phi) * r;
    double Vy = cos(phi) * r;
    double Vz = sqrt(2.0 - z);
    /* Compute the row vector S = Transpose(V) * R, where R is a simple */
    /* rotation by theta about the z-axis.  No need to compute Sz since */
    /* it's just Vz.                                                    */
    double st = sin(theta);
    double ct = cos(theta);
    double Sx = Vx * ct - Vy * st;
    double Sy = Vx * st + Vy * ct;
    /* Construct the rotation matrix  (V Transpose(V) - I) R, which   */
    /* is equivalent to V S - R.                                        */
    out[0 +3* 0] = Vx * Sx - ct;
    out[0 +3* 1] = Vx * Sy - st;
    out[0 +3* 2] = Vx * Vz;
    out[1 +3* 0] = Vy * Sx + st;
    out[1 +3* 1] = Vy * Sy - ct;
    out[1 +3* 2] = Vy * Vz;
    out[2 +3* 0] = Vz * Sx;
    out[2 +3* 1] = Vz * Sy;
    out[2 +3* 2] = 1.0 - z;   /* This equals Vz * Vz - 1.0 */
}

void rotationWrappedTrivariateNormal(double *center, double kappa, double *out) {
  double v[3], w[9], r[9];
  double stddev = 1.0 / kappa;
  v[0] = doubleNormal(0, stddev);
  v[1] = doubleNormal(0, stddev);
  v[2] = doubleNormal(0, stddev);
  antiFromVector(v, w);
  rotationExp(w, r);
  matrixMatrixMultiply(center, r, out);
}

// Warning: The rejection sampling here is not efficient as kappa gets large.
// And it downright fails if kappa >= 26.7, due to overflow of exp(kappa^2) I think.
// !! Replace with sampling from matrix Fisher distribution?
void rotationIsotropicMatrixFisher(double *center, double kappa, double *out) {
	// Rejection-sample to get the angle alpha.
	double maximum, alpha, cosAlpha, f, b;
	if (kappa < 1.0 / sqrt(2.0))
		maximum = 2.0 * exp(-kappa * kappa);
	else
		maximum = exp(kappa * kappa - 1.0) / (kappa * kappa);
	do {
		alpha = doubleUniform(-M_PI, M_PI);
    cosAlpha = cos(alpha);
		f = doubleUniform(0.0, maximum);
    b = (1.0 - cosAlpha) * exp(kappa * kappa * cosAlpha);
	} while (f > b);
	// Rotate away from center through alpha about a random axis.
	double axis[3], r[9];
	unitUniform(axis);
	rotationFromAngleAxis(alpha, axis, r);
	matrixMatrixMultiply(center, r, out);
}



/// BOOTSTRAPPING ///

// Inputs:
// * n is the number of data.
// * rots is 9 * n doubles storing n rotation matrices (in column-major order).
// * numBoots is the number of bootstrap samples; 10,000 is typical.
// Outputs:
// * center is the mean of the bootstrapped means.
// * covarInv is the 3x3 inverse covariance matrix of the bootstrapped means about their center.
// * percentiles are 101 percentiles of Mahalanobis norm: 0.900, 0.901, 0.902, ..., 0.999, 1.000.
// * boots is 9 * numBoots doubles storing the bootstrapped means for later plotting.
// Returns:
// * 0 iff no error occurred (in allocating memory).
int pvalueRotationBootstrapping(unsigned int n, double *rots, unsigned int numBoots,
    double *center, double *covarInv, double *percentiles, double *boots) {
  unsigned int i, j, used;
  double error, variance, *aux;
  if (n > numBoots)
    aux = (double *)malloc(9 * n * sizeof(double));
  else
    aux = (double *)malloc(9 * numBoots * sizeof(double));
  if (aux == NULL) {
    fprintf(stderr, "error: pvalueRotationBootstrapping: malloc failed; try decreasing numBoots\n");
    return 1;
  } else {
    fprintf(stderr, "pvalueRotationBootstrapping: beginning bootstraps at "); printTime(stderr); fprintf(stderr, "\n");
    for (i = 0; i < numBoots; i += 1) {
      //fprintf(stderr, "%f ", i / (double)numBoots);
      for (j = 0; j < n; j += 1)
        matrixCopy(&(rots[9 * unsignedUniform(0, n - 1)]), &(aux[9 * j]));
      rotationMeanMany(n, aux, 10, 0.0000001, 1000, &(boots[9 * i]), &error, &used, &variance);
    }
    fprintf(stderr, "\npvalueRotationBootstrapping: beginning cleanup at "); printTime(stderr); fprintf(stderr, "\n");
    rotationMeanMany(numBoots, boots, 10, 0.0000001, 1000, center, &error, &used, &variance);
    pvalueRotationMahalanobis(numBoots, boots, center, aux, covarInv, percentiles);
    return 0;
  }
}


