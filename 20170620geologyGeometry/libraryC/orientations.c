


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



/// GEODESIC METHODS ///

double orientationDistance(double *q, double *r, unsigned int g, double *group) {
  double rT[9], qRT[9], jQRT[9], newTrace, maxTrace = -2.0;
  int i;
	matrixTranspose(r, rT);
	matrixMatrixMultiply(q, rT, qRT);
  // Maximize tr(J Q R^T) over J in G.
  for (i = 0; i < g; i += 1) {
    matrixMatrixMultiply(&(group[i * 9]), qRT, jQRT);
    newTrace = trace(jQRT);
    if (newTrace > maxTrace)
      maxTrace = newTrace;
  }
	// Distance = arccos((tr J Q R^T - 1) / 2).
	return arcCos((maxTrace - 1.0) * 0.5);
}

// Input and output may safely alias each other.
void orientationNearest(double *q, double *center, unsigned int g, double *group, double *out) {
  double cT[9], qCT[9], jQCT[9], newTrace, trJQCTBest = -2.0;
  int i, iBest = -1;
  matrixTranspose(center, cT);
	matrixMatrixMultiply(q, cT, qCT);
  // Maximize tr(J Q C^T) over J in G.
  for (i = 0; i < g; i += 1) {
    matrixMatrixMultiply(&(group[i * 9]), qCT, jQCT);
    newTrace = trace(jQCT);
    if (newTrace > trJQCTBest) {
      trJQCTBest = newTrace;
      iBest = i;
    }
  }
  matrixMatrixMultiply(&(group[iBest * 9]), q, out);
}

double orientationVariance(unsigned int n, double *rots, double *center, unsigned int groupSize, double *group) {
  double dist, var = 0.0;
  int i;
  for (i = 0; i < n; i += 1) {
    dist = orientationDistance(&(rots[9 * i]), center, groupSize, group);
    var += dist * dist;
  }
  return var / (2.0 * n);
}

void orientationMean(unsigned int n, double *rots, double *seed,
  	double epsilon, unsigned int iterations,
		double *mean, double *error, unsigned int *used,
    unsigned int groupSize, double *group) {
	double r[9], q[9], qT[9], qTR[9], logQTR[9], w[9], changeSquared;
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
    matrixZero(w);
    matrixTranspose(q, qT);
		for (i = 0; i < n; i += 1) {
      orientationNearest(&(rots[9 * i]), q, groupSize, group, r);
			matrixMatrixMultiply(qT, r, qTR);
			rotationLog(qTR, logQTR);
			matrixAdd(w, logQTR, w);
		}
		matrixScale(1.0 / n, w, w);
		changeSquared = frobenius(w, w);
	} while (changeSquared >= epsilon * epsilon && *used < iterations);
	*error = sqrt(changeSquared);
  matrixCopy(q, mean);
}

// Before calling this function, ensure that initializeRandom has been called.
// calls orientationMeanMany somewhere between numSeeds and n times, and picks best answer.
void orientationMeanMany(unsigned int n, double *rots, unsigned int numSeeds,
    double epsilon, unsigned int iterations,
  	double *mean, double *error, unsigned int *used, double *variance,
    unsigned int groupSize, double *group) {
  unsigned int i, j, newUsed;
  double newVar, newMean[9], newError, bestMean[9];
  // 5 exceeds all variances that will ever happen.
  *variance = 5.0;
  if (n / numSeeds < 2)
    numSeeds = n;
  for (i = 0; i < numSeeds; i += 1) {
    if (numSeeds == n)
      j = i;
    else
      j = unsignedUniform(0, n - 1);
    orientationMean(n, rots, &(rots[j * 9]), epsilon, iterations, newMean, &newError, &newUsed, groupSize, group);
    newVar = orientationVariance(n, rots, newMean, groupSize, group);
    if (newVar < *variance) {
      *variance = newVar;
      matrixCopy(newMean, bestMean);
      *error = newError;
      *used = newUsed;
    }
  }
  matrixCopy(bestMean, mean);
}

// Inputs:
// * n is the number of data.
// * rots is 9 * n doubles storing n rotation matrices (in column-major order).
// * g is the size of the group.
// * group is 9 * g doubles storing g rotation matrices.
// * numBoots is the number of bootstrap samples; 10,000 is typical.
// Outputs:
// * center is the mean of the bootstrapped means.
// * covarInv is the 3x3 inverse covariance matrix of the bootstrapped means about their center.
// * percentiles are 101 percentiles of Mahalanobis norm: 0.900, 0.901, 0.902, ..., 0.999, 1.000.
// * boots is 9 * numBoots doubles storing the bootstrapped means for later plotting.
// Returns:
// * 0 iff no error occurred (in allocating memory).
int pvalueOrientationBootstrapping(unsigned int n, double *rots, unsigned int g, double *group, unsigned int numBoots,
    double *center, double *covarInv, double *percentiles, double *boots) {
  unsigned int i, j, used;
  double error, variance, *aux;
  if (n > numBoots)
    aux = (double *)malloc(9 * n * sizeof(double));
  else
    aux = (double *)malloc(9 * numBoots * sizeof(double));
  if (aux == NULL) {
    fprintf(stderr, "error: pvalueOrientationBootstrapping: malloc failed; try decreasing numBoots\n");
    return 1;
  } else {
    fprintf(stderr, "pvalueOrientationBootstrapping: beginning bootstraps at "); printTime(stderr); fprintf(stderr, "\n");
    for (i = 0; i < numBoots; i += 1) {
      //fprintf(stderr, "%f ", i / (double)numBoots);
      for (j = 0; j < n; j += 1)
        matrixCopy(&(rots[9 * unsignedUniform(0, n - 1)]), &(aux[9 * j]));
      // This orientationMeanMany used to use 10 seeds, but now uses only 5, because it's so slow.
      orientationMeanMany(n, aux, 5, 0.0000001, 1000, &(boots[9 * i]), &error, &used, &variance, g, group);
    }
    fprintf(stderr, "\npvalueOrientationBootstrapping: beginning cleanup at "); printTime(stderr); fprintf(stderr, "\n");
    orientationMeanMany(numBoots, boots, 10, 0.0000001, 1000, center, &error, &used, &variance, g, group);
    for (i = 0; i < numBoots; i += 1)
      orientationNearest(&(boots[9 * i]), center, g, group, &(boots[9 * i]));
    pvalueRotationMahalanobis(numBoots, boots, center, aux, covarInv, percentiles);
    free(aux);
    return 0;
  }
}



/// MARKOV CHAIN MONTE CARLO ///

// I'm leaving out the 1/g factor, because it will cancel itself anyway.
double likelihoodOrientationWrappedTrivariateNormal(unsigned int n, double *data, double *s, double kappa, 
    unsigned int terms, unsigned int g, double *group) {
  double sT[9], rST[9], rSTJ[9], alpha, tr, factor, like = 1.0;
  int i, j;
  matrixTranspose(s, sT);
  for (i = 0; i < n; i += 1) {
    factor = 0.0;
    matrixMatrixMultiply(&(data[9 * i]), sT, rST);
    for (j = 0; j < g; j += 1) {
      // alpha == rotationDistance(s, j * datum) and 3 - tr == 2 - 2 cos(alpha).
      matrixMatrixMultiply(rST, &(group[9 * j]), rSTJ);
      tr = trace(rSTJ);
      alpha = arcCos((tr - 1.0) * 0.5);
      factor += wtndAngularDensity(alpha, kappa, terms) / (3.0 - tr);
    }
    like *= factor;
  }
  return like;
}

// Returns 1 if challenger won, and 0 if incumbent won.
// Input s can safely alias output out.
int sOrientationWrappedTrivariateNormal(unsigned int n, double *data, double *s, double kappa, 
    unsigned int terms, double nu, unsigned int g, double *group, double *out) {
  double numer, denom, r, sNew[9];
  rotationWrappedTrivariateNormal(s, exp(-nu), sNew);
  // The J factor cancels in the ratio, so just use the likelihood part.
  numer = likelihoodOrientationWrappedTrivariateNormal(n, data, sNew, kappa, terms, g, group);
  denom = likelihoodOrientationWrappedTrivariateNormal(n, data, s, kappa, terms, g, group);
  r = numer / denom;
  if (r >= 1.0 || doubleUniform(0.0, 1.0) <= r) {
    matrixCopy(sNew, out);
    return 1;
  } else {
    matrixCopy(s, out);
    return 0;
  }
}

// Returns 1 if challenger won, and 0 if incumbent won.
int etaOrientationWrappedTrivariateNormal(unsigned int n, double *data, double *s, double eta, unsigned int terms,
    double logGamma, unsigned int g, double *group, double *out) {
  double numer, denom, r, etaNew;
  etaNew = doubleNormal(eta, exp(logGamma));
  numer = likelihoodOrientationWrappedTrivariateNormal(n, data, s, exp(-etaNew), terms, g, group) * wtndJ(etaNew);
  denom = likelihoodOrientationWrappedTrivariateNormal(n, data, s, exp(-eta), terms, g, group) * wtndJ(eta);
  r = numer / denom;
  if (r >= 1.0 || doubleUniform(0.0, 1.0) <= r) {
    *out = etaNew;
    return 1;
  } else {
    *out = eta;
    return 0;
  }
}

// data is an array of 9 * n doubles, to store the data rotations in column-major order.
// group is 9 * g doubles storing g rotation matrices.
// terms controls certain asymptotic expansions; 10 seems like a good value.
// tuning is the number of MCMC iterations between tunings; 1,000 or 10,000 seems good.
// burnin bounds the number of tunings undertaken in the burn-in phase;
// 100 or 1,000 seems good, for a total of about 10^6 MCMC iterations.
// collection is the number of tunings in the collection phase;
// 10^7 / tuning is common, for a total of 10^7 MCMC iterations.
// meta is an array of 8 doubles reporting nu, nuRate, logGamma, and gammaRate
// at the end of burn-in and the end of collection.
// sBuffer is an array of 9 * tuning * collection doubles, to store the Ss in column-major order.
// etaBuffer is an array of tuning * collection doubles, to store the etas.
// Before calling this function, ensure that initializeRandom has been called.
void mcmcOrientationWrappedTrivariateNormal(unsigned int n, double *data, unsigned int g, double *group,
    unsigned int terms, unsigned int tuning, unsigned int burnin, unsigned int collection,
    double *meta, double *sBuffer, double *etaBuffer) {
  double error, var, center[9], s[9], eta, nu, gamma, logGamma;
  int i, j, stable = 0, nuCount = 0, gammaCount = 0;
  unsigned int used;
  fprintf(stderr, "mcmcOrientationWrappedTrivariateNormal: beginning seeding at "); printTime(stderr); fprintf(stderr, "\n");
  orientationMeanMany(n, data, 10, 0.0000001, 1000, center, &error, &used, &var, g, group);
  for (i = 0; i < n; i += 1)
    orientationNearest(&(data[9 * i]), center, g, group, &(data[9 * i]));
  mcmcRotationSeeding(n, data, tuning, s, &eta, &nu, &gamma);
  logGamma = log(gamma);
  fprintf(stderr, "mcmcOrientationWrappedTrivariateNormal: beginning burn-in at "); printTime(stderr); fprintf(stderr, "\n");
  double nuRateOld = -1.0, nuOld, gammaRateOld = -1.0, logGammaOld;
  i = 0;
  while (i < burnin && stable < 3) {
    i += 1;
    for (j = 0; j < tuning; j += 1) {
      nuCount += sOrientationWrappedTrivariateNormal(n, data, s, exp(-eta), terms, nu, g, group, s);
      gammaCount += etaOrientationWrappedTrivariateNormal(n, data, s, eta, terms, logGamma, g, group, &eta);
    }
    wtndTuning(&nuRateOld, &nuOld, tuning, &nuCount, &nu);
    wtndTuning(&gammaRateOld, &logGammaOld, tuning, &gammaCount, &logGamma);
    if (nu == nuOld && logGamma == logGammaOld)
      stable += 1;
    else
      stable = 0;
    fprintf(stderr, "tuning %d, last nu-rate %f, last gamma-rate %f, nu %f, logGamma %f\n",
      i, nuRateOld, gammaRateOld, nu, logGamma);
  }
  // Store the first chunk of metadata.
  meta[0] = nu;
  meta[1] = nuRateOld;
  meta[2] = logGamma;
  meta[3] = gammaRateOld;
  if (0.3 <= nuRateOld && nuRateOld <= 0.4 && 0.3 <= gammaRateOld && gammaRateOld <= 0.4) {
    fprintf(stderr, "mcmcOrientationWrappedTrivariateNormal: beginning collection at "); printTime(stderr); fprintf(stderr, "\n");
    for (i = 0; i < collection; i += 1) {
      for (j = 0; j < tuning; j += 1) {
        nuCount += sOrientationWrappedTrivariateNormal(n, data, s, exp(-eta), terms, nu, g, group, s);//!!!
        gammaCount += etaOrientationWrappedTrivariateNormal(n, data, s, eta, terms, logGamma, g, group, &eta);//!!!
        matrixCopy(s, &(sBuffer[9 * (i * tuning + j)]));
        etaBuffer[i * tuning + j] = eta;
      }
      wtndTuning(&nuRateOld, &nuOld, tuning, &nuCount, &nu);
      wtndTuning(&gammaRateOld, &logGammaOld, tuning, &gammaCount, &logGamma);
      fprintf(stderr, "tuning %d, last nu-rate %f, last gamma-rate %f, nu %f, logGamma %f\n",
        i, nuRateOld, gammaRateOld, nu, logGamma);
    }
    // Store the second chunk of metadata.
    meta[4] = nu;
    meta[5] = nuRateOld;
    meta[6] = logGamma;
    meta[7] = gammaRateOld;
  } else
    fprintf(stderr, "mcmcOrientationWrappedTrivariateNormal: aborting after failed burn-in\n");
}


