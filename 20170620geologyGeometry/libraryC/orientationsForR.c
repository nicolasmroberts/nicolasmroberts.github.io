


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



// Compile with "R CMD SHLIB orientationsForR.c".

#include <R.h>
#include <Rdefines.h>
#include <stdio.h>

#include "linear.c"
#include "miscellany.c"
#include "rotations.c"
#include "rotationsmcmc.c"
#include "orientations.c"



// Warning: May have a side-effect on sexp?
int intFromSEXP(SEXP sexp) {
  SEXP temp = coerceVector(sexp, INTSXP);
  return INTEGER(temp)[0];
}

double doubleFromSEXP(SEXP sexp) {
  return REAL(sexp)[0];
}

double *doubleArrayFromSEXP(SEXP sexp) {
  return REAL(sexp);
}

int lengthFromSEXP(SEXP sexp) {
  return length(sexp);
}

// numCollection * numTuning should not exceed about 100,000,000 or we might overflow some 32-bit integers.
// If burn-in fails, then returns 4-dimensional double array nu, nuRateOld, logGamma, gammaRateOld.
// If burn-in succeeds, then returns big array of doubles in this format:
// 0-8: mean of MCMC means
// 9-17: inverse of covariance of MCMC means in tangent space at that mean
// 18-118: 101 percentiles (0.900, 0.901, 0.902, ..., 0.999, 1.000) of Mahalanobis norm of MCMC means
// 119-126: 8 meta-data: nu, nuRate, logGamma, gammaRate at end of burn-in and end of collection
// 127-(127 + numReport * 9 - 1): numReport MCMC means, evenly sampled from the MCMC means
// (127 + numReport * 9)-(127 + numReport * 9 + numReport - 1): the etas corresponding to the reported means
SEXP mcmcOrientationWrappedTrivariateNormalC(SEXP rs, SEXP group, SEXP numTerms, SEXP numTuning, SEXP numBurnin,
    SEXP numCollection, SEXP numReport) {
  initializeRandom();
  // Massage the inputs into C.
  double *rots = doubleArrayFromSEXP(rs);
  int num = lengthFromSEXP(rs) / 9;
  double *syms = doubleArrayFromSEXP(group);
  int gro = lengthFromSEXP(group) / 9;
  int ter = intFromSEXP(numTerms);
  int tun = intFromSEXP(numTuning);
  int bur = intFromSEXP(numBurnin);
  int col = intFromSEXP(numCollection);
  int rep = intFromSEXP(numReport);
  if (col * tun < rep) {
    fprintf(stderr, "warning: mcmcOrientationWrappedTrivariateNormalC: clamping numReported to numCollected * numTuning");
    rep = col * tun;
  }
  // Allocate memory: sBuffer, then etaBuffer, then three other scalar buffers.
  double *sBuffer, *etaBuffer, *aux;
  sBuffer = (double *)malloc(13 * col * tun * sizeof(double));
  if (sBuffer == NULL) {
    fprintf(stderr, "error: mcmcOrientationWrappedTrivariateNormalC: could not allocate memory; try decreasing numCollection");
    return R_NilValue;
  } else {
    double meta[8];
    etaBuffer = &(sBuffer[9 * col * tun]);
    aux = &(sBuffer[10 * col * tun]);
    // Do the MCMC to fill meta, sBuffer, and etaBuffer.
    mcmcOrientationWrappedTrivariateNormal(num, rots, gro, syms, ter, tun, bur, col, meta, sBuffer, etaBuffer);
    SEXP out;
    int i;
    if (0.3 <= meta[1] && meta[1] <= 0.4 && 0.3 <= meta[3] && meta[3] <= 0.4) {
      // Compute mean, covariance, percentiles 0.900, 0.901, 0.902, ..., 0.999 of Mahalanobis distance.
      double center[9], error, var, covarInv[9], percentiles[101];
      unsigned int used;
      fprintf(stderr, "mcmcOrientationWrappedTrivariateNormalC: beginning mean at "); printTime(stderr); fprintf(stderr, "\n");
      // This orientationMeanMany used to use 10 seeds, but now uses only 5, because it's so slow.
      orientationMeanMany(col * tun, sBuffer, 5, 0.0000001, 1000, center, &error, &used, &var, gro, syms);
      fprintf(stderr, "mcmcOrientationWrappedTrivariateNormalC: beginning cleanup at "); printTime(stderr); fprintf(stderr, "\n");
      for (i = 0; i < col * tun; i += 1)
        orientationNearest(&(sBuffer[9 * i]), center, gro, syms, &(sBuffer[9 * i]));
      pvalueRotationMahalanobis(col * tun, sBuffer, center, aux, covarInv, percentiles);
      // Start reporting results: center, covarInv, 101 percentiles, meta-data.
      PROTECT(out = NEW_NUMERIC(9 + 9 + 101 + 8 + rep * 9 + rep));
      matrixCopy(center, REAL(out));
      matrixCopy(covarInv, &(REAL(out)[9]));
      for (i = 0; i < 101; i += 1)
        REAL(out)[18 + i] = percentiles[i];
      for (i = 0; i < 8; i += 1)
        REAL(out)[119 + i] = meta[i];
      // Report some of sBuffer and etaBuffer. Be careful about integer overflow.
      if (1 <= rep) {
        int delta = col * tun / rep;
        for (i = 0; i < rep; i += 1) {
          matrixCopy(&(sBuffer[9 * (i * delta)]), &(REAL(out)[127 + 9 * i]));
          REAL(out)[127 + 9 * rep + i] = etaBuffer[i * delta];
        }
      }
    } else {
      // Burn-in failed, so just report its statistics.
      PROTECT(out = NEW_NUMERIC(4));
      for (i = 0; i < 4; i += 1)
        REAL(out)[i] = meta[i];
    }
    free(sBuffer);
    UNPROTECT(1);
    return out;
	}
}

SEXP pvalueOrientationBootstrappingC(SEXP rs, SEXP group, SEXP numBoots) {
  initializeRandom();
  // Massage the inputs into C.
  double *rots = doubleArrayFromSEXP(rs);
  int num = lengthFromSEXP(rs) / 9;
  double *syms = doubleArrayFromSEXP(group);
  int gro = lengthFromSEXP(group) / 9;
  int boo = intFromSEXP(numBoots);
  // Allocate memory.
  double *boots;
  boots = (double *)malloc(9 * boo * sizeof(double));
  if (boots == NULL) {
    fprintf(stderr, "error: pvalueOrientationBootstrappingC: could not allocate memory; try decreasing numBoots");
    return R_NilValue;
  } else {
    // Call the underlying function.
    double center[9], covarInv[9], percentiles[101];
    int error = pvalueOrientationBootstrapping(num, rots, gro, syms, boo, center, covarInv, percentiles, boots);
    if (error) {
      free(boots);
      return R_NilValue;
    } else {
      // Pack center, covarInv, percentiles, boots into one big vector.
      SEXP result;
      PROTECT(result = NEW_NUMERIC(9 + 9 + 101 + 9 * boo));
      matrixCopy(center, REAL(result));
      matrixCopy(covarInv, &(REAL(result)[9]));
      int i;
      for (i = 0; i < 101; i += 1)
        REAL(result)[18 + i] = percentiles[i];
      for (i = 0; i < boo; i += 1)
        matrixCopy(&(boots[9 * i]), &(REAL(result)[119 + 9 * i]));
      free(boots);
      UNPROTECT(1);
      return result;
    }
  }
}


