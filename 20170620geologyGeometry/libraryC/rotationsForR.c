


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



// Compile with "R CMD SHLIB rotationsForR.c".

#include <R.h>
#include <Rdefines.h>
#include <stdio.h>

#include "linear.c"
#include "miscellany.c"
#include "rotations.c"
#include "rotationsmcmc.c"
#include "rotationsplot.c"



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
SEXP mcmcRotationWrappedTrivariateNormalC(SEXP rs, SEXP numTerms, SEXP numTuning, SEXP numBurnin,
    SEXP numCollection, SEXP numReport) {
  initializeRandom();
  // Massage the inputs into C.
  double *rots = doubleArrayFromSEXP(rs);
  int num = lengthFromSEXP(rs) / 9;
  int ter = intFromSEXP(numTerms);
  int tun = intFromSEXP(numTuning);
  int bur = intFromSEXP(numBurnin);
  int col = intFromSEXP(numCollection);
  int rep = intFromSEXP(numReport);
  if (col * tun < rep) {
    fprintf(stderr, "warning: mcmcRotationWrappedTrivariateNormalC: clamping numReported to numCollected * numTuning");
    rep = col * tun;
  }
  // Allocate memory: sBuffer, then etaBuffer, then three other scalar buffers.
  double *sBuffer, *etaBuffer, *aux;
  sBuffer = (double *)malloc(13 * col * tun * sizeof(double));
  if (sBuffer == NULL) {
    fprintf(stderr, "error: mcmcRotationWrappedTrivariateNormalC: could not allocate memory; try decreasing numCollection");
    return R_NilValue;
	} else {
    double meta[8];
    etaBuffer = &(sBuffer[9 * col * tun]);
    aux = &(sBuffer[10 * col * tun]);
    // Do the MCMC to fill meta, sBuffer, and etaBuffer.
    mcmcRotationWrappedTrivariateNormal(num, rots, ter, tun, bur, col, meta, sBuffer, etaBuffer);
    fprintf(stderr, "mcmcRotationWrappedTrivariateNormalC: beginning cleanup at ");
    printTime(stderr);
	  fprintf(stderr, "\n");
    SEXP out;
    int i;
    if (0.3 <= meta[1] && meta[1] <= 0.4 && 0.3 <= meta[3] && meta[3] <= 0.4) {
      // Compute mean, covariance, percentiles 0.900, 0.901, 0.902, ..., 0.999 of Mahalanobis distance.
      double center[9], error, var, covarInv[9], percentiles[101];
      unsigned int used;
      rotationMeanMany(col * tun, sBuffer, 10, 0.0000001, 1000, center, &error, &used, &var);
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

SEXP mcmcRotationWrappedSurgicalC(SEXP rs, SEXP numTuning, SEXP etaStep) {
  initializeRandom();
  // Massage the inputs into C.
  double *rots = doubleArrayFromSEXP(rs);
  int num = lengthFromSEXP(rs) / 9;
  int tun = intFromSEXP(numTuning);
  double eSt = doubleFromSEXP(etaStep);
  // Compute the acceptance rate.
  double accRate = mcmcRotationSurgical(num, rots, tun, eSt);
  // Return the result to R.
  SEXP result;
  PROTECT(result = NEW_NUMERIC(1));
  REAL(result)[0] = accRate;
  UNPROTECT(1);
  return result;
}

SEXP pvalueRotationBootstrappingC(SEXP rs, SEXP numBoots) {
  initializeRandom();
  // Massage the inputs into C.
  double *rots = doubleArrayFromSEXP(rs);
  int num = lengthFromSEXP(rs) / 9;
  int boo = intFromSEXP(numBoots);
  // Allocate memory.
  double *boots;
  boots = (double *)malloc(9 * boo * sizeof(double));
  if (boots == NULL) {
    fprintf(stderr, "error: pvalueRotationBootstrappingC: could not allocate memory; try decreasing numBoots");
    return R_NilValue;
  } else {
    // Call the underlying function.
    double center[9], covarInv[9], percentiles[101];
    int error = pvalueRotationBootstrapping(num, rots, boo, center, covarInv, percentiles, boots);
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

SEXP volumetricRawTrianglesFromSurface(int numBoxes, int *ns, double *poly) {
  // Count the triangles represented by ns.
  int b, d, t, v, numTris = 0;
  for (b = 0; b < numBoxes; b += 1)
    if (ns[b] >= 3)
      numTris += ns[b] - 2;
  // Expand the polygons into raw triangles.
  SEXP raw;
  PROTECT(raw = NEW_NUMERIC(9 * numTris));
  v = 0;
  d = 0;
  for (b = 0; b < numBoxes; b += 1)
    if (ns[b] >= 3) {
      for (t = 1; t <= ns[b] - 2; t += 1) {
        // Copy this box's vertices 0, t, t + 1 from poly into raw.
        REAL(raw)[0 * numTris + v + 0] = poly[d + 0 * 3 + 0];
        REAL(raw)[3 * numTris + v + 0] = poly[d + 0 * 3 + 1];
        REAL(raw)[6 * numTris + v + 0] = poly[d + 0 * 3 + 2];
        REAL(raw)[0 * numTris + v + 1] = poly[d + t * 3 + 0];
        REAL(raw)[3 * numTris + v + 1] = poly[d + t * 3 + 1];
        REAL(raw)[6 * numTris + v + 1] = poly[d + t * 3 + 2];
        REAL(raw)[0 * numTris + v + 2] = poly[d + (t + 1) * 3 + 0];
        REAL(raw)[3 * numTris + v + 2] = poly[d + (t + 1) * 3 + 1];
        REAL(raw)[6 * numTris + v + 2] = poly[d + (t + 1) * 3 + 2];
        v += 3;
      }
        d += ns[b] * 3;
    }
        UNPROTECT(1);
        return raw;
}
        
// rs is an array of 9 * n doubles representing n rotation matrices in column-major order.
// but these n rotations are actually n / groupSize orientations, if groupSize != 1.
// k and multiple are doubles. degree and depth are integers.
// Returns either NULL or a list of 9 * m numbers suitable for passing to plotBall as trianglesRaw.
SEXP volumetricKambTrianglesRawC(SEXP rs, SEXP multiple, SEXP k, SEXP degree, SEXP nonAdapt,
    SEXP adapt, SEXP groupSize) {
  // Massage the inputs into C.
  double *rots = doubleArrayFromSEXP(rs);
  int num = lengthFromSEXP(rs) / 9;
  double kay = doubleFromSEXP(k);
  double mul = doubleFromSEXP(multiple);
  int deg = intFromSEXP(degree);
  int non = intFromSEXP(nonAdapt);
  int ada = intFromSEXP(adapt);
  int gro = intFromSEXP(groupSize);
  // Allocate memory.
  double *aux, *out, *poly;
  int *ns;
  volumetricLevelMalloc(non, ada, &aux, &out, &ns, &poly);
  if (aux == NULL)
    return R_NilValue;
  else {
    int numBoxes = rotationKambSurface(num, rots, kay, mul, deg, non, ada, gro, aux, out, ns, poly);
    SEXP raw = volumetricRawTrianglesFromSurface(numBoxes, ns, poly);
    free(aux);
    return raw;
  }
}
        
// rs is an array of 9 * n doubles representing n rotation matrices in column-major order.
// k and multiple are doubles. degree and depth are integers.
// Returns either NULL or a list of 9 * m numbers suitable for passing to plotBall as trianglesRaw.
SEXP volumetricEllipsoidTrianglesRawC(SEXP center, SEXP covarInv, SEXP level, SEXP nonAdapt, SEXP adapt) {
  // Massage the inputs into C.
  double *cen = doubleArrayFromSEXP(center);
  double *inv = doubleArrayFromSEXP(covarInv);
  double lev = doubleFromSEXP(level);
  int non = intFromSEXP(nonAdapt);
  int ada = intFromSEXP(adapt);
  // Allocate memory.
  double *aux, *out, *poly;
  int *ns;
  volumetricLevelMalloc(non, ada, &aux, &out, &ns, &poly);
  if (aux == NULL)
    return R_NilValue;
  else {
    int numBoxes = rotationEllipsoidSurface(cen, inv, lev, non, ada, aux, out, ns, poly);
    SEXP raw = volumetricRawTrianglesFromSurface(numBoxes, ns, poly);
    free(aux);
    return raw;
  }
}


