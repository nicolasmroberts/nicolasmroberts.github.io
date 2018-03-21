


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



/// CONVERSIONS AMONG REPRESENTATIONS ///

// Given a 3x3 log-ellipsoid tensor, returns a 5d (if normalized) or 6D (if not)
// vector description of that tensor.
void ellVectorFromLog(double *logEll, int isNormalized, double *vec) {
  vec[0] = sqrt(2.0) * logEll[0 +3* 1];
  vec[1] = sqrt(2.0) * logEll[0 +3* 2];
  vec[2] = sqrt(2.0) * logEll[1 +3* 2];
  if (isNormalized) {
    vec[3] = (logEll[1 +3* 1] + logEll[0 +3* 0]) * sqrt(1.5);
    vec[4] = (logEll[1 +3* 1] - logEll[0 +3* 0]) / sqrt(2);
  } else {
    vec[3] = logEll[0 +3* 0];
    vec[4] = logEll[1 +3* 1];
    vec[5] = logEll[2 +3* 2];
  }
}

// Computes the matrix logarithm of a symmetric positive-definite matrix.
// Returns: 0 on success, nonzero on failure.
int ellipsoidLog(double *m, double *out) {
    double eigvecs[9], eigvals[3];
    symmetricEigensystem(3, m, eigvecs, eigvals);
    if (eigvals[0] <= 0.0 || eigvals[1] <= 0.0 || eigvals[2] <= 0.0)
        return 1;
    else {
        matrixSetDiagonal(log(eigvals[0]), log(eigvals[1]), log(eigvals[2]), out);
        matrixMatrixMultiply(eigvecs, out, out);
        matrixTranspose(out, out);
        matrixMatrixMultiply(eigvecs, out, out);
        return 0;
    }
}

// Computes the matrix exponential of a symmetric matrix.
void symmetricExp(double *m, double *out) {
    double eigvecs[9], eigvals[3];
    symmetricEigensystem(3, m, eigvecs, eigvals);
    matrixSetDiagonal(exp(eigvals[0]), exp(eigvals[1]), exp(eigvals[2]), out);
    matrixMatrixMultiply(eigvecs, out, out);
    matrixTranspose(out, out);
    matrixMatrixMultiply(eigvecs, out, out);
}



/// SHAPE ///

void semiaxesFromShape(double B12, double B13, double *a) {
    double c12, c13;
    c12 = (1.0 - B12) / (1.0 + B12);
    c13 = (1.0 - B13) / (1.0 + B13);
    a[1] = pow(c12 * c12 / c13, 1.0 / 6.0);
    a[2] = pow(c13 * c13 / c12, 1.0 / 6.0);
    a[0] = 1.0 / (a[1] * a[2]);
}

void shapeFromSemiaxes(double *a, double *B12, double *B13) {
    *B12 = (a[0] * a[0] - a[1] * a[1]) / (a[0] * a[0] + a[1] * a[1]);
    *B13 = (a[0] * a[0] - a[2] * a[2]) / (a[0] * a[0] + a[2] * a[2]);
}



/// STATISTICAL OPERATIONS ///

// Computes the geometric mean ellipsoid of a set of n ellipsoids.
void ellipsoidMean(int n, double *ells, double *out) {
    int k;
    double ellLog[9];
    double mean[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // Sum the logarithms.
    for (k = 0; k < n; k += 1) {
        ellipsoidLog(&(ells[9 * k]), ellLog);
        mean[0] += ellLog[0 +3* 0];
        mean[1] += ellLog[0 +3* 1];
        mean[2] += ellLog[0 +3* 2];
        mean[3] += ellLog[1 +3* 1];
        mean[4] += ellLog[1 +3* 2];
        mean[5] += ellLog[2 +3* 2];
    }
    // Divide the sum by the count to get the mean.
    mean[0] /= n;
    mean[1] /= n;
    mean[2] /= n;
    mean[3] /= n;
    mean[4] /= n;
    mean[5] /= n;
    // Exponentiate back to the ellipsoid.
    ellLog[0 +3* 0] = mean[0];
    ellLog[0 +3* 1] = mean[1];
    ellLog[0 +3* 2] = mean[2];
    ellLog[1 +3* 0] = mean[1];
    ellLog[1 +3* 1] = mean[3];
    ellLog[1 +3* 2] = mean[4];
    ellLog[2 +3* 0] = mean[2];
    ellLog[2 +3* 1] = mean[4];
    ellLog[2 +3* 2] = mean[5];
    symmetricExp(ellLog, out);
}

// Computes the variance of a set of n ellipsoids.
// Input:
// * Number n of ellipsoids.
// * The n ellipsoid tensors, packed into an array of 9 n doubles.
// * An array of 6 n doubles, used as storage by the function.
//   If NULL, then a dynamically allocated array is used instead.
// Output:
// * The variance, measured in the six-dimensional log space.
double ellipsoidSpread(int n, double *ells, double *storage) {
    int k, l;
    double ellLog[9];
    // Allocate temporary storage if needed.
    double *vecs;
    if (storage == NULL)
        vecs = malloc(6 * n * sizeof(double));
    else
        vecs = storage;
    if (vecs == NULL) {
        printf("error: ellipsoidSpread: malloc failed\n");
        return 0.0;
    }
    else {
        // Convert each ellipsoid into a vector (6 entries of the log).
        for (k = 0; k < n; k += 1) {
            ellipsoidLog(&(ells[9 * k]), ellLog);
            vecs[0 +6* k] = ellLog[0 +3* 0];
            vecs[1 +6* k] = ellLog[0 +3* 1];
            vecs[2 +6* k] = ellLog[0 +3* 2];
            vecs[3 +6* k] = ellLog[1 +3* 1];
            vecs[4 +6* k] = ellLog[1 +3* 2];
            vecs[5 +6* k] = ellLog[2 +3* 2];
        }
        // Compute the mean of the vectors.
        double mean[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (k = 0; k < n; k += 1)
            for (l = 0; l < 6; l += 1)
                mean[l] += vecs[l +6* k];
        for (l = 0; l < 6; l += 1)
            mean[l] /= n;
        // Shift the vectors by the mean.
        for (k = 0; k < n; k += 1)
            for (l = 0; l < 6; l += 1)
                vecs[l +6* k] -= mean[l];
        // Compute the trace of the sum of the outer products of the vectors.
        double trace = 0.0;
        for (k = 0; k < n; k += 1)
            for (l = 0; l < 6; l += 1)
                trace += vecs[l +6* k] * vecs[l +6* k];
        // Dividing by n - 1 gives an unbiased estimator.
        if (storage == NULL)
            free(vecs);
        return trace / (n - 1);
    }
}



/// RIGID SPHEROIDS AND ELLIPSOIDS ///

// Given ambient velocity gradient tensor L.
// And Bretherton shape factor describing the shape of a spheroid.
// Computes B D + W, which controls the Jeffery (1922) rigid dynamics.
// Output bDPlusW can safely alias input l.
void spheroidBDPlusW(double *l, double b, double *bDPlusW) {
  int i, j;
  double d[9], w[9];
  for (i = 0; i < 3; i += 1)
    for (j = 0; j < 3; j += 1) {
      d[i +3* j] = 0.5 * (l[i +3* j] + l[j +3* i]);
      w[i +3* j] = 0.5 * (l[i +3* j] - l[j +3* i]);
    }
  for (i = 0; i < 3; i += 1)
    for (j = 0; j < 3; j += 1)
      bDPlusW[i +3* j] = b * d[i +3* j] + w[i +3* j];
}

// Computes yDot = M y - (y^T M y) y, where M = B D + W.
void spheroidAxisVelocity(double *y, double *bDPlusW, double *yDot) {
    double my[3], ymyy[3];
    vectorMatrixMultiply(bDPlusW, y, my);
    vectorScale(scalarProduct(y, my), y, ymyy);
    vectorSubtract(my, ymyy, yDot);
}

// Simulates a rigid spheroid under the Jeffery (1922) dynamics.
// Input u0 can safely alias output u1.
void spheroidRigidRK4(
        double *u0,     // spheroid axis at t=0 (unit 3D Cartesian)
        double *BDPlusW,// B D + W, where B is shape factor, D symmetric part of VGT, W antisymmetric part
        int n,          // number of steps to use in approximation
        double *u1) {   // output spheroid axis at t=1 (unit 3D Cartesian)
    int i;
    // Use RK4 on uDot = M u - (ut M u) u, where y is u and M is B D + W.
    double y[3], h, v[3], y2[3], y3[3], y4[3];
    double k1[3], k2[3], k3[3], k4[3];
    vectorCopy(u0, y);
    h = 1.0 / n;
    for (i = 0; i < n; i += 1) {
        // k1 = f(y).
        spheroidAxisVelocity(y, BDPlusW, k1);
        // k2 = f(y + (h / 2) k1).
        vectorScale(0.5 * h, k1, v);
        vectorAdd(y, v, y2);
        vectorNormalize(y2, y2);
        spheroidAxisVelocity(y2, BDPlusW, k2);
        // k3 = f(y + (h / 2) k2).
        vectorScale(0.5 * h, k2, v);
        vectorAdd(y, v, y3);
        vectorNormalize(y3, y3);
        spheroidAxisVelocity(y3, BDPlusW, k3);
        // k4 = f(y + h k3).
        vectorScale(h, k3, v);
        vectorAdd(y, v, y4);
        vectorNormalize(y4, y4);
        spheroidAxisVelocity(y4, BDPlusW, k4);
        // y = y + (h / 6) (k1 + 2 k2 + 2 k3 + k4).
        vectorScale(2.0, k2, k2);
        vectorScale(2.0, k3, k3);
        vectorAdd(k1, k2, v);
        vectorAdd(v, k3, v);
        vectorAdd(v, k4, v);
        vectorScale(h / 6.0, v, v);
        vectorAdd(y, v, y);
        vectorNormalize(y, y);
    }
    vectorCopy(y, u1);
}

// Computes Qdot = Wtilde Q - Q W.
void triaxialOrientationVelocity(double *Q, double *D, double *W, double *b, double *Qdot) {
    double Dtilde[9], Wtilde[9], WtildeQ[9], QW[9];
    // Compute Dtilde = Q D Qt.
    matrixMatrixMultiply(Q, D, Dtilde);
    matrixTranspose(Dtilde, Dtilde);
    matrixMatrixMultiply(Q, Dtilde, Dtilde);
    // Compute Wtilde.
    Wtilde[0 +3* 0] = 0.0;
    Wtilde[0 +3* 1] = b[0] * Dtilde[0 +3* 1];
    Wtilde[0 +3* 2] = -b[2] * Dtilde[0 +3* 2];
    Wtilde[1 +3* 0] = -Wtilde[0 +3* 1];
    Wtilde[1 +3* 1] = 0.0;
    Wtilde[1 +3* 2] = b[1] * Dtilde[1 +3* 2];
    Wtilde[2 +3* 0] = -Wtilde[0 +3* 2];
    Wtilde[2 +3* 1] = -Wtilde[1 +3* 2];
    Wtilde[2 +3* 2] = 0.0;
    // Compute Qdot = Wtilde Q - Q W.
    matrixMatrixMultiply(Wtilde, Q, WtildeQ);
    matrixMatrixMultiply(Q, W, QW);
    matrixSubtract(WtildeQ, QW, Qdot);
    //printf("B12 = %f, B13 = %f, Wtilde = \n", b[0], -b[2]);
    //printMatrix(Wtilde);
}

// Simulates a rigid triaxial ellipsoid under the Jeffery (1922) dynamics.
void triaxialRigidRK4(
        double *Q0,   // ellipsoid orientation tensor at t=0
        double *a,    // three ellipsoid semi-axes
        double *L,    // ambient velocity gradient tensor
        int n,        // number of steps to use in approximation
        double *Q1) { // output ellipsoid tensor at t=1
    double b[3], D[9], W[9];
    int i, j;
    // Encode shape factors in b[0] = B12, b[1] = B23, b[2] = B31.
    b[0] = (a[0] * a[0] - a[1] * a[1]) / (a[0] * a[0] + a[1] * a[1]);
    b[1] = (a[1] * a[1] - a[2] * a[2]) / (a[1] * a[1] + a[2] * a[2]);
    b[2] = (a[2] * a[2] - a[0] * a[0]) / (a[2] * a[2] + a[0] * a[0]);
    // D and W are the symmetric and antisymmetric parts of L.
    for (i = 0; i < 3; i += 1)
        for (j = 0; j < 3; j += 1) {
            D[i +3* j] = 0.5 * (L[i +3* j] + L[j +3* i]);
            W[i +3* j] = 0.5 * (L[i +3* j] - L[j +3* i]);
        }
    // Perform Runge-Kutta on Qdot = Wtilde Q - Q W, where y is Q.
    double y[9], h, v[9], y2[9], y3[9], y4[9];
    double k1[9], k2[9], k3[9], k4[9];
    matrixCopy(Q0, y);
    h = 1.0 / n;
    for (i = 0; i < n; i += 1)
    {
        // k1 = f(y).
        triaxialOrientationVelocity(y, D, W, b, k1);
        // k2 = f(y + (h / 2) k1).
        matrixScale(0.5 * h, k1, v);
        matrixAdd(y, v, y2);
        triaxialOrientationVelocity(y2, D, W, b, k2);
        // k3 = f(y + (h / 2) k2).
        matrixScale(0.5 * h, k2, v);
        matrixAdd(y, v, y3);
        triaxialOrientationVelocity(y3, D, W, b, k3);
        // k4 = f(y + h k3).
        matrixScale(h, k3, v);
        matrixAdd(y, v, y4);
        triaxialOrientationVelocity(y4, D, W, b, k4);
        // y = y + (h / 6) (k1 + 2 k2 + 2 k3 + k4).
        matrixScale(2.0, k2, k2);
        matrixScale(2.0, k3, k3);
        matrixAdd(k1, k2, v);
        matrixAdd(v, k3, v);
        matrixAdd(v, k4, v);
        matrixScale(h / 6.0, v, v);
        matrixAdd(y, v, y);
    }
    matrixCopy(y, Q1);
}

// Simulates a rigid ellipsoid (of any kind) under the Jeffery (1922) dynamics.
void ellipsoidRigidRK4(
        double *Q0,   // ellipsoid orientation tensor at t=0
        double *a,    // three ellipsoid semi-axes
        double *L,    // ambient velocity gradient tensor
        int n,        // number of steps to use in approximation
        double *Q1) { // output ellipsoid tensor at t=1
    double u[3], B, D[9], W[9], BDPlusW[9];
    int i, j;
    if (a[0] == a[1] && a[1] == a[2])
        // Handle the sphere case.
        matrixCopy(Q0, Q1);
    else if (a[0] == a[1]) {
        // Handle a spheroid case.
        B = (a[2] * a[2] - a[0] * a[0]) / (a[2] * a[2] + a[0] * a[0]);
        for (i = 0; i < 3; i += 1)
          for (j = 0; j < 3; j += 1) {
            D[i +3* j] = 0.5 * (L[i +3* j] + L[j +3* i]);
            W[i +3* j] = 0.5 * (L[i +3* j] - L[j +3* i]);
            BDPlusW[i +3* j] = B * D[i +3* j] + W[i +3* j];
          }
        vectorSet(Q0[2], Q0[5], Q0[8], u);
        spheroidRigidRK4(u, BDPlusW, n, &(Q1[6]));
        vectorsOrthogonalToVector(&(Q1[6]), &(Q1[0]), &(Q1[3]));
        matrixTranspose(Q1, Q1);
    }
    else if (a[1] == a[2]) {
        // Handle a spheroid case.
        B = (a[0] * a[0] - a[1] * a[1]) / (a[0] * a[0] + a[1] * a[1]);
        for (i = 0; i < 3; i += 1)
          for (j = 0; j < 3; j += 1) {
            D[i +3* j] = 0.5 * (L[i +3* j] + L[j +3* i]);
            W[i +3* j] = 0.5 * (L[i +3* j] - L[j +3* i]);
            BDPlusW[i +3* j] = B * D[i +3* j] + W[i +3* j];
          }
        vectorSet(Q0[0], Q0[3], Q0[6], u);
        spheroidRigidRK4(u, BDPlusW, n, &(Q1[0]));
        vectorsOrthogonalToVector(&(Q1[0]), &(Q1[3]), &(Q1[6]));
        matrixTranspose(Q1, Q1);
    }
    else if (a[2] == a[0]) {
        // Handle a spheroid case.
        B = (a[1] * a[1] - a[2] * a[2]) / (a[1] * a[1] + a[2] * a[2]);
        for (i = 0; i < 3; i += 1)
          for (j = 0; j < 3; j += 1) {
            D[i +3* j] = 0.5 * (L[i +3* j] + L[j +3* i]);
            W[i +3* j] = 0.5 * (L[i +3* j] - L[j +3* i]);
            BDPlusW[i +3* j] = B * D[i +3* j] + W[i +3* j];
          }
        vectorSet(Q0[1], Q0[4], Q0[7], u);
        spheroidRigidRK4(u, BDPlusW, n, &(Q1[3]));
        vectorsOrthogonalToVector(&(Q1[3]), &(Q1[6]), &(Q1[0]));
        matrixTranspose(Q1, Q1);
    }
    else 
        // Handle the triaxial case.
        triaxialRigidRK4(Q0, a, L, n, Q1);
}


