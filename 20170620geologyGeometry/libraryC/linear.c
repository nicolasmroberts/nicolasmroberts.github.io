


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



// various routines for matrices, mostly 3x3
// matrices are stored in column-major order as arrays of doubles
// so the (i, j)th entry of matrix M is m[i +3* j],
// where i, j = 0, 1, 2, not 1, 2, 3

#include <stdio.h>
#include <float.h>                           // required for DBL_EPSILON
#include <math.h>                            // required for fabs()



// WARNING: If you store matrices as I do --- flat double arrays, in 
// column-major order --- then this function stores the eigenvectors as the 
// ROWS of the matrix eigenvectors, not the columns!
// from http://www.mymathlib.com/matrices/eigen/symmetric.html
////////////////////////////////////////////////////////////////////////////////
// File: jacobi_cyclic_method.c                                               //
// Routines:                                                                  //
//    Jacobi_Cyclic_Method                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Jacobi_Cyclic_Method                                                 //
//            (double eigenvalues[], double *eigenvectors, double *A, int n)  //
//                                                                            //
//  Description:                                                              //
//     Find the eigenvalues and eigenvectors of a symmetric n x n matrix A    //
//     using the Jacobi method. Upon return, the input matrix A will have     //
//     been modified.                                                         //
//     The Jacobi procedure for finding the eigenvalues and eigenvectors of a //
//     symmetric matrix A is based on finding a similarity transformation     //
//     which diagonalizes A.  The similarity transformation is given by a     //
//     product of a sequence of orthogonal (rotation) matrices each of which  //
//     annihilates an off-diagonal element and its transpose.  The rotation   //
//     effects only the rows and columns containing the off-diagonal element  //
//     and its transpose, i.e. if a[i][j] is an off-diagonal element, then    //
//     the orthogonal transformation rotates rows a[i][] and a[j][], and      //
//     equivalently it rotates columns a[][i] and a[][j], so that a[i][j] = 0 //
//     and a[j][i] = 0.                                                       //
//     The cyclic Jacobi method considers the off-diagonal elements in the    //
//     following order: (0,1),(0,2),...,(0,n-1),(1,2),...,(n-2,n-1).  If the  //
//     the magnitude of the off-diagonal element is greater than a treshold,  //
//     then a rotation is performed to annihilate that off-diagnonal element. //
//     The process described above is called a sweep.  After a sweep has been //
//     completed, the threshold is lowered and another sweep is performed     //
//     with the new threshold. This process is completed until the final      //
//     sweep is performed with the final threshold.                           //
//     The orthogonal transformation which annihilates the matrix element     //
//     a[k][m], k != m, is Q = q[i][j], where q[i][j] = 0 if i != j, i,j != k //
//     i,j != m and q[i][j] = 1 if i = j, i,j != k, i,j != m, q[k][k] =       //
//     q[m][m] = cos(phi), q[k][m] = -sin(phi), and q[m][k] = sin(phi), where //
//     the angle phi is determined by requiring a[k][m] -> 0.  This condition //
//     on the angle phi is equivalent to                                      //
//               cot(2 phi) = 0.5 * (a[k][k] - a[m][m]) / a[k][m]             //
//     Since tan(2 phi) = 2 tan(phi) / (1.0 - tan(phi)^2),                    //
//               tan(phi)^2 + 2cot(2 phi) * tan(phi) - 1 = 0.                 //
//     Solving for tan(phi), choosing the solution with smallest magnitude,   //
//       tan(phi) = - cot(2 phi) + sgn(cot(2 phi)) sqrt(cot(2phi)^2 + 1).     //
//     Then cos(phi)^2 = 1 / (1 + tan(phi)^2) and sin(phi)^2 = 1 - cos(phi)^2 //
//     Finally by taking the sqrts and assigning the sign to the sin the same //
//     as that of the tan, the orthogonal transformation Q is determined.     //
//     Let A" be the matrix obtained from the matrix A by applying the        //
//     similarity transformation Q, since Q is orthogonal, A" = Q'AQ, where Q'//
//     is the transpose of Q (which is the same as the inverse of Q).  Then   //
//         a"[i][j] = Q'[i][p] a[p][q] Q[q][j] = Q[p][i] a[p][q] Q[q][j],     //
//     where repeated indices are summed over.                                //
//     If i is not equal to either k or m, then Q[i][j] is the Kronecker      //
//     delta.   So if both i and j are not equal to either k or m,            //
//                                a"[i][j] = a[i][j].                         //
//     If i = k, j = k,                                                       //
//        a"[k][k] =                                                          //
//           a[k][k]*cos(phi)^2 + a[k][m]*sin(2 phi) + a[m][m]*sin(phi)^2     //
//     If i = k, j = m,                                                       //
//        a"[k][m] = a"[m][k] = 0 =                                           //
//           a[k][m]*cos(2 phi) + 0.5 * (a[m][m] - a[k][k])*sin(2 phi)        //
//     If i = k, j != k or m,                                                 //
//        a"[k][j] = a"[j][k] = a[k][j] * cos(phi) + a[m][j] * sin(phi)       //
//     If i = m, j = k, a"[m][k] = 0                                          //
//     If i = m, j = m,                                                       //
//        a"[m][m] =                                                          //
//           a[m][m]*cos(phi)^2 - a[k][m]*sin(2 phi) + a[k][k]*sin(phi)^2     //
//     If i= m, j != k or m,                                                  //
//        a"[m][j] = a"[j][m] = a[m][j] * cos(phi) - a[k][j] * sin(phi)       //
//                                                                            //
//     If X is the matrix of normalized eigenvectors stored so that the ith   //
//     column corresponds to the ith eigenvalue, then AX = X Lamda, where     //
//     Lambda is the diagonal matrix with the ith eigenvalue stored at        //
//     Lambda[i][i], i.e. X'AX = Lambda and X is orthogonal, the eigenvectors //
//     are normalized and orthogonal.  So, X = Q1 Q2 ... Qs, where Qi is      //
//     the ith orthogonal matrix,  i.e. X can be recursively approximated by  //
//     the recursion relation X" = X Q, where Q is the orthogonal matrix and  //
//     the initial estimate for X is the identity matrix.                     //
//     If j = k, then x"[i][k] = x[i][k] * cos(phi) + x[i][m] * sin(phi),     //
//     if j = m, then x"[i][m] = x[i][m] * cos(phi) - x[i][k] * sin(phi), and //
//     if j != k and j != m, then x"[i][j] = x[i][j].                         //
//                                                                            //
//  Arguments:                                                                //
//     double  eigenvalues                                                    //
//        Array of dimension n, which upon return contains the eigenvalues of //
//        the matrix A.                                                       //
//     double* eigenvectors                                                   //
//        Matrix of eigenvectors, the ith column of which contains an         //
//        eigenvector corresponding to the ith eigenvalue in the array        //
//        eigenvalues.                                                        //
//     double* A                                                              //
//        Pointer to the first element of the symmetric n x n matrix A. The   //
//        input matrix A is modified during the process.                      //
//     int     n                                                              //
//        The dimension of the array eigenvalues, number of columns and rows  //
//        of the matrices eigenvectors and A.                                 //
//                                                                            //
//  Return Values:                                                            //
//     Function is of type void.                                              //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], double eigenvalues[N], double eigenvectors[N][N]       //
//                                                                            //
//     (your code to initialize the matrix A )                                //
//                                                                            //
//     Jacobi_Cyclic_Method(eigenvalues, (double*)eigenvectors,               //
//                                                          (double *) A, N); //
////////////////////////////////////////////////////////////////////////////////

void Jacobi_Cyclic_Method(double eigenvalues[], double *eigenvectors, double *A, int n)
{
   //int row;
   int i, j, k, m;
   double *pAk, *pAm, *p_r, *p_e;
   double threshold_norm, threshold;
   double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
   double sin_2phi, cos_2phi, cot_2phi;
   double dum1, dum2, dum3;
   //double r;
   double max;
   // Take care of trivial cases
   if ( n < 1) return;
   if ( n == 1) {
      eigenvalues[0] = *A;
      *eigenvectors = 1.0;
      return;
   }
   // Initialize the eigenvalues to the identity matrix.
   for (p_e = eigenvectors, i = 0; i < n; i++)
      for (j = 0; j < n; p_e++, j++)
         if (i == j) *p_e = 1.0; else *p_e = 0.0;
   // Calculate the threshold and threshold_norm.
   for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++) 
      for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
   threshold = sqrt(threshold + threshold);
   threshold_norm = threshold * DBL_EPSILON;
   max = threshold + 1.0;
   while (threshold > threshold_norm) {
      threshold /= 10.0;
      if (max < threshold) continue;
      max = 0.0;
      for (pAk = A, k = 0; k < (n-1); pAk += n, k++) {
         for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++) {
            if ( fabs(*(pAk + m)) < threshold ) continue;
            // Calculate the sin and cos of the rotation angle which annihilates A[k][m].
            cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
            dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
            if (cot_2phi < 0.0) dum1 = -dum1;
            tan_phi = -cot_2phi + dum1;
            tan2_phi = tan_phi * tan_phi;
            sin2_phi = tan2_phi / (1.0 + tan2_phi);
            cos2_phi = 1.0 - sin2_phi;
            sin_phi = sqrt(sin2_phi);
            if (tan_phi < 0.0) sin_phi = - sin_phi;
            cos_phi = sqrt(cos2_phi); 
            sin_2phi = 2.0 * sin_phi * cos_phi;
            cos_2phi = cos2_phi - sin2_phi;
            // Rotate columns k and m for both the matrix A and the matrix of eigenvectors.
            p_r = A;
            dum1 = *(pAk + k);
            dum2 = *(pAm + m);
            dum3 = *(pAk + m);
            *(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
            *(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
            *(pAk + m) = 0.0;
            *(pAm + k) = 0.0;
            for (i = 0; i < n; p_r += n, i++) {
               if ( (i == k) || (i == m) ) continue;
               if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
               if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
               dum3 = dum1 * cos_phi + dum2 * sin_phi;
               if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
               dum3 = - dum1 * sin_phi + dum2 * cos_phi;
               if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
            }
            for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++) {
               dum1 = *(p_e + k);
               dum2 = *(p_e + m);
               *(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
               *(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
            }
         }
         for (i = 0; i < n; i++)
            if ( i == k ) continue;
            else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
      }
   }
   for (pAk = A, k = 0; k < n; pAk += n, k++) eigenvalues[k] = *(pAk + k); 
}



/// MATRIX FUNCTIONS ///

double trace(double *m) {
	return m[0 +3* 0] + m[1 +3* 1] + m[2 +3* 2];
}

double determinant(double *m) {
  return -m[2] * m[4] * m[6] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] 
    - m[0] * m[5] * m[7] - m[1] * m[3] * m[8] + m[0] * m[4] * m[8];
}

void matrixIdentity(double *m) {
	m[0] = 1.0;
	m[1] = 0.0;
	m[2] = 0.0;
	m[3] = 0.0;
	m[4] = 1.0;
	m[5] = 0.0;
	m[6] = 0.0;
	m[7] = 0.0;
	m[8] = 1.0;
}

void matrixZero(double *m) {
	int i;
	for (i = 0; i < 9; i += 1)
		m[i] = 0.0;
}

void matrixSetDiagonal(double m11, double m22, double m33, double *out) {
    out[0 +3* 0] = m11;
    out[1 +3* 0] = 0.0;
    out[2 +3* 0] = 0.0;
    out[0 +3* 1] = 0.0;
    out[1 +3* 1] = m22;
    out[2 +3* 1] = 0.0;
    out[0 +3* 2] = 0.0;
    out[1 +3* 2] = 0.0;
    out[2 +3* 2] = m33;
}

void matrixCopy(double *m, double *out) {
    int i, j;
    for (i = 0; i < 3; i += 1)
        for (j = 0; j < 3; j += 1)
            out[i +3* j] = m[i +3* j];
}

void matrixScale(double c, double *m, double *out) {
    int i, j;
    for (i = 0; i < 3; i += 1)
        for (j = 0; j < 3; j += 1)
            out[i +3* j] = c * m[i +3* j];
}

void matrixAdd(double *m, double *n, double *out) {
    int i, j;
    for (i = 0; i < 3; i += 1)
        for (j = 0; j < 3; j += 1)
            out[i +3* j] = m[i +3* j] + n[i +3* j];
}

void matrixSubtract(double *m, double *n, double *out) {
    int i, j;
    for (i = 0; i < 3; i += 1)
        for (j = 0; j < 3; j += 1)
            out[i +3* j] = m[i +3* j] - n[i +3* j];
}

// Input and output may safely alias each other.
void matrixTranspose(double *m, double *out) {
    double m12, m13, m23;
    m12 = m[0 +3* 1];
    m13 = m[0 +3* 2];
    m23 = m[1 +3* 2];
    out[0 +3* 0] = m[0 +3* 0];
    out[0 +3* 1] = m[1 +3* 0];
    out[0 +3* 2] = m[2 +3* 0];
    out[1 +3* 0] = m12;
    out[1 +3* 1] = m[1 +3* 1];
    out[1 +3* 2] = m[2 +3* 1];
    out[2 +3* 0] = m13;
    out[2 +3* 1] = m23;
    out[2 +3* 2] = m[2 +3* 2];
}

// Input and output may safely alias each other.
void matrixMatrixMultiply(double *m, double *n, double *out) {
    double p[9];
    int i, j;
    for (i = 0; i < 3; i += 1)
        for (j = 0; j < 3; j += 1)
            p[i +3* j] = m[i +3* 0] * n[0 +3* j] + m[i +3* 1] * n[1 +3* j] + m[i +3* 2] * n[2 +3* j];
    matrixCopy(p, out);
}

// Returns the Frobenius inner product. Take the square root to get the norm.
double frobenius(double *m, double *n) {
	double mT[9], mTn[9];
	matrixTranspose(m, mT);
	matrixMatrixMultiply(mT, n, mTn);
	return trace(mTn);
}

// Sets out to the outer product of column vectors, v w^T.
void matrixOuter(double *v, double *w, double *out) {
  int i, j;
  for (i = 0; i < 3; i += 1)
    for (j = 0; j < 3; j += 1)
      out[i +3* j] = v[i] * w[j];
}

// Computes the eigenvalues and eigenvectors of a symmetric matrix.
// Input: dimension n, n x n symmetric matrix X.
// Output: n x n matrix eigvec, whose columns are the eigenvectors. n-array eigval holding eigenvalues.
// Returns: 0 on success, nonzero on various kinds of error.
// Notes: All matrices are stored in column-major order. 
// The contents of X are not changed, unless X aliases eigvec, which is otherwise okay.
// The user must pre-allocate eigvec, eigval before calling this function.
// Uses Jacobi cyclic method.
// Allocates scratch space to hold a copy of X. Memory allocation errors are the only kind possible.
int symmetricEigensystem(int n, double *X, double *eigvec, double *eigval) {
  if (n == 3) {
    double x[9];
    matrixCopy(X, x);
    Jacobi_Cyclic_Method(eigval, eigvec, x, 3);
    matrixTranspose(eigvec, eigvec);
    return 0;
  } else {
    double *x;
    x = (double *)malloc(n * n * sizeof(double));
    if (x == NULL) {
      fprintf(stderr, "error: symmetricEigensystem: could not allocate memory");
      return 1;
    } else {
      int i, j;
      for (i = 0; i < n; i += 1)
        for (j = 0; j < n; j += 1)
          x[i +n* j] = X[i +n* j];
      Jacobi_Cyclic_Method(eigval, eigvec, x, n);
      // Transpose eigvec.
      double swap;
      for (i = 0; i < n; i += 1)
        for (j = 0; j < n; j += 1) {
          swap = eigvec[i +n* j];
          eigvec[i +n* j] = eigvec[j +n* i];
          eigvec[j +n* i] = swap;
        }
      free(x);
      return 0;
    }
  }
}

int symmetricInverse(double *ell, double *inv) {
  double vecs[9], vals[3], vecsT[9];
  if (symmetricEigensystem(3, ell, vecs, vals) != 0)
  	return 1;
	else if (vals[0] == 0.0 || vals[1] == 0.0 || vals[2] == 0.0)
    return 1;
  else {
		matrixSetDiagonal(1.0 / vals[0], 1.0 / vals[1], 1.0 / vals[2], inv);
		matrixMatrixMultiply(vecs, inv, inv);
		matrixTranspose(vecs, vecsT);
		matrixMatrixMultiply(inv, vecsT, inv);
		return 0;
	}
}

void printMatrix(double *m) {
    printf("%f %f %f\n", m[0], m[3], m[6]);
    printf("%f %f %f\n", m[1], m[4], m[7]);
    printf("%f %f %f\n", m[2], m[5], m[8]);
}

void printMatrixMathematica(double *m) {
    printf("{{%f, %f, %f}, {%f, %f, %f}, {%f, %f, %f}}", m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]);
}

void printMatrixR(double *m) {
    printf("matrix(c(%f, %f, %f, %f, %f, %f, %f, %f, %f), 3, 3)", m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
}



/// 3-DIMENSIONAL VECTOR FUNCTIONS ///

void vectorSet(double x, double y, double z, double *out) {
    out[0] = x;
    out[1] = y;
    out[2] = z;
}

void vectorCopy(double *v, double *out) {
    vectorSet(v[0], v[1], v[2], out);
}

void vectorAdd(double *v, double *w, double *out) {
    vectorSet(v[0] + w[0], v[1] + w[1], v[2] + w[2], out);
}

void vectorSubtract(double *v, double *w, double *out) {
    vectorSet(v[0] - w[0], v[1] - w[1], v[2] - w[2], out);
}

void vectorScale(double c, double *v, double *out) {
    vectorSet(c * v[0], c * v[1], c * v[2], out);
}

void vectorMatrixMultiply(double *m, double *v, double *out) {
    out[0] = m[0 +3* 0] * v[0] + m[0 +3* 1] * v[1] + m[0 +3* 2] * v[2];
    out[1] = m[1 +3* 0] * v[0] + m[1 +3* 1] * v[1] + m[1 +3* 2] * v[2];
    out[2] = m[2 +3* 0] * v[0] + m[2 +3* 1] * v[1] + m[2 +3* 2] * v[2];
}

double scalarProduct(double *v, double *w) {
    return v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
}

int vectorNormalize(double *v, double *out) {
    double norm;
    norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (norm != 0.0) {
        vectorScale(1.0 / norm, v, out);
        return 0;
    }
    else
        return 1;
}

void vectorCross(double *v, double *w, double *out) {
    vectorSet(v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2], v[0] * w[1] - v[1] * w[0], out);
}

void vectorsOrthogonalToVector(double *u, double *v, double *w) {
    double proj[3];
    if (u[0] > 0.5 || u[0] < -0.5)
        vectorSet(0.0, 1.0, 0.0, v);
    else
        vectorSet(1.0, 0.0, 0.0, v);
    vectorScale(scalarProduct(v, u) / scalarProduct(u, u), u, proj);
    vectorSubtract(v, proj, v);
    vectorNormalize(v, v);
    vectorCross(u, v, w);
}

void vectorMean(int n, double *vecs, double *out) {
	int i;
	vectorSet(0.0, 0.0, 0.0, out);
	for (i = 0; i < n; i += 1)
		vectorAdd(&(vecs[3 * i]), out, out);
	vectorScale(1.0 / n, out, out);
}



/// N-DIMENSIONAL VECTOR FUNCTIONS ///

void vectorSubtractN(int dim, double *v, double *w, double *out) {
  int i;
  for (i = 0; i < dim; i += 1)
    out[i] = v[i] - w[i];
}

double scalarProductN(int dim, double *v, double *w) {
  int i;
  double dot = 0.0;
  for (i = 0; i < dim; i += 1)
    dot += v[i] * w[i];
  return dot;
}



/// DIRECTIONAL STATISTICS ///

void binghamCenter(int n, double *vecs, double *out) {
  double sym[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double *v;
	int i;
	// Add up the outer products of the vectors.
	for (i = 0; i < n; i += 1) {
		v = &(vecs[3 * i]);
		sym[0] += v[0] * v[0];
		sym[1] += v[1] * v[0];
		sym[2] += v[2] * v[0];
		sym[4] += v[1] * v[1];
		sym[5] += v[2] * v[1];
		sym[8] += v[2] * v[2];
	}
	sym[3] = sym[1];
	sym[6] = sym[2];
	sym[7] = sym[5];
	// Compute the eigensystem.
	double eigvecs[9], eigvals[3];
	symmetricEigensystem(3, sym, eigvecs, eigvals);
	// Return the eigenvector associated to the largest eigenvalue.
	if (eigvals[0] >= eigvals[1] && eigvals[0] >= eigvals[2])
		vectorNormalize(&(eigvecs[0]), out);
	else if (eigvals[1] >= eigvals[0] && eigvals[1] >= eigvals[2])
		vectorNormalize(&(eigvecs[3]), out);
	else
		vectorNormalize(&(eigvecs[6]), out);
}


