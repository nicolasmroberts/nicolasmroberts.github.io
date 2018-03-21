


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This file offers many functions for computing with rotations of 3D space, borrowed from all over the orientation statistics literature. In this file, 'rotation matrix' means 'real 3x3 matrix that is special orthogonal'. That's our most common representation of rotations, although that may change in a future release.



### CONVERSIONS AMONG REPRESENTATIONS OF ROTATIONS ###

# The radius of the equal-volume plot. Roughly 0.62.
rotEqualVolumeRadius <- (3 / (4 * pi))^(1 / 3)

#' Conversion from anti-symmetric matrix (infinitesimal rotation) to a vector of its non-redundant entries.
#'
#' @param w A rotation matrix.
#' @return A real 3D vector of length between 0 and pi, inclusive.
rotVectorFromAntisymmetric <- function(w) {
	c(w[3, 2], w[1, 3], w[2, 1])
}

#' Conversion from vector representation to anti-symmetric matrix (infinitesimal rotation).
#'
#' @param v A real 3D vector.
#' @return A rotation matrix.
rotAntisymmetricFromVector <- function(v) {
	matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3)
}

#' Conversion from angle-axis representation to matrix.
#'
#' @param au A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians.
#' @return A rotation matrix.
rotMatrixFromAxisAngle <- function(ua) {
	m <- rotAntisymmetricFromVector(ua[1:3])
	diag(c(1, 1, 1)) + sin(ua[[4]]) * m + (1 - cos(ua[[4]])) * (m %*% m)
}

#' Conversion from matrix to angle-axis representation.
#'
#' @param r A rotation matrix.
#' @return A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians.
rotAxisAngleFromMatrix <- function(r) {
	cosine <- (tr(r) - 1) / 2
	if (cosine >= 1)
		c(0, 0, -1, 0)
	else {
		u <- squareRoot(1 + (r[1, 1] - 1) / (1 - cosine))
		if (r[3, 2] < r[2, 3])
			u <- -u
		v <- squareRoot(1 + (r[2, 2] - 1) / (1 - cosine))
		if (r[1, 3] < r[3, 1])
			v <- -v
		w <- squareRoot(1 + (r[3, 3] - 1) / (1 - cosine))
		if (r[2, 1] < r[1, 2])
			w <- -w
		nrm <- sqrt(u^2 + v^2 + w^2)
		a <- arcCos(cosine)
		c(u / nrm, v / nrm, w / nrm, a)
	}
}

#' Matrix exponentiation of infinitesimal rotation to finite rotation.
#'
#' @param w A 3x3 real matrix (anti-symmetric).
#' @return A rotation matrix.
rotExp <- function(w) {
	th <- sqrt(tr(crossprod(w, w)) / 2)
	if (th == 0)
    diag(c(1, 1, 1))
	else
    diag(c(1, 1, 1)) + (sin(th) / th) * w + ((1 - cos(th)) / th^2) * (w %*% w)
}

#' Matrix logarithm to produce infinitesimal from finite rotation.
#'
#' @param r A 3x3 rotation matrix.
#' @return A 3x3 real matrix (anti-symmetric). The principal logarithm.
rotLog <- function(r) {
	ua <- rotAxisAngleFromMatrix(r)
	ua[[4]] * rotAntisymmetricFromVector(ua[1:3])
}

#' Maps a tangent space into the space of rotations.
#' 
#' Converts tangent vector v into rotation matrix R = exp(V) C.
#' @param v A 3D real vector.
#' @param center A rotation matrix.
#' @return A rotation matrix.
rotMatrixFromRightTangent <- function(v, center) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    center
  else
    rotMatrixFromAngleAxis(c(v / nrm, nrm)) %*% center
}

#' Maps the space of rotations into one of its tangent spaces.
#' 
#' Converts rotation matrix R = exp(V) C into vector v. Appropriate only if the two given rotations are close to each other.
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @return A 3D real vector.
rotRightTangentFromMatrix <- function(r, center) {
  ua <- rotAxisAngleFromMatrix(r %*% t(center))
  ua[1:3] * ua[[4]]
}

#' Maps a tangent space into the space of rotations.
#' 
#' Converts tangent vector v into rotation matrix R = C exp(V).
#' @param v A 3D real vector.
#' @param center A rotation matrix.
#' @return A 3x3 rotation matrix.
rotMatrixFromLeftTangent <- function(v, center) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    center
  else
    center %*% rotMatrixFromAxisAngle(c(v / nrm, nrm))
}

#' Maps the space of rotations into one of its tangent spaces.
#' 
#' Converts rotation matrix R = exp(V) C into vector v. Appropriate only if the two given rotations are close to each other.
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @return A 3D real vector.
rotLeftTangentFromMatrix <- function(r, center) {
  ua <- rotAxisAngleFromMatrix(t(center) %*% r)
  ua[1:3] * ua[[4]]
}

#' Rotation matrix about the x-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
rotMatrixAboutX <- function(a) {
	cosine <- cos(a)
	sine <- sin(a)
	matrix(c(1, 0, 0, 0, cosine, sine, 0, -sine, cosine), 3, 3)
}

#' Rotation matrix about the y-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
rotMatrixAboutY <- function(a) {
	cosine <- cos(a)
	sine <- sin(a)
	matrix(c(cosine, 0, -sine, 0, 1, 0, sine, 0, cosine), 3, 3)
}

#' Rotation matrix about the z-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
rotMatrixAboutZ <- function(a) {
	cosine <- cos(a)
	sine <- sin(a)
	matrix(c(cosine, sine, 0, -sine, cosine, 0, 0, 0, 1), 3, 3)
}

#' Rotation matrix from xzx Euler angles.
#'
#' The resulting matrix represents rotation about the global x-axis through the angle cba[[3]], followed by rotation about the global z-axis through cba[[2]], followed by rotation about the global x-axis through cba[[1]].
#' @param cba A 3D real vector. The Euler angles in radians.
#' @return A rotation matrix.
rotMatrixFromXZXAngles <- function(cba) {
  rotMatrixAboutX(cba[[1]]) %*% rotMatrixAboutZ(cba[[2]]) %*% rotMatrixAboutX(cba[[3]])
}

#' Extraction of xzx Euler angles from rotation matrix.
#'
#' A vector cba is produced, such that the matrix represents rotation about the global x-axis through the angle cba[[3]], followed by rotation about the global z-axis through cba[[2]], followed by rotation about the global x-axis through cba[[1]].
#' @param r A rotation matrix.
#' @return A 3D real vector. The Euler angles in radians. The middle angle is always in [0, pi].
rotXZXAnglesFromMatrix <- function(r) {
  bb <- arcCos(r[1, 1])
  aa <- atan2(r[1, 3], -r[1, 2])
  cc <- atan2(r[3, 1], r[2, 1])
  c(cc, bb, aa)
}

#' Conversion of angle-axis representation to equal-volume representation.
#'
#' @param ua A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians.
#' @return A 3D real vector.
rotEqualVolumeFromAxisAngle <- function(ua) {
	rhoCubed <- (ua[[4]] - sin(ua[[4]])) * 3 / (4 * pi^2)
	rhoCubed^(1 / 3) * ua[1:3]
}

#' Solution of x - sin(x) == d, within tolerance of epsilon.
#'
#' Helper function for rotAxisAngleFromEqualVolume, etc.
#' @param d A real number.
#' @return A real number, hopefully in [0, pi].
rotXMinusSinXSolution <- function(d) {
	# On [0, pi], x^3 / 9 is close to x - sin(x). Use that to choose the seed.
	x0 <- (9.0 * d)^(1.0 / 3.0)
	# Proceed by custom Newton's method. (See also optim, nlm, nlminb.)
	x1 <- x0 - (x0 - sin(x0) - d) / (1.0 - cos(x0))
	while (abs(x1 - x0) > epsilon) {
		x0 <- x1
		x1 <- x0 - (x0 - sin(x0) - d) / (1.0 - cos(x0))
	}
	x1
}

#' Conversion of equal-volume representation to angle-axis representation.
#'
#' @param v A 3D real vector.
#' @return A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians.
rotAxisAngleFromEqualVolume <- function(v) {
	rho <- sqrt(dot(v, v))
	if (rho == 0)
    c(0, 0, -1, 0)
	else {
		a <- rotXMinusSinXSolution(rho^3 * 4 * pi^2 / 3)
    c(v / rho, a)
	}
}

#' Conversion of rotation matrix to its equal-volume representation.
#'
#' @param r A rotation matrix.
#' @return A 3D real vector.
rotEqualVolumeFromMatrix <- function(r) {
	rotEqualVolumeFromAxisAngle(rotAxisAngleFromMatrix(r))
}

#' Conversion of equal-volume representation to rotation matrix.
#'
#' @param v A 3D real vector.
#' @return A rotation matrix.
rotMatrixFromEqualVolume <- function(v) {
	rotMatrixFromAxisAngle(rotAxisAngleFromEqualVolume(v))
}

#' Conversion of rotation matrix to its equal-angle representation.
#'
#' @param r A rotation matrix.
#' @return A 3D real vector.
rotEqualAngleFromMatrix <- function(r) {
  ua <- rotAxisAngleFromMatrix(r)
  ua[1:3] * tan(ua[[4]] / 4)
}

#' Conversion of equal-angle representation to rotation matrix.
#'
#' @param v A 3D real vector.
#' @return A rotation matrix.
rotMatrixFromEqualAngle <- function(v) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    diag(c(1, 1, 1))
  else
    rotMatrixFromAxisAngle(c(v / nrm, 4 * atan(nrm)))
}

#' Conversion of rotation matrix to its Rodrigues representation.
#'
#' @param r A rotation matrix.
#' @return A 3D real vector.
rotRodriguesFromMatrix <- function(r) {
  ua <- rotAxisAngleFromMatrix(r)
  ua[1:3] * tan(ua[[4]] / 2)
}

#' Conversion of Rodrigues representation to rotation matrix.
#'
#' @param v A 3D real vector.
#' @return A rotation matrix.
rotMatrixFromRodrigues <- function(v) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    diag(c(1, 1, 1))
  else
    rotMatrixFromAxisAngle(c(v / nrm, 2 * atan(nrm)))
}

#' Conversion of axis-angle representation to quaternion.
#'
#' @param ua A 4D real vector. Unit vector u followed by angle a in radians.
#' @return A 4D real vector (unit length).
rotQuaternionFromAxisAngle <- function(ua) {
  c(cos(ua[[4]] / 2), sin(ua[[4]] / 2) * ua[1:3])
}

#' Conversion of quaternion representation to axis-angle.
#'
#' @param q A 4D real vector (unit length).
#' @return A 4D real vector. Unit vector u followed by angle a in radians.
rotAxisAngleFromQuaternion <- function(q) {
	a <- 2 * arcCos(q[[1]])
	sine <- sin(a / 2)
	v <- q[2:4]
	if (sine == 0)
		if (abs(v[[1]]) > abs(v[[2]]) && abs(v[[1]]) > abs(v[[3]]))
			if (v[[1]] > 0)
        c(1, 0, 0, a)
			else
        c(-1, 0, 0, a)
		else if (abs(v[[2]]) > abs(v[[1]]) && abs(v[[2]]) > abs(v[[3]]))
			if (v[[2]] > 0)
        c(0, 1, 0, a)
			else
        c(0, -1, 0, a)
		else
			if (v[[3]] > 0)
        c(0, 0, 1, a)
			else
        c(0, 0, -1, a)
	else
    c(v / sine, a)
}

#' Conversion of quaternion representation to its rotation matrix.
#'
#' @param q A 4D real vector (unit length).
#' @return A rotation matrix.
rotMatrixFromQuaternion <- function(q) {
	r11 <- q[[1]]^2 + q[[2]]^2 - q[[3]]^2 - q[[4]]^2
	r12 <- 2 * (q[[2]] * q[[3]] - q[[1]] * q[[4]])
	r13 <- 2 * (q[[1]] * q[[3]] + q[[2]] * q[[4]])
	r21 <- 2 * (q[[1]] * q[[4]] + q[[2]] * q[[3]])
	r22 <- q[[1]]^2 + q[[3]]^2 - q[[2]]^2 - q[[4]]^2
	r23 <- 2 * (q[[3]] * q[[4]] - q[[1]] * q[[2]])
	r31 <- 2 * (q[[2]] * q[[4]] - q[[1]] * q[[3]])
	r32 <- 2 * (q[[1]] * q[[2]] + q[[3]] * q[[4]])
	r33 <- q[[1]]^2 + q[[4]]^2 - q[[2]]^2 - q[[3]]^2
	matrix(c(r11, r21, r31, r12, r22, r32, r13, r23, r33), 3, 3)
}

#' Conversion of rotation matrix to its quaternion representation.
#'
#' @param r A rotation matrix.
#' @return A 4D real vector (unit length).
rotQuaternionFromMatrix <- function(r) {
	rotQuaternionFromAxisAngle(rotAxisAngleFromMatrix(r))
}



### GEODESIC METHODS ###

#' The distance between two rotations as points in SO(3).
#'
#' @param r A rotation matrix.
#' @param q A rotation matrix.
#' @return A real number (in the interval [0, pi]).
rotDistance <- function(r, q) {
  arcCos((tr(crossprod(r, q)) - 1) / 2)
}

#' Diameter of a set of rotations.
#' 
#' @param rs A list of rotation matrices.
#' @return A real number (in the interval [0, pi]).
rotDiameter <- function(rs) {
  f <- function(i) {
    max(sapply(1:(i - 1), function(j) rotDistance(rs[[i]], rs[[j]])))
  }
  max(sapply(2:(length(rs)), f))
}

#' The smallest rotation matrix R such that R u = v.
#' 
#' @param u A ray (unit real 3D vector).
#' @param v A ray.
#' @return A rotation matrix.
rotSmallestRotationFromTwoRays <- function(u, v) {
  axis <- rayNormalized(cross(u, v))
  angle <- arcCos(dot(u, v))
  rotMatrixFromAxisAngle(c(axis, angle))
}

#' The smallest rotation matrix R such that R u = v or R u = -v.
#' 
#' @param u A line (unit real 3D vector).
#' @param v A line.
#' @return A rotation matrix.
rotSmallestRotationFromTwoLines <- function(u, v) {
  if (dot(u, v) < 0)
    rotSmallestRotationFromTwoRays(u, -v)
  else
    rotSmallestRotationFromTwoRays(u, v)
}

#' The minimum distance between two elements of a set of rotations.
#'
#' @param rs A list of rotation matrices.
#' @return A real number (in the interval [0, pi]).
rotSeparation <- function(rs) {
  f <- function(i) {
    min(sapply(1:(i - 1), function(j) rotDistance(rs[[i]], rs[[j]])))
  }
  min(sapply(2:(length(rs)), f))
}

#' The Frechet (geodesic L^2) variance of a set of rotations about a given rotation.
#'
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix.
#' @return A real number (between 0 and pi^2 / 2).
rotVariance <- function(rs, center) {
  dists <- sapply(rs, rotDistance, center)
  sum(dists^2) / (2 * length(rs))
}

#' The Frechet (geodesic L^2) mean of a set of rotations, starting from a seed.
#'
#' An interative algorithm for computing the Frechet mean --- the rotation that minimizes the Frechet variance. The iterations continue until error squared of epsilon is achieved or numSteps iterations have been used. Try multiple seeds, to improve your chances of finding the global optimum.
#' @param rs A list of rotation matrices.
#' @param numSeeds A real number (positive integer). How many rs to try as seeds.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $mean (a special orthogonal real 3x3 matrix), $variance (a real number), $changeSquared (a real number), and $numSteps (a non-negative integer). changeSquared is the square of the size of the final step. numSteps is the number of iterations used.
rotMeanVariance <- function(rs, numSeeds=1, numSteps=100) {
  seeds <- sample(rs, numSeeds)
  # No variance is ever as large as 5.
  best <- c(5)
  for (seed in seeds) {
    rBar <- seed
    changeSquared <- epsilon + 1.0
    k <- 0
    while (changeSquared >= epsilon && k < numSteps) {
      w <- arithmeticMean(lapply(rs, function(r) rotLog(crossprod(rBar, r))))
      rBar <- rBar %*% rotExp(w)
      changeSquared <- tr(crossprod(w, w))
      k <- k + 1
    }
    var <- rotVariance(rs, rBar)
    if (var < best[[1]])
      best <- list(var, rBar, changeSquared, k)
  }
  list(variance=best[[1]], mean=best[[2]], changeSquared=best[[3]], numSteps=best[[4]])
}

#' The Frechet (geodesic L^2) mean. Convenience shortcut for rotMeanVariance.
#' 
#' @param rs A list of rotation matrices.
#' @param numSeeds A real number (positive integer). How many rs to try as seeds.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A rotation matrix.
rotFrechetMean <- function(rs, numSeeds=1, numSteps=100) {
  rotMeanVariance(rs, numSeeds=numSeeds, numSteps=numSteps)$mean
}

#' Point on a geodesic that is closest to a given point.
#' 
#' Returns the point on the geodesic (exp a M) B that is closest to Q.
#' @param q A rotation matrix.
#' @param m A 3x3 real matrix (antisymmetric).
#' @param b A rotation matrix.
#' @return A rotation matrix.
rotNearestPointOnRightGeodesic <- function(q, m, b) {
  v <- rotVectorFromAntisymmetric(m)
  vNorm <- sqrt(dot(v, v))
  if (vNorm == 0)
    b
  else {
    u <- m / vNorm
    uBQT <- u %*% b %*% t(q)
    alpha <- atan2(-tr(uBQT), tr(u %*% uBQT))
    r <- rotExp(alpha * u) %*% b
    s <- rotExp((alpha + pi) * u) %*% b
    if (rotDistance(q, r) < rotDistance(q, s))
      r
    else
      s
  }
}

#' Point on a geodesic that is closest to a given point.
#' 
#' Returns the point on the geodesic B (exp a M) that is closest to Q.
#' @param q A rotation matrix.
#' @param m A 3x3 real matrix (antisymmetric).
#' @param b A rotation matrix.
#' @return A rotation matrix.
rotNearestPointOnLeftGeodesic <- function(q, m, b) {
  v <- rotVectorFromAntisymmetric(m)
  vNorm <- sqrt(dot(v, v))
  if (vNorm == 0)
    b
  else {
    u <- m / vNorm
    qTBU <- t(q) %*% b %*% u
    alpha <- atan2(-tr(qTBU), tr(qTBU %*% u))
    r <- b %*% rotExp(alpha * u)
    s <- b %*% rotExp((alpha + pi) * u)
    if (rotDistance(q, r) < rotDistance(q, s))
      r
    else
      s
  }
}
# Here's some test code.
#b <- rotUniform()
#u <- rayUniform()
#rs <- lapply(0:360, function(i) b %*% rotMatrixFromAxisAngle(c(u, i * degree)))
#m <- rotAntisymmetricFromVector(1.34 * u) # arbitrary number
#q <- rotUniform()
#s <- rotNearestPointOnLeftGeodesic(q, m, b)
#min(sapply(rs, function(r) rotDistance(q, r)))
#rotDistance(q, s)
#rotEqualVolumePlot(list(b, s, q), list(rs, rotGeodesicPoints(q, s, 100)))

#' Points evenly spaced on a geodesic from one rotation to another.
#'
#' Doesn't work well if rotations are pi away from each other. Give it an intermediate point, to help it out.
#' @param r A rotation matrix.
#' @param q A rotation matrix.
#' @param numSteps A real number (positive integer).
#' @return A list of rotation matrices, of length numSteps + 1. The first one is r and the last one is q.
rotGeodesicPoints <- function(r, q, numSteps) {
  ua <- rotAxisAngleFromMatrix(r %*% t(q))
  lapply(0:numSteps, function(i) rotMatrixFromAxisAngle(c(ua[1:3], ua[[4]] * i / numSteps)) %*% q)
}



### PROJECTED ARITHMETIC MEAN AND RELATED COMPUTATIONS ###

#' Projects matrices near SO(3) onto SO(3).
#'
#' @param m A real 3x3 matrix. Presumed to be nearly special orthogonal.
#' @return A rotation matrix.
rotProjectedMatrix <- function(m) {
  # Q D^(1 / 2) Q^T = (M^T M)^(1 / 2).
	valsvecs <- eigen(crossprod(m, m), symmetric=TRUE)
	q <- valsvecs$vectors
	dSqrtInv <- valsvecs$values^(-1 / 2)
	# projection = M Q D^(-1 / 2) Q^T.
	m %*% q %*% diag(dSqrtInv) %*% t(q)
}

#' The projected arithmetic mean of a set of rotations.
#'
#' This function is appropriate only if the rotations are already known to be clustered about a central tendency (rather than girdled, say). In this case the projected arithmetic mean equals the MLE of the matrix Fisher mean (rotFisherMLE) and the quaternionic axial mean (rotMeanScatter). If the data are not clustered, then the projected arithmetic mean may equal the negation of the quaternionic axial mean (so orthogonal but not special orthogonal).
#' @param rs A list of rotation matrices.
#' @return A rotation matrix.
rotProjectedMean <- function(rs) {
  rotProjectedMatrix(arithmeticMean(rs))
}

#' Mean and dispersion of a sample, computed via axial treatment of quaternions.
#' 
#' See Prentice (1986), p. 218.
#' @param rs A list of rotation matrices.
#' @return A list consisting of $values (real 4-vector, non-negative, descending order, sum to 1) and $rotations (list of four special orthogonal real 3x3 matrices). If val1 + val4 > 0.5, then the sample is bipolar and $rotations[[1]] equals the projected arithmetic mean and the MLE of the matrix Fisher mean. If val1 + val4 < 0.5, then the sample is equatorial and $rotations[[4]] equals the MLE of the matrix Fisher mean.
rotMeanScatter <- function(rs) {
  qs <- lapply(rs, rotQuaternionFromMatrix)
  ts <- lapply(qs, function(q) {outer(q, q)})
  tMatrix <- arithmeticMean(ts)
  eig <- eigen(tMatrix, symmetric=TRUE)
  rots <- lapply(1:4, function(j) rotMatrixFromQuaternion(eig$vectors[,j]))
  list(values=eig$values, rotations=rots)
}

# Helper function for rotationSVD. Given permutation P of 1:n,
# returns matrix M of 0s and 1s such that for all V, (M V)[i] == V[P[i]].
rotOrthogonalFromPermutation <- function(p) {
	n <- length(p)
	m <- matrix(0, n, n)
	for (i in p)
		m[i, p[i]] <- 1
	m
}

#' Signed singular value decomposition, using only rotations but allowing one negative singular value.
#'
#' @param m A real 3x3 matrix.
#' @return A list of three 3x3 real matrices $u, $d, $v, such that M = U D V^T. U and V are rotation matrices and D is diagonal with |D_11| >= D_22 >= D_33 >= 0. If M is an arithmetic mean of rotation matrices, then furthermore 1 >= |D_11|.
rotSingularValueDecomposition <- function(m) {
	# M = U D V^T.
	duv <- svd(m)
	u <- duv$u
	vT <- t(duv$v)
	# Ensure that the singular values are in decreasing order.
	p <- rotOrthogonalFromPermutation(order(duv$d, decreasing=TRUE))
	u <- u %*% t(p)
	d <- p %*% diag(duv$d) %*% t(p)
	vT <- p %*% vT
	# Ensure that det U = det V = 1.
	if (det(u) < 0.0) {
		u[,1] <- -u[,1]
		d[1, 1] <- -d[1, 1]
	}
	if (det(vT) < 0.0) {
		vT[1,] <- -vT[1,]
		d[1, 1] <- -d[1, 1]
	}
	list(u=u, d=d, v=t(vT))
}

#' Maximum likelihood estimation of the matrix Fisher parameters M and K.
#'
#' The method proceeds by minimizing a certain function using numerical methods. Things can go wrong. The results include two kinds of diagnostic information. !!Warning: Don't be surprised if this doesn't work very well yet. The seeds need better choosing. Also, we should constrain the optimization to prevent the integrals from getting insanely big.
#' @param rs A list of rotation matrices.
#' @param seeds A list of 3-dimensional real vectors, or NULL. Seeds for the minimization. The details are complicated. If NULL is given, then this function picks 6 of them automatically.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $mHat (rotation matrix), $kHat (symmetric 3x3 real matrix), $error (real number), $minEigenvalue (real number). $error is 0 if and only if the minimization succeeds. If $error is 1, then try increasing numSteps. $minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive.
rotFisherMLE <- function(rs, seeds=NULL, numSteps=1000) {
  # Shockingly, the data are boiled down to three numbers for all of the hard stuff.
	udv <- rotSingularValueDecomposition(arithmeticMean(rs))
  d <- diag(udv$d)
	# Here is the function to be minimized.
	f <- function(logP) {
    p <- exp(logP)
	  s <- sort(p, decreasing=TRUE)
	  f <- function(u) {
	    besselI((s[1] - s[2]) * u, 0) * besselI((s[1] + s[2]) * (1 - u), 0) * exp(s[3] * (1 - 2 * u))
	  }
		log(integrate(f, 0, 1)$value) - p %*% d
	}
  # Choose seeds based on the first octant of the 24 chambers of Sei et al. (2013, p. 448).
  if (is.null(seeds))
    seeds <- list(c(-1, 0, 1), c(-1, 1, 0), c(0, -1, 1), c(0, 1, -1), c(1, -1, 0), c(1, 0, -1))
  # Try all of the seeds and report the best answer found.
	solutions <- lapply(seeds, function(seed) optim(seed, f, hessian=TRUE, control=list(maxit=numSteps)))
  values <- sapply(solutions, function(solution) solution$value)
  solution <- solutions[[which.min(values)]]
	mHat <- udv$u %*% t(udv$v)
	kHat <- udv$v %*% diag(exp(solution$par)) %*% t(udv$v)
	eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
	list(mHat=mHat, kHat=kHat, error=solution$convergence, minEigenvalue=min(eigvals))
}

#' P-value function for whether the data come from a matrix Fisher distribution with the hypothesized mean.
#'
#' Based on Downs (1972; Eq. 5.7, p = 2 case). Assumes large sample size and tightly clustered sample. Uses the Stiefel manifold version of the matrix Fisher distribution, rather than the SO(3) version (Sei et al. (2013)). Only the left 3x2 submatrices of the given matrices are used. Indeed, the matrices may safely be given as 3x2.
#' @param rs A list of rotation matrices.
#' @return An R function from {rotation matrices} to {real numbers}. For any given hypothesized mean, this function produces the p-value for that hypothesis.
rotDownsInference <- function(rs) {
  n <- length(rs)
  r32s <- lapply(rs, function(a) {a[1:3, 1:2]})
  # hHat = (rBar^T rBar)^(1 / 2).
  rBar <- arithmeticMean(r32s)
  eig <- eigen(crossprod(rBar, rBar), symmetric=TRUE)
  hHat <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  top <- 2.0 * (diag(c(1.0, 1.0)) - hHat)
  top <- sqrt(tr(crossprod(top, top)))
  f <- function(r) {
    r32 <- r[1:3, 1:2]
    bottom <- diag(c(2, 2)) - crossprod(rBar, r32) - crossprod(r32, rBar)
    xi <- top / sqrt(tr(crossprod(bottom, bottom)))
    q <- (1 / sqrt(xi) - 1) * (n - 5 / 3)
    1 - pf(q, 3, 3 * n - 5)
  }
  f
}

#' P-value function for whether the data come from a distribution with the hypothesized mean.
#'
#' Based on sampling theory of the moment of inertia matrix (Prentice, 1986).
#' @param rs A list of rotation matrices.
#' @return An R function from {rotation matrices} to {real numbers union NULL}. If the data are unsuitable (because they are not bipolar), then NULL is always output. Otherwise, for any given hypothesized mean, this function produces the p-value for that hypothesis.
rotPrenticeInference <- function(rs) {
  # Convert to quaternions in unusual order.
  qs <- lapply(rs, rotQuaternionFromMatrix)
  qs <- lapply(qs, function(q) c(q[[2]], q[[3]], q[[4]], q[[1]]))
  # Compute the lambdai and Ahat in ascending order.
  qqTs <- lapply(qs, function(q) outer(q, q))
  eig <- eigen(arithmeticMean(qqTs), symmetric=TRUE)
  vals <- rev(eig$values)
  vecs <- cbind(eig$vectors[,4], eig$vectors[,3], eig$vectors[,2], eig$vectors[,1])
  if (vals[[1]] + vals[[4]] <= 0.5) {
    # The data are not bipolar, so we return this garbage p-value function.
    f <- function(r) {NULL}
  }
  else {
    # The data are bipolar. Ensure that vecs has determinant 1.
    if (det(vecs) < 0)
      vecs[,4] <- -vecs[,4]
    # Estimate the center M of the distribution.
    #mHat <- vecs[,4]
    #mHat <- rotationFromQuaternion(c(mHat[[4]], mHat[[1]], mHat[[2]], mHat[[3]]))
    # Approximate the confidence region metric.
    x2s <- lapply(qqTs, function(qqT) diag(t(vecs) %*% qqT %*% vecs))
    cHat <- arithmeticMean(lapply(x2s, function(x2) outer(x2, x2)))
    aHat4By3 <- vecs[c(1, 2, 3, 4), c(1, 2, 3)]
    fHat <- diag(c(
      (vals[[4]] - vals[[1]])^2 / cHat[4, 1],
      (vals[[4]] - vals[[2]])^2 / cHat[4, 2],
      (vals[[4]] - vals[[3]])^2 / cHat[4, 3]))
    fHat <- length(rs) * aHat4By3 %*% fHat %*% t(aHat4By3)
    # Put the metric back into my preferred quaternion order.
    fHat <- rbind(fHat[4,], fHat[1,], fHat[2,], fHat[3,])
    fHat <- cbind(fHat[,4], fHat[,1], fHat[,2], fHat[,3])
    f <- function(r) {
      q <- rotQuaternionFromMatrix(r)
      1 - pchisq(as.numeric(q %*% fHat %*% q), 3)
    }
  }
  f
}

#' Test of uniformity, based on the Rayleigh statistic.
#'
#' @param rs A list of rotation matrices.
#' @return A list with members $p, $rayleigh. $p is a real number (between 0 and 1), the p-value for the test. Low values indicate that the sample was not drawn from the uniform distribution. $rayleigh is the Rayleigh statistic that produced $p.
rotRayleighInference <- function(rs) {
  rBar <- arithmeticMean(rs)
  rayleigh <- 3 * length(rs) * tr(crossprod(rBar, rBar))
	p <- 1 - pchisq(rayleigh, 9)
  list(p=p, rayleigh=rayleigh)
}



### TANGENT SPACE METHODS ###

#' Sample covariance matrix, approximated in the tangent space at a given rotation.
#'
#' Appropriate only if the sample is tightly concentrated near the center.
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Typically the Frechet mean of the rs.
#' @return A 3x3 real matrix (symmetric, non-negative eigenvalues).
rotLeftCovariance <- function(rs, center) {
	vs <- lapply(rs, rotLeftTangentFromMatrix, center)
	ms <- lapply(vs, function(v) {outer(v, v)})
	arithmeticMean(ms)
}

#' Principal geodesic analysis in the tangent space at a given rotation.
#'
#' Appropriate only if the sample is tightly concentrated near the center.
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Typically the Frechet mean of the rs.
#' @param numPoints A real number (integer, 0 or >= 3). The number of points to return on each of the six geodesics through the center.
#' @return A list consisting of $magnitudes (3D real vector, nonnegative) and $directions (3x3 real matrix, whose columns are unit-length vectors). The $magnitudes are in decreasing order. The $directions are the corresponding directions, suitable for use in rotMatrixFromLeftTangent. If numPoints >= 1, then there is also a $curves field (list of three lists of (2 numPoints + 1) rotation matrices).
rotLeftPrincipalComponentAnalysis <- function(rs, center, numPoints=0) {
	eig <- eigen(rotLeftCovariance(rs, center), symmetric=TRUE)
  mags <- sqrt(eig$values)
  dirs <- eig$vectors
	result <- list(magnitudes=mags, directions=dirs)
  if (numPoints >= 1) {
    f <- function(s, magDir) {
      rotMatrixFromLeftTangent(s / numPoints * magDir, center)
    }
    curve1 <- lapply(-numPoints:numPoints, f, mags[[1]] * dirs[,1])
    curve2 <- lapply(-numPoints:numPoints, f, mags[[2]] * dirs[,2])
    curve3 <- lapply(-numPoints:numPoints, f, mags[[3]] * dirs[,3])
    result$curves <- list(curve1, curve2, curve3)
  }
  result
}

#' Mahalanobis distance of a rotation relative to a sample.
#'
#' Appropriate only if the sample is tightly clustered and the given rotation is close to its center.
#' @param r A rotation matrix.
#' @param center A rotation matrix. Typically the Frechet mean of a sample.
#' @param covarInv A 3x3 real matrix (symmetric, positive-definite). Typically the inverse of the covariance obtained from rotLeftCovariance.
#' @return A real number.
rotMahalanobisNorm <- function(r, center, leftCovarInv) {
	v <- rotLeftTangentFromMatrix(r, center)
	sqrt(v %*% leftCovarInv %*% v)
}

#' P-value function for any hypothesized mean, using tangent space approximation and Mahalanobis distances.
#' 
#' Appropriate only if the rotations are tightly clustered. May fail if the rotations live on a geodesic curve or surface.
#' @param rs A list of rotation matrices. Typically the result of bootstrapping, MCMC, etc.
#' @param center A rotation matrix. Typically the Frechet mean of the rs.
#' @return A list containing elements $pvalue (an R function from {rotation matrices} to {real numbers}; for any given hypothesized mean, this function produces the p-value for that hypothesis), $center (as passed to this function), $leftCovarInv (real symmetric 3x3 matrix), and $q000, $q025, $q050, $q075, $q095, $q099, $q100 (real numbers; quantiles of Mahalanobis distance). For example, a rotation r with tangent vector v at center is in the 95% confidence region if sqrt(v %*% covarInv %*% v) < q095.
rotMahalanobisInference <- function(rs, center) {
  vs <- lapply(rs, rotLeftTangentFromMatrix, center)
  covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
  covarInv <- solve(covar)
  norms <- sapply(vs, function(v) {sqrt(v %*% covarInv %*% v)})
  empiricalCDF <- ecdf(norms)
  # Build the p-value function.
  f <- function(r) {
    v <- rotLeftTangentFromMatrix(r, center)
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  # Compute a few popular percentiles.
  qs <- quantile(norms, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(pvalue=f, center=center, leftCovarInv=covarInv,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' Ellipsoidal surface from Mahalanobis inference.
#' 
#' @param center A rotation matrix. Typically the $center from the inference.
#' @param leftCovarInv A 3x3 real matrix (symmetric, positive-definite). Typically the $leftCovarInv from the inference.
#' @param level A real number. Typically $q095^2 from the inference.
#' @param numNonAdapt A real number (non-negative integer). The number of refinements to the sphere that is deformed into the ellipsoid. Incrementing numNonAdapt improves visual quality but increases time and memory requirements by a factor of four.
#' @return A list of triangles, where each triangle is a list of three rotation matrices.
rotEllipsoidTriangles <- function(center, leftCovarInv, level, numNonAdapt=3) {
  # Diagonalize the inverse covariance.
  eig <- eigen(leftCovarInv, symmetric=TRUE)
  q <- eig$vectors
  a <- sqrt(level) * eig$values^(-0.5)
  # Make the ellipsoid in the left-invariant tangent space and transfer it to SO(3).
  sphere <- rayTetrahedralSphere(numNonAdapt)
  ellipsoid <- lapply(sphere, function(tri) lapply(tri, function(v) rotMatrixFromLeftTangent(q %*% (a * v), center)))
  ellipsoid
}

#' Approximation of one rotation in the tangent space at the other rotation.
#'
#' This function is based on Rancourt et al. (2000). It produces a somewhat different approximation from that of rotLeftTangentFromMatrix. Appropriate only if r and center are close to each other. Definitely inappropriate if the distance between r and center is greater than pi / 2.
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @return A 3D real vector.
rotRancourtFromMatrix <- function(r, center) {
	w <- (crossprod(center, r) - crossprod(r, center)) / 2.0
	rotVectorFromAntisymmetric(w)
}

#' Mapping a tangent space into the space of rotations.
#'
#' Inverse to rotRancourtFromMatrix.
#' @param v A 3D real vector.
#' @param center A rotation matrix.
#' @return A rotation matrix.
rotMatrixFromRancourt <- function(v, center) {
	cosine <- squareRoot(1 - dot(v, v))
	w <- rotAntisymmetricFromVector(v)
	r <- diag(c(1, 1, 1)) + w + (w %*% w) / (1 + cosine)
	center %*% r
}

#' P-value function for any hypothesized mean, using tangent space approximation of Rancourt et al. (2000).
#' 
#' Appropriate only if the rotations are tightly clustered. Definitely inappropriate if any rotation is more than pi / 2 away from the mean. Uses Eq. (2.6) of Rancourt et al. (2000), 'Using orientation statistics to investigate variations in human kinematics'. Inserts an extra factor of n on the left side of that equation. Without this extra factor, coverage rates are much too large. Also the corresponding author confirmed via e-mail on 2015/08/31 that the extra factor should be there.
#' @param rs A list of rotation matrices. The data.
#' @return An R function from {rotation matrices} to {real numbers union NA}. For any given hypothesized mean, returns NA if that mean is more than pi / 2 away from the mean of the rs, or the p-value if not.
rotRancourtInference <- function(rs) {
  udv <- rotSVD(arithmeticMean(rs))
  mHat <- udv$u %*% t(udv$v)
  vs <- lapply(rs, rotRancourtFromMatrix, mHat)
  ms <- lapply(vs, function(v) {outer(v, v)})
  n <- length(rs)
  s <- arithmeticMean(ms) * n / (n - 1)
  sInv <- solve(s)
  # This line has an extra n factor, compared to Eq. (2.6) of Rancourt et al. (2000).
  g <- n * (n - 3) / (3 * (n - 1)) * sInv
  f <- function(r) {
    if (rotDistance(r, mHat) > pi / 2)
      NA
    else {
      v <- rotRancourtFromMatrix(r, mHat)
      1 - pf(v %*% g %*% v, 3, n - 3)
    }
  }
  f
}



### SAMPLING FROM DISTRIBUTIONS ###

#' Sampling from the uniform distribution on the space of rotations.
#'
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotUniform <- function(n=NULL) {
	if (is.null(n)) {
		x <- rayNormalized(c(1, 0, 0) - rayUniform())
		refl <- diag(c(1, 1, 1)) - 2 * (x %o% x)
		a <- runif(1, 0, 2 * pi)
		-refl %*% rotMatrixAboutX(a)
	}
	else
    replicate(n, rotUniform(), simplify=FALSE)
}

#' Sampling from the wrapped trivariate normal distribution.
#'
#' @param s A rotation matrix. The center of the distribution.
#' @param kappa A real number (positive). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotWrappedTrivariateNormal <- function(s, kappa, n=NULL) {
	if (is.null(n)) {
		v <- rnorm(3, 0, 1 / kappa)
		s %*% rotExp(rotAntisymmetricFromVector(v))
	}
	else
    replicate(n, rotationWrappedTrivariateNormal(s, kappa), simplify=FALSE)
}

#' One rotation drawn from the matrix Fisher distribution on SO(3).
#'
#' This function is a helper function for rotFisher. Probably you do not want to call this function yourself.
#' @param lambda A 4D real vector. lambda from Kent et al. (2013).
#' @param omega A 4x4 real matrix. Omega from Kent et al. (2013).
#' @param sigma A 4x4 real matrix. Inverse of omega. Sigma from Kent et al. (2013).
#' @param bound A real number. M* from Kent et al. (2013).
#' @return A rotation matrix.
rotFisherHelper <- function(lambda, omega, sigma, bound) {
	x <- rayNormalized(mvrnorm(1, c(0, 0, 0, 0), sigma))
	w <- runif(1, 0, 1)
	# Repeat until w < f*(x) / (M* g*(x)).
	while (bound * w >= exp(-x %*% lambda %*% x) * (x %*% omega %*% x)^2) {
		x <- rayNormalized(mvrnorm(1, c(0, 0, 0, 0), sigma))
		w <- runif(1, 0, 1)
	}
	rotMatrixFromQuaternion(x)
}

#' Sampling from the matrix Fisher distribution on SO(3).
#'
#' Follows Kent et al. (2013), 'A new method to simulate the Bingham and related distributions in directional data analysis with applications'.
#' @param m A rotation matrices. The center of the distribution.
#' @param k A 3x3 real matrix (symmetric, positive-definite). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotFisher <- function(m, k, n=NULL) {
	if (is.null(n))
		rotFisher(m, k, 1)[[1]]
	else {
		# Compute the Lambdas for the Bingham distribution.
		eig <- eigen(k, symmetric=TRUE)
		d <- eig$values
		l1 <- 2 * (d[[2]] + d[[3]])
		l2 <- 2 * (d[[1]] + d[[3]])
		l3 <- 2 * (d[[1]] + d[[2]])
		lams <- c(0, l1, l2, l3)
		# Compute b as the greatest real root of a certain quartic polynomial.
		a0 <- 8 * l1 * l2 * l3
		a1 <- 8 * (l1 * l2 + l1 * l3 + l2 * l3 - l1 * l2 * l3)
		a2 <- 6 * (l1 + l2 + l3) - 4 * (l1 * l2 + l1 * l3 + l2 * l3)
		a3 <- 4 - 2 * (l1 + l2 + l3)
		a4 <- -1
		roots <- polyroot(c(a0, a1, a2, a3, a4))
		roots <- Filter(function (x) (Im(x)^2 < epsilon), roots)
		roots <- sapply(roots, function(x) Re(x))
		b <- max(roots)
		# Compute bound M*, Omega for ACG, and Sigma = Omega^-1 for normal.
		bound <- exp((b - 4) / 2) * (4 / b)^2
		oms <- c(1, 1, 1, 1) + 2 * lams / b
		sigs <- 1 / oms
		# Generate the random sample.
		v <- eig$vectors
		rs <- replicate(n, rotFisherHelper(diag(lams), diag(oms), diag(sigs), bound), simplify=FALSE)
		lapply(rs, function(r) (m %*% v %*% r %*% t(v)))
	}
}

#' Sampling from the isotropic matrix Fisher distribution.
#'
#' @param s A rotation matrix. The center of the distribution.
#' @param kappa A real number (positive). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotIsotropicFisher <- function(s, kappa, n=NULL) {
	rotFisher(s, diag((kappa^2 / 2) * c(1, 1, 1)), n)
}

#' Jeffreys prior for the concentration parameter of the isotropic matrix Fisher distribution on SO(3).
#'
#' See Bingham et al. (2010, Section 3.2). To avoid a notational conflict, rename their kappa to lambda. The relationship between lambda and the kappa of Qiu et al. (2013) is lambda == kappa^2 / 2. Then also let eta == -log kappa, so that kappa = exp(-eta). We work in terms of eta.
rotIsotropicFisherJeffreysPrior <- function(negLogKappa) {
  if (177 < negLogKappa)
    0
  else if (negLogKappa <= -3.433049)
    sqrt(6)
  else if (negLogKappa < -2.94) {
    # Linearly interpolate the function between -3.433049 and -2.94.
    y1 <- 2.451246 # == rotIsotropicFisherJeffreysPrior(-2.93)
    y0 <- 2.451212 # == rotIsotropicFisherJeffreysPrior(-2.94)
    m <- (y1 - y0) /  (-2.93 - -2.94)
    m * (negLogKappa + 2.94) + y0
  } else {
    #lam <- exp(-2 * negLogKappa) / 2
    #i0 <- besselI(2 * lam, 0)
    #i1 <- besselI(2 * lam, 1)
    #numer <- i0^2 * 2 / lam - i0 * i1 * 2 / lam^2 + i1^2 * (1 / lam^2 - 2 / lam)
    #denom <- (i0 - i1)^2
    #sqrt(numer / denom) * 2 * lam
    kSq <- exp(-2 * negLogKappa)
    i0 <- besselI(kSq, 0)
    i1 <- besselI(kSq, 1)
    numer <- (kSq - 1) * i0^2 - kSq * i1^2
    denom <- (i0 - i1)^2
    jeff <- 2 * sqrt(1 + numer / denom)
    jeff
  }
}

# This code reproduces the Fisher panel of Fig. 2 of Qiu et al. (2013).
#etas <- seq(from=-2, to=5, by=0.01)
#priors <- sapply(etas, rotIsotropicFisherJeffreysPrior)
#plot(x=etas, y=(priors / sqrt(6)))

# This code checks the behavior around 177. Looks good.
#etas <- seq(from=170, to=180, by=0.01)
#priors <- sapply(etas, rotIsotropicFisherJeffreysPrior)
#plot(x=etas, y=priors)

# This code checks the behavior around -3.433 and -2.94. Not smooth, but okay.
#etas <- seq(from=-4, to=-2, by=0.01)
#priors <- sapply(etas, rotIsotropicFisherJeffreysPrior)
#plot(x=etas, y=priors)



### BOOTSTRAPPING ###

#' Inference about the population mean, based on non-parametric bootstrapping.
#' 
#' This function bootstraps the rotation mean, returning two pieces of information. The first is a list of the bootstrapped means. The user should rotEqualAnglePlot or rotEqualVolumePlot this list, to make sure that the means form a fairly tight ellipsoidal cluster. If so, then the second piece of information may be used: An R function that takes a rotation R0 as input, as produces as output a p-value for the hypothesis that the mean of the population is R0. This p-value is the fraction of means that are farther from the mean of means than R0 is, based on Mahalanobis distance in the set of means.
#' @param rs A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisInference.
rotBootstrapInference <- function(rs, numBoots, func=rotFrechetMean) {
  boots <- replicate(numBoots, func(sample(rs, length(rs), replace=TRUE)), simplify=FALSE)
  bootMean <- func(boots)
  infer <- rotMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}

#' Inference about the population mean, based on parametric bootstrapping with the Fisher distribution.
#' 
#' Identical to rotBootstrapInference, but draws its bootstrap samples from the MLE matrix Fisher distribution for the data, rather than from the data set itself.
#' @param rs A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisInference.
rotFisherBootstrapInference <- function(rs, numBoots, func=rotFrechetMean) {
  mle <- rotFisherMLE(rs, numSteps=10000)
  boots <- replicate(numBoots, func(rotFisher(mle$mHat, mle$kHat, length(rs)), ...), simplify=FALSE)
  bootMean <- func(boots, ...)
  infer <- rotMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}

#' Two-sample inference about the difference in population means, based on non-parametric bootstrapping.
#' 
#' This function bootstraps the difference of the means of the two samples, returning two pieces of information. The first is a list of the bootstrapped differences. The user should rotEqualAnglePlot or rotEqualVolumePlot this list, to make sure that the differences form a fairly tight ellipsoidal cluster. If so, then the second piece of information may be used: An R function that takes a rotation R0 as input, as produces as output a p-value for the hypothesis that the difference of means is R0. This p-value is the fraction of differences that are farther from the mean of differences than R0 is, based on the Mahalanobis distance of the differences.
#' @param firsts A list of rotation matrices.
#' @param seconds A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisInference.
rotTwoSampleBootstrapInference <- function(firsts, seconds, numBoots, func=rotFrechetMean) {
  f <- function() {
    firstMean <- func(sample(firsts, length(firsts), replace=TRUE))
    secondMean <- func(sample(seconds, length(seconds), replace=TRUE))
    secondMean %*% t(firstMean)
  }
  boots <- replicate(numBoots, f(), simplify=FALSE)
  bootMean <- func(boots, ...)
  infer <- rotMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}

#' Two-sample inference about the difference in population means, based on parametric bootstrapping using the Fisher distribution.
#' 
#' Identical to rotTwoSampleBootstrapInference, but draws its bootstrap samples from matrix Fisher distributions MLE-fit to the data, rather than from the data set itself.
#' @param firsts A list of rotation matrices.
#' @param seconds A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisInference.
rotTwoSampleFisherBootstrapInference <- function(firsts, seconds, numBoots, func=rotFrechetMean) {
  firstMLE <- rotFisherMLE(firsts, numSteps=10000)
  secondMLE <- rotFisherMLE(seconds, numSteps=10000)
  f <- function() {
    firstMean <- func(rotationFisher(firstMLE$mHat, firstMLE$kHat, length(firsts)))
    secondMean <- func(rotationFisher(secondMLE$mHat, secondMLE$kHat, length(seconds)))
    secondMean %*% t(firstMean)
  }
  boots <- replicate(numBoots, f(), simplify=FALSE)
  bootMean <- func(boots, ...)
  infer <- rotMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}



### REGRESSION ###

#' Regression to fit a geodesic curve to data (x, R), with no pre-processing of the independent variable.
#'
#' Returns the best-fit geodesic R(x) = (exp (x m)) b for the given (x, r) pairs. Minimizes the sum of squared distances in the r-direction only, as in ordinary linear regression. Usually you want to try rotGeodesicRegression first.
#' @param xs A list of real numbers.
#' @param rs A list of rotation matrices of same length as xs.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $m (anti-symmetric 3x3 real matrix), $b (rotation matrix), $error (real number), $minEigenvalue (real number), $rSquared (real number), $prediction (R function). $error is 0 if and only if the minimization succeeds. If $error == 1, try increasing steps, to 1000 say. $minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive. $rSquared is the R^2 statistic measuring the amount of variance captured, from none (0) to all (1). $prediction takes an x as input and returns the predicted rotation R(x).
rotNativeGeodesicRegression <- function(xs, rs, numSteps=100) {
  # Let Q be the R whose x is closest to zero.
  q <- rs[[which.min(as.numeric(xs)^2)]]
  # Define the function to be minimized.
  n <- length(rs)
  e <- function(mw) {
    b <- rotExp(rotAntisymmetricFromVector(mw[4:6])) %*% q
    m <- rotAntisymmetricFromVector(mw[1:3])
    f <- function(i) {
      rotDistance(rs[[i]], rotExp(xs[[i]] * m) %*% b)^2
    }
    sum(sapply(1:n, f)) / (2 * n)
  }
  # Find the minimum, using the constant geodesic Q as the seed.
  seed <- c(0, 0, 0, 0, 0, 0)
  solution <- optim(seed, e, lower=c(-pi, -pi, -pi, 0), upper=c(pi, pi, pi, pi), hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  m <- rotAntisymmetricFromVector(solution$par[1:3])
  b <- rotExp(rotAntisymmetricFromVector(solution$par[4:6])) %*% q
  rSq <- 1 - solution$value / rotVariance(rs, rotFrechetMean(rs))
  prediction <- function(x) {
    rotExp(x * m) %*% b
  }
  list(b=b, m=m, error=solution$convergence, minEigenvalue=min(eigvals), rSquared=rSq, prediction=prediction)
}

#' Geodesic curve regression, with a convenient behind-the-scenes rescaling.
#' 
#' In theory, this function is equivalent to rotGeodesicRegression. In practice, it uses a behind-the-scenes rescaling of the xs onto the interval [0, 1], which I find improves the performance of the optimization. You should use this version of the function, unless you have some specific reason to use rotNativeGeodesicRegression.
#' @param xs A list of real numbers.
#' @param rs A list of rotation matrices of same length as xs.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $m (anti-symmetric 3x3 real matrix), $b (rotation matrix), $error (real number), $minEigenvalue (real number), $rSquared (real number), $prediction (R function). $error is 0 if and only if the minimization succeeds. If $error == 1, try increasing steps, to 1000 say. $minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive. $rSquared is the R^2 statistic measuring the amount of variance captured, from none (0) to all (1). $prediction takes an x as input and returns the predicted rotation R(x).
rotGeodesicRegression <- function(xs, rs, numSteps=100) {
  x0 <- min(xs)
  x1 <- max(xs)
  regr <- rotNativeGeodesicRegression(scales(xs), rs, numSteps)
  mNew <- regr$m / (x1 - x0)
  bNew <- rotExp(regr$m * -x0 / (x1 - x0)) %*% regr$b
  prediction <- function(x) {
    regr$prediction((x - x0) / (x1 - x0))
  }
  list(m=mNew, b=bNew, error=regr$error, minEigenvalue=regr$minEigenvalue, rSquared=regr$rSquared, prediction=prediction)
}

#' Permutation test for significance of a geodesic regression.
#' 
#' Returns up to numPerms R^2 values. May be fewer than numPerms, because of some regressions failing. Let n be the dimension of this vector and g the number of R^2 values greater than the R^2 for the original regression of the data. Let p = g / n. Small values of p indicate that the dependency detected by the regression is meaningful. Uses an internal rescaling, just like rotGeodesicRegression.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices, of same length as xs.
#' @param numPerms A real number (positive integer). The number of permutations to try.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A vector of real numbers, of dimension at most numPerms.
rotGeodesicRegressionPermutations <- function(xs, rs, numPerms, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- rotNativeGeodesicRegression(sample(ys, size=length(ys)), rs, numSteps)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

#' Kernel regression to fit a rotation R(x) to (x, R) data, with no pre-processing of the independent variable.
#' 
#' Theoretically identical to rotKernelRegression, but doesn't rescale the xs to the interval [0, 1].
#' @param x A real number. The x-value at which to predict the rotation R.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $r (the rotation matrix R(x)), $error (which should be 0; if it is 1, then increase numSteps), and $minEigenvalue (which should be positive; if not, then the minimization has failed).
rotNativeKernelRegression <- function(x, xs, rs, h, k=dnorm, numSteps=100) {
  # Let Q be the Ri whose xi is closest to x.
  q <- rs[[which.min((xs - x)^2)]]
  # Define the function to be minimized.
  kh <- function(x) {k(x / h) / h}
  e <- function(w) {
    r <- rotExp(rotAntisymmetricFromVector(w)) %*% q
    f <- function(i) {
      kh(x - xs[[i]]) * rotDistance(rs[[i]], r)^2
    }
    sum(sapply(1:length(xs), f)) / sum(sapply(xs, function(y) kh(x - y)))
  }
  # Minimize the function with respect to R, using Q as the seed.
  seed <- c(0, 0, 0)
  solution <- optim(seed, e, hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  r <- rotExp(rotAntisymmetricFromVector(solution$par)) %*% q
  list(r=r, error=solution$convergence, minEigenvalue=min(eigvals))
}

#' Kernel regression to fit a rotation R(x) to (x, R) data, with no pre-processing of the independent variable.
#' 
#' This function interpolates/extrapolates a rotation R for a given x-value, based on a given set of (x, R) data. The chosen bandwidth h may have a substantial effect on the results. See rotBandwidthForKernelRegression. If this function doesn't work as you expect, then try rotNativeKernelRegression.
#' @param x A real number. The x-value at which to predict the rotation R.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $r (the rotation matrix R(x)), $error (which should be 0; if it is 1, then increase numSteps), and $minEigenvalue (which should be positive; if not, then the minimization has failed).
rotKernelRegression <- function(x, xs, rs, h, k=dnorm, numSteps=100) {
  x0 <- min(xs)
  x1 <- max(xs)
  rotNativeKernelRegression((x - x0) / (x1 - x0), scales(xs), rs, h, k, numSteps)
}

#' Cross-validation algorithm for choosing the bandwidth for kernel regression.
#' 
#' See Davis et al. (2010).
#' @param xs A vector of real numbers. The x-values at which R(x) is known. Assumed to have length >= 3.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A real number (positive). The bandwidth to use in kernel regression.
rotBandwidthForKernelRegression <- function(xs, rs, k=dnorm, numSteps=100) {
  # Build jackknifed lists ahead of time.
  n <- length(xs)
  xjs <- listOmitting(xs)
  rjs <- listOmitting(rs)
  # Define the function to minimize.
  g <- function(h) {
    rjhs <- lapply(1:n, function(j) rotKernelRegression(xs[[j]], xjs[[j]], rjs[[j]], h, k, numSteps))
    sum(sapply(1:n, function(j) rotDistance(rs[[j]], rjhs[[j]]$r)^2))
  }
  # Minimize the function with respect to h.
  solution <- optimize(g, interval=c(0, pi))
  solution$minimum
}

# What is the point of this? It might be orphaned.
rotKernelRegressionBandwidthMisfit <- function(xs, rs, h, k=dnorm, numSteps=100) {
  n <- length(xs)
  xjs <- listOmitting(xs)
  rjs <- listOmitting(rs)
  rjhs <- lapply(1:n, function(j) rotKernelRegression(xs[[j]], xjs[[j]], rjs[[j]], h, k, numSteps))
  sum(sapply(1:n, function(j) rotDistance(rs[[j]], rjhs[[j]]$r)^2))
}

#' R^2 statistic for kernel regression.
#' 
#' As always, make sure the errors are 0 and the minEigenvalues are positive.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list with elements $rSquared, $errors, $minEigenvalues.
rotRsquaredForKernelRegression <- function(xs, rs, h, k=dnorm, numSteps=100) {
  n <- length(xs)
  regrs <- lapply(xs, function(x) rotKernelRegression(x, xs, rs, h, k, numSteps))
  e <- sum(sapply(1:n, function(i) rotDistance(rs[[i]], regrs[[i]]$r)^2)) / (2 * n)
  rSquared <- 1 - e / rotVariance(rs, rotFrechetMean(rs))
  errors <- sapply(regrs, function(regr) regr$error)
  minEigenvalues <- sapply(regrs, function(regr) regr$minEigenvalue)
  list(rSquared=rSquared, errors=errors, minEigenvalues=minEigenvalues)
}

#' Permutation test for significance of a kernel regression.
#' 
#' Returns up to numPerms R^2 values. May be fewer than numPerms, because of some regressions failing. Let n be the dimension of this vector and g the number of R^2 values greater than the R^2 for the original regression of the data. Let p = g / n. Small values of p indicate that the dependency detected by the regression is meaningful.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param numPerms A real number (positive integer). The number of permutations to try.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A vector of length up to numPerms, containing R^2 values for successful permutations.
rotKernelRegressionPermutations <- function(xs, rs, h, numPerms, k=dnorm, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    rSquaredErrorsMins <- rotRsquaredForKernelRegression(sample(ys, size=length(ys)), rs, h, numSteps=numSteps)
    error <- max(abs(rSquaredErrorsMins$errors))
    minEigenvalue <- min(rSquaredErrorsMins$minEigenvalues)
    c(error, minEigenvalue, rSquaredErrorsMins$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}



### PLOTTING ###

#' x-z-x Euler angle plot of the space of rotations.
#'
#' The plot represents each rotation as a point in Euler angle space, which is a (2 pi) x (pi) x (2 pi) box. This plot does not have an equal-angle or equal-volume property. Indeed, it is extremely distorted near certain gimbal lock lines. Curves and triangles are not supported.
#' @param points A list of rotation matrices.
#' @param showBoundary Logical. Whether to show the bounding box.
#' @param ... Additional plotting options to be passed to plot3D.
#' @return NULL.
rotEulerAnglePlot <- function(points, showBoundary=FALSE, ...) {
  xzxs <- lapply(points, rotXZXAnglesFromMatrix)
  if (showBoundary)
    curves <- list(
      list(c(pi, pi, pi), c(-pi, pi, pi), c(-pi, 0, pi), c(pi, 0, pi), c(pi, pi, pi)),
      list(c(pi, pi, -pi), c(-pi, pi, -pi), c(-pi, 0, -pi), c(pi, 0, -pi), c(pi, pi, -pi)),
      list(c(pi, pi, pi), c(pi, pi, -pi)),
      list(c(-pi, pi, pi), c(-pi, pi, -pi)),
      list(c(-pi, 0, pi), c(-pi, 0, -pi)),
      list(c(pi, 0, pi), c(pi, 0, -pi)))
  else
    curves <- list()
  plot3D(radius=pi, points=xzxs, curves=curves, ...)
}

#' Breaking a curve into pieces that don't cross the boundary of a ball plot.
#'
#' Probably you don't want to use this. It's used internally in our plotting functions. Anyway, if two successive points are very far apart (more than 1.5 times the radius of the plot), then the curve is broken between them. Two antipodal boundary points are inserted, so that one piece of the curve ends on the boundary, and the next piece of the curve starts at the antipode. To clarify, if the original curve never needs breaking, then the returned list will have length 1.
#' @param vs A list of 3D real vectors.
#' @param radius A real number (positive). The radius of the plot.
#' @param ballFromRotation A function from the set of special orthogonal 3x3 real matrices to the set of 3D real vectors. A right inverse to rotationFromBall.
#' @param rotationFromBall A function from the set of 3D real vectors to the set of special orthogonal 3x3 real matrices. A left inverse to ballFromRotation.
#' @return A list of lists of 3D real vectors.
rotBallCurves <- function(vs, radius, ballFromRotation, rotationFromBall) {
  n <- length(vs)
  if (n <= 1)
    list(vs)
  else {
    curves <- list()
    curve <- list(vs[[1]])
    boundSquared <- (1.5 * radius)^2
    for (i in 1:(n - 1))
      if (sum((vs[[i]] - vs[[i + 1]])^2) <= boundSquared)
        # This segment of the curve does not need to be broken.
        curve[[length(curve) + 1]] <- vs[[i + 1]]
    else {
      # The curve must be broken between the ith and (i + 1)th elements.
      q <- rotationFromBall(vs[[i]])
      r <- rotationFromBall(vs[[i + 1]])
      u <- rotAxisAngleFromMatrix(r %*% t(q))[1:3]
      # Finding the boundary means finding t such that aa t^2 + bb t + cc == 0,
      # where aa, bb, cc have been derived using Mathematica.
      cc <- 1 + q[1, 1] + q[2, 2] + q[3, 3]
      bb <- u[1] * (q[2, 3] - q[3, 2]) + u[2] * (q[3, 1] - q[1, 3]) + u[3] * (q[1, 2] - q[2, 1])
      bb <- 2 * bb
      aa <- u[1] * u[2] * (q[1, 2] + q[2, 1]) - u[3]^2 * (q[1, 1] + q[2, 2])
      aa <- aa + u[1] * u[3] * (q[1, 3] + q[3, 1]) - u[2]^2 * (q[1, 1] + q[3, 3])
      aa <- aa + u[2] * u[3] * (q[2, 3] + q[3, 2]) - u[1]^2 * (q[2, 2] + q[3, 3])
      aa <- 2 * aa + cc
      sols <- realQuadraticSolutions(aa, bb, cc)
      if (length(sols) == 0 || length(sols) == 3) {
        # Degenerate case. Just break the segment here.
        curves[[length(curves) + 1]] <- curve
        curve <- list(vs[[i + 1]])
      } else {
        # Typical case. Find a ball vector at the boundary.
        s <- sols[[1]]
        a <- arcCos((1 - s^2) / (1 + s^2))
        r <- rotMatrixFromAxisAngle(c(u, a)) %*% q
        v <- ballFromRotation(r)
        # Determine how to pair v and -v with vs[[i]] and vs[[i + 1]].
        if (sum((vs[[i]] - v)^2) < sum((vs[[i]] + v)^2))
          w <- -v
        else {
          w <- v
          v <- -v
        }
        # Finish this curve and start the next.
        curve[[length(curve) + 1]] <- v
        curves[[length(curves) + 1]] <- curve
        curve <- list(w, vs[[i + 1]])
      }
    }
    curves[[length(curves) + 1]] <- curve
    curves
  }
}

#' Infrastucture for equal-volume and similar plots.
#'
#' Unless you are designing your own system of plotting, you probably do not want to call this function. Call rotEqualVolumePlot, etc. instead.
#' @param radius A real number (positive). If NULL, then a radius is chosen to contain all points, curves, and triangles.
#' @param points A list of 3D real vectors.
#' @param curves A list of lists of 3D real vectors.
#' @param triangles A list of length-3 lists of 3D real vectors.
#' @param colors A list of strings (colors). Used to color the points only. Colors are recycled, just as in rgl.points and rgl.spheres.
#' @param ballFromRotation A function from the set of rotation matrices to the set of 3D real vectors. A right inverse to rotationFromBall.
#' @param rotationFromBall A function from the set of 3D real vectors to the set of rotation matrices. A left inverse to ballFromRotation.
#' @param boundaryAlpha A real number (between 0 and 1). The opacity of the boundary sphere. A value of 0 turns off the sphere entirely.
#' @param simplePoints A logical. Whether to plot points as points or as spheres. Spheres are higher-quality and slower.
#' @param backgroundColor String (color). Color of background and fog.
#' @param curveColor A string (color). Color of curves.
#' @param curveWidth A real number (positive). Width of curves in pixels.
#' @param fogStyle A string, either "none", "linear", "exp", or "exp2". The style of fog used to suggest depth. See rgl.bg.
#' @param pointSize A real number (positive). The size of the points or spheres used to depict points. If simplePoints, then measured in pixels with default 3. If not simplePoints, then measured in the same units as the radius of the plot, with default 0.02 * radius.
#' @param trianglesRaw Either NULL or a vector of 9 * m real numbers. An alternative way to pass m triangles: All of the x-coordinates of vertices, then all of the ys, then all of the zs.
#' @param axesColors A character vector of length 3. The colors for the axes.
#' @return NULL.
rotBallPlot <- function(radius=NULL, points, curves, triangles, colors, ballFromRotation, rotationFromBall, boundaryAlpha=0.1, simplePoints=FALSE, backgroundColor="black", curveColor="white", curveWidth=1, fogStyle="linear", pointSize=NULL, trianglesRaw=NULL, axesColors=c("red", "green", "blue")) {
  # Initialize the window. The user will close it later.
  rgl.open()
  rgl.bg(color=backgroundColor, fogtype=fogStyle)
  # Figure out the radius and pointSize if necessary.
  if (is.null(radius)) {
    pointRadius <- 0
    curveRadius <- 0
    triangleRadius <- 0
    f <- function(pointList) {max(abs(simplify2array(pointList)))}
    if (length(points) >= 1)
      pointRadius <- f(points)
    if (length(curves) >= 1)
      curveRadius <- max(sapply(curves, f))
    if (length(triangles) >= 1)
      triangleRadius <- max(sapply(triangles, f))
    radius <- max(pointRadius, curveRadius, triangleRadius)
  }
  if (is.null(pointSize)) {
    if (simplePoints)
      pointSize <- 3
    else
      pointSize <- 0.02 * radius
  }
  # Draw the three coordinate axes.
  xs <- radius * c(0, 1, 1, 1, 1, 1, -1, -1, -1, -1)
  ys <- radius * c(0, 0, -0.1, 0.1, 0, 0, -0.1, 0.1, 0, 0)
  zs <- radius * c(0, 0, 0, 0, -0.1, 0.1, 0, 0, -0.1, 0.1)
  if (radius > 1) {
    xs <- c(xs, as.numeric(sapply(1:radius, function(x) c(x, x, x, x))))
    ys <- c(ys, replicate(floor(radius), c(-0.1, 0.1, 0, 0)))
    zs <- c(zs, replicate(floor(radius), c(0, 0, -0.1, 0.1)))
  }
  rgl.lines(x=xs, y=ys, z=zs, color=axesColors[[1]], lwd=3)
  rgl.lines(x=zs, y=xs, z=ys, color=axesColors[[2]], lwd=3)
  rgl.lines(x=ys, y=zs, z=xs, color=axesColors[[3]], lwd=3)
  # Draw the points.
  if (length(points) >= 1) {
    xs <- sapply(points, function(p) {p[1]})
    ys <- sapply(points, function(p) {p[2]})
    zs <- sapply(points, function(p) {p[3]})
    if (simplePoints)
      rgl.points(x=xs, y=ys, z=zs, color=colors, size=pointSize)
    else
      rgl.spheres(x=xs, y=ys, z=zs, radius=pointSize, color=colors, lit=FALSE)
  }
  # Draw the curves.
  if (length(curves) >= 1)
    for (curve in curves) {
      curvesNew <- rotBallCurves(curve, radius, ballFromRotation, rotationFromBall)
      for (cur in curvesNew) {
        xs <- sapply(cur, function(p) {p[1]})
        ys <- sapply(cur, function(p) {p[2]})
        zs <- sapply(cur, function(p) {p[3]})
        rgl.linestrips(x=xs, y=ys, z=zs, color=curveColor, lwd=curveWidth)
      }
    }
  # Draw the triangles.
  if (length(triangles) >= 1) {
    xs <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[1]})})
    ys <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[2]})})
    zs <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[3]})})
    rgl.triangles(x=xs, y=ys, z=zs)#, alpha=1.0, back="cull")
  }
  if (!is.null(trianglesRaw)) {
    numVert <- length(trianglesRaw) / 3
    xs <- trianglesRaw[1:numVert]
    ys <- trianglesRaw[(numVert + 1):(2 * numVert)]
    zs <- trianglesRaw[(2 * numVert + 1):(3 * numVert)]
    rgl.triangles(x=xs, y=ys, z=zs)
  }
  # Draw the container sphere.
  if (boundaryAlpha > 0)
    rgl.spheres(x=0, y=0, z=0, radius=radius, alpha=boundaryAlpha, back="cull")
  NULL
}

#' Axis-angle plot of the space of rotations.
#'
#' The plot represents each rotation as a point, whose direction from the origin indicates the axis of rotation, and whose distance from the origin indicates the angle of rotation. Angles between 0 and pi are used. So SO(3) plots as a ball of radius pi. This plot does not have an equal-angle or equal-volume property. Curves that cross the boundary of the plot are clipped automatically. On the other hand, triangles are not clipped.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotAxisAnglePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  aaFromRotation <- function(r) rotLeftTangentFromMatrix(r, diag(c(1, 1, 1)))
  rotationFromAA <- function(v) rotMatrixFromLeftTangent(v, diag(c(1, 1, 1)))
  pointsNew <- lapply(points, aaFromRotation)
  curvesNew <- lapply(curves, function(curve) lapply(curve, aaFromRotation))
  trianglesNew <- lapply(triangles, function(tri) lapply(tri, aaFromRotation))
  rotBallPlot(pi, pointsNew, curvesNew, trianglesNew, colors, aaFromRotation, rotationFromAA, ...)
}

#' Equal-angle plot of the space of rotations.
#'
#' The plot represents each rotation as a point, whose direction from the origin indicates the axis of rotation, and whose distance from the origin indicates the tangent quarter-angle of rotation. All positive angles are potentially used. So SO(3) plots as a ball of radius 1.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotEqualAnglePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  pointsNew <- lapply(points, rotEqualAngleFromMatrix)
  curvesNew <- lapply(curves, function(curve) lapply(curve, rotEqualAngleFromMatrix))
  trianglesNew <- lapply(triangles, function(tri) lapply(tri, rotEqualAngleFromMatrix))
  rotBallPlot(1, pointsNew, curvesNew, trianglesNew, colors, rotEqualAngleFromMatrix, rotMatrixFromEqualAngle, ...)
}

#' Rodrigues plot of the space of rotations.
#'
#' The plot represents each rotation as a point, whose direction from the origin indicates the axis of rotation, and whose distance from the origin indicates the tangent half-angle of rotation. All positive angles are potentially used. So SO(3) plots as a ball of infinite radius.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotRodriguesPlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  pointsNew <- lapply(points, rotRodriguesFromMatrix)
  curvesNew <- lapply(curves, function(curve) lapply(curve, rotRodriguesFromMatrix))
  trianglesNew <- lapply(triangles, function(tri) lapply(tri, rotRodriguesFromMatrix))
  rotBallPlot(NULL, pointsNew, curvesNew, trianglesNew, colors, rotRodriguesFromMatrix, rotMatrixFromRodrigues, ...)
}

#' Equal-volume plot of the space of rotations, presented in equal-volume coordinates.
#'
#' Probably you don't want to use this. Look at rotEqualVolumePlot first. Each rotation is passed to this function as an equal-volume vector, as produced by rotEqualVolumeFromMatrix, for example. The function simply plots those vectors, in a ball of radius approximately 0.62. The plot is equal-volume; it accurately represents volumes of regions of SO(3). Curves that cross the boundary of the plot are clipped automatically. On the other hand, triangles are not clipped.
#' @param points A list of 3D real vectors.
#' @param curves A list of lists of 3D real vectors.
#' @param triangles A list of length-3 lists of 3D real vectors.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotNativeEqualVolumePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  rotBallPlot(rotEqualVolumeRadius, points, curves, triangles, colors, rotEqualVolumeFromMatrix, rotMatrixFromEqualVolume, ...)
}

#' Equal-volume plot of the space of rotations.
#'
#' Each rotation is passed to this function as a special orthogonal matrix. The function converts to volumetric representation and then simply calls volumetricPlot.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotEqualVolumePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  pointsNew <- lapply(points, rotEqualVolumeFromMatrix)
  curvesNew <- lapply(curves, function(curve) {lapply(curve, rotEqualVolumeFromMatrix)})
  trianglesNew <- lapply(triangles, function(tri) {lapply(tri, rotEqualVolumeFromMatrix)})
  rotNativeEqualVolumePlot(points=pointsNew, curves=curvesNew, triangles=trianglesNew, colors, ...)
}



### HIGHER-LEVEL PLOTS ###

#' Equal-angle plot of geodesic regression results.
#' 
#' The user may choose to draw extra spurs, joining the Rs to their corresponding predictions on the curve.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices.
#' @param regr The output of rotGeodesicRegression with those xs and rs.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param spurs A list of rotation matrices. When not NULL, it's usually equal to Rs.
#' @param ... Other plotting options to be passed to rotEqualAnglePlot.
#' @return NULL.
rotGeodesicRegressionPlot <- function(xs, rs, regr, colors=hues(xs), spurs=NULL, ...) {
  rotA <- rotExp(min(xs) * regr$m) %*% regr$b
  rotB <- rotExp(mean(range(xs)) * regr$m) %*% regr$b
  rotC <- rotExp(max(xs) * regr$m) %*% regr$b
  curves <- list(rotGeodesicPoints(rotA, rotB, 20), rotGeodesicPoints(rotB, rotC, 20))
  if (!is.null(spurs)) {
    spurCurves <- lapply(1:length(xs), function(i) {rotGeodesicPoints(spurs[[i]], regr$prediction(xs[[i]]), 10)})
    curves <- c(curves, spurCurves)
  }
  rotEqualAnglePlot(points=rs, curves=curves, colors=colors, ...)
}

#' Equal-angle plot of Mahalanobis inference results.
#' 
#' @param points A list of rotation matrices.
#' @param center A rotation matrix. Typically the $center from the inference.
#' @param leftCovarInv A 3x3 real matrix (symmetric, positive-definite). Typically the $leftCovarInv from the inference.
#' @param level A real number. Typically $q095^2 from the inference.
#' @param numNonAdapt A real number (non-negative integer). The number of refinements to the sphere that is deformed into the ellipsoid. Incrementing numNonAdapt improves visual quality but increases time and memory requirements by a factor of four.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Additional plotting options to be passed to rotEqualAnglePlot.
#' @return NULL.
rotEllipsoidPlot <- function(points, center, leftCovarInv, level, numNonAdapt=3, colors=c("white"), ...) {
  tris <- rotEllipsoidTriangles(center, leftCovarInv, level, numNonAdapt=numNonAdapt)
  rotEqualAnglePlot(points=points, triangles=tris, colors=colors, ...)
}

# Helper function for rotEqualAreaTickPlot.
rotEqualAreaTick <- function(v, tickSize) {
  r <- sqrt(v[[1]]^2 + v[[2]]^2)
  if (r <= 0 || 1 - r^2 <= 0)
    tick <- c(0, 0, 0)
  else {
    scalar <- sqrt(1 - r^2) / r
    tick <- c(scalar * v[[1]], scalar * v[[2]], v[[3]] / -scalar)
  }
  list(lower(v), rayNormalized(lower(v) + tickSize * tick))
}

#' Equal-area plot, attempting to show all degrees of freedom.
#' 
#' In each rotation, regards the first row as +-down-pole and the third row as +-ray. Depicts each rotation as a great circle (perpendicular to the pole), a point on that great circle (corresponding to the ray), and a tick mark (indicating the direction of the ray, either toward the center if up or away if down). Really intended for representing faults-with-slip in geology, where the second row is then the vorticity of fault slip. But could be useful in more abstract rotations?
#' @param rs A list of rotation matrices.
#' @param numSteps A real number (positive integer). The number of steps used in drawing each great circle.
#' @param tickSize A real number (non-negative). The length of the ticks.
#' @param ... Additional options to be passed to the underlying lineEqualAreaPlot.
#' @return NULL.
rotEqualAreaTickPlot <- function(rs, numSteps=60, tickSize=0.1, ...) {
  greats <- lapply(rs, function(r) rayGreatCircle(r[1,], numSteps))
  hangs <- lapply(rs, function(r) {if (r[[1, 3]] < 0) r[3,] else -r[3,]})
  ticks <- lapply(hangs, rotEqualAreaTick, tickSize)
  lineEqualAreaPlot(points=hangs, curves=c(greats, ticks), ...)
}


