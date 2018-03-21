


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This file offers a few functions dealing with deformations, mostly homogeneous and steady. Most computations involve velocity and position gradient tensors (VGTs and PGTs, respectively). In the context of homogeneous deformations, PGTs are also called finite deformation tensors. This library also offers dynamics for rigid (Jeffery, 1922) and deformable (Eshelby, 1957; Bilby et al., 1975) ellipsoids in slow viscous flows.



### MISCELLANY ###

#' Position gradient tensor for a given velocity gradient tensor.
#' 
#' Given a VGT L, returns the corresponding PGT F = exp L. This F represents the finite deformation at time t = 1 of a steady, homogeneous progressive deformation that started at t = 0. Requires R package 'expm'.
defExp <- expm

#' Velocity gradient tensor deduced from a position gradient tensor.
#' 
#' Given a PGT F, returns the 'principal' logarithm L = log F. This L is a VGT representing a steady, homogeneous progressive deformation, whose finite deformation from t = 0 to t = 1 is F. Requires R package 'expm'.
defLog <- logm

#' Kinematic vorticity of a steady, homogeneous deformation.
#' 
#' The kinematic vorticity, commonly denoted w_k, is the ratio of rotation to distortion in the deformation. If there is no distortion, then it is undefined. Otherwise, it can take on any non-negative real value. If the VGT has only real eigenvalues, then w_k <= 1.
#' @param vgt A real 3x3 matrix. The VGT of the deformation.
#' @return A real number (non-negative) or Inf.
defKinematicVorticity <- function(vgt) {
  d <- 0.5 * (vgt + t(vgt))
  w <- 0.5 * (vgt - t(vgt))
  dNorm <- sqrt(tr(crossprod(d, d)))
  wNorm <- sqrt(tr(crossprod(w, w)))
  wNorm / dNorm
}

#' Finite strain ellipsoid of a homogeneous deformation.
#' 
#' @param pgt A real 3x3 matrix, with positive determinant. The PGT of the deformation.
#' @param doNormalize Logical. Whether to scale the ellipsoid to have the same volume as the unit sphere, or not.
#' @return An ellipsoid (a list with members $tensor, $vector, $a, $logA, $rotation; see ellipsoid.R).
defFiniteStrainEllipsoid <- function(pgt, doNormalize=TRUE) {
  eig <- eigen(pgt %*% t(pgt), symmetric=TRUE)
  # Form a rotation matrix with the finite strain axes as the rows.
  rot <- t(eig$vectors)
  if (det(rot) < 0)
    rot[3,] <- -rot[3,]
  # The Finger tensor has eigenvalues ai^2.
  logA <- 0.5 * log(eig$values)
  ellEllipsoidFromRotationLogA(rot, logA, doNormalize=doNormalize)
}



### HOMOGENEOUS SIMPLE SHEAR ###

#' Velocity gradient tensor for homogeneous simple shear along the x-z-plane.
#' 
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @return A real 3x3 matrix.
defSimpleVGT <- function(gamma) {
  rbind(c(0, gamma, 0), c(0, 0, 0), c(0, 0, 0))
}

#' Position gradient tensor for homogeneous simple shear along the x-z-plane.
#' 
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @return A real 3x3 matrix.
defSimplePGT <- function(gamma) {
  rbind(c(1, gamma, 0), c(0, 1, 0), c(0, 0, 1))
}



### HOMOGENEOUS MONOCLINIC TRANSPRESSION ###

#' Velocity gradient tensor for homogeneous monoclinic transpression along the x-z-plane.
#' 
#' See Fossen and Tikoff (1993).
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real 3x3 matrix.
defMonoclinicVGT <- function(gamma, logK) {
  matrix(c(0, 0, 0, gamma, logK, 0, 0, 0, -logK), 3, 3)
}

#' Position gradient tensor for homogeneous monoclinic transpression along the x-z-plane.
#' 
#' See Fossen and Tikoff (1993).
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real 3x3 matrix.
defMonoclinicPGT <- function(gamma, logK) {
  if (logK == 0)
    rbind(c(1, gamma, 0), c(0, 1, 0), c(0, 0, 1))
  else {
    k <- exp(logK)
    rbind(c(1, gamma * (k - 1) / logK, 0), c(0, k, 0), c(0, 0, 1 / k))
  }
}

defMonoclinicInversePGT <- function(gamma, logK) {
  if (logK == 0)
    rbind(c(1, -gamma, 0), c(0, 1, 0), c(0, 0, 1))
  else {
    k <- exp(logK)
    rbind(c(1, gamma * (1 - k) / (k * logK), 0), c(0, 1 / k, 0), c(0, 0, k))
  }
}

#' Critical curve for homogeneous monoclinic transpression.
#' 
#' For certain combinations of (gamma, logK), homogeneous monoclinic transpression produces a finite strain ellipsoid that is an oblate spheroid. Thus the long axis direction, which is frequently interpreted as the lineation direction, is undefined. For any logK, this function returns the corresponding gamma >= 0 to make this happen.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real number, non-negative. The gamma corresponding to logK.
defMonoclinicCriticalGammaFromLogK <- function(logK) {
  k <- exp(logK)
  -(k + 1) * sqrt(k^2 + 1) * logK / k
}

#' Angle of oblique convergence for homogeneous monoclinic transpression.
#' 
#' This is the horizontal angle between the movement direction and the shear zone. For dextral transpression, it is between 0 and pi / 2. Together with magnitude, it forms polar coordinates on the space of homogeneous monoclinic transpressions.
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real number. The angle of oblique convergence, in radians.
defMonoclinicAOCFromGammaLogK <- function(gamma, logK) {
  atan2(-logK, gamma)
}

#' Magnitude of homogeneous monoclinic transpression.
#' 
#' Because of the nature of steady homogeneous deformations, the magnitude can be interpreted in various ways. On the one hand, it measures the intensity of a deformation run for a fixed time interval, say from t = 0 to t = 1. On the other hand, it measures the amount of time that a deformation runs, if we fix that deformation's intensity. Anyway, this is a good measure of 'how far away from no deformation at all' the transpression is. Together with angle of oblique convergence, it forms polar coordinates on the space of homogeneous monoclinic transpressions.
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real number. The angle of oblique convergence, in radians.
defMonoclinicMagnitudeFromGammaLogK <- function(gamma, logK) {
  sqrt(gamma^2 + logK^2)
}

#' Gamma parameter for homogeneous monoclinic transpression, from other parameters.
#' 
#' This function and defMonoclinicLogKFromAOCMagnitude are inverse to defMonoclinicCriticalGammaFromLogK and defMonoclinicAOCFromGammaLogK.
#' @param aoc A real number.  The angle of oblique convergence, in radians.
#' @param mag A real number. The magnitude.
#' @return A real number, gamma.
defMonoclinicGammaFromAOCMagnitude <- function(aoc, mag) {
  mag * cos(aoc)
}

#' log k parameter for homogeneous monoclinic transpression, from other parameters.
#' 
#' This function and defMonoclinicGammaFromAOCMagnitude are inverse to defMonoclinicCriticalGammaFromLogK and defMonoclinicAOCFromGammaLogK.
#' @param aoc A real number.  The angle of oblique convergence, in radians.
#' @param mag A real number. The magnitude.
#' @return A real number, logK.
defMonoclinicLogKFromAOCMagnitude <- function(aoc, mag) {
  -mag * sin(aoc)
}



### HOMOGENEOUS TRICLINIC TRANSPRESSION (JONES AND HOLDSWORTH, 1998; LIN ET AL., 1998) ###

#' Velocity gradient tensor for homogeneous triclinic transpression along a vertical, EW-striking shear plane.
#' 
#' This deformation is homogeneous, so there are effectively only two positions: inside the shear zone versus outside. If you just pass the vector u, then the shear zone is assumed to be infinitely thick, so the station in question is inside the shear zone. If you pass any more parameters than u, then pass all of them. They are used to determine whether the station is inside the shear zone at time t or not.
#' @param u A 3-dimensional real vector. The movement vector of the bounding rigid blocks relative to the shear plane.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param xOfT A 3-dimensional real vector. The location of the point, at which we wish to compute the VGT. Because the deformation is homogeneous, all points inside the zone have the same VGT, and all points outside have VGT = 0.
#' @param tt A real number. The time t at which x(t) = xOfT.
#' @return A real 3x3 matrix.
defTriclinicVGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(0, 0, 0))
  else
    cbind(c(0, 0, 0), u, c(0, 0, -u[[2]]))
}

#' Position gradient tensor for homogeneous triclinic transpression along a vertical, EW-striking shear plane.
#' 
#' Equals defExp(defTriclinicVGT(...)), but faster and more robustly. Same parameters as defTriclinicVGT.!!
defTriclinicPGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(1, 1, 1))
  else if (u[[2]] == 0)
    cbind(c(1, 0, 0), c(u[[1]], 1, u[[3]]), c(0, 0, 1))
  else {
    expU2 <- exp(u[[2]])
    rbind(
      c(1, (expU2 - 1) * u[[1]] / u[[2]], 0),
      c(0, expU2, 0),
      c(0, (expU2 - 1 / expU2) * u[[3]] / (2 * u[[2]]), 1 / expU2))
  }
}



### HOMOGENEOUS TRANSPORT TRANSPRESSION ###

#' Velocity gradient tensor for homogeneous 'transport' transpression along a vertical, EW-striking shear plane.
#' 
#' This deformation is homogeneous, so there are effectively only two positions: inside the shear zone versus outside. If you just pass the vector u, then the shear zone is assumed to be infinitely thick, so the station in question is inside the shear zone. If you pass any more parameters than u, then pass all of them. They are used to determine whether the station is inside the shear zone at time t or not.
#' @param u A 3-dimensional real vector. The movement vector of the bounding rigid blocks relative to the shear plane.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param xOfT A 3-dimensional real vector. The location of the point, at which we wish to compute the VGT. Because the deformation is homogeneous, all points inside the zone have the same VGT, and all points outside have VGT = 0.
#' @param tt A real number. The time t at which x(t) = xOfT.
#' @return A real 3x3 matrix.
defTransportVGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(0, 0, 0))
  else
    rbind(
      c(0, u[[1]], 0),
      c(0, u[[2]], 0),
      c(0, 2 * u[[3]], -u[[2]]))
}

#' Position gradient tensor for homogeneous 'transport' transpression along a vertical, EW-striking shear plane.
#' 
#' Equals defExp(defTransportVGT(...)), but more quickly and robustly. The parameters are the same as in defTransportVGT.
#' @param u A 3-dimensional real vector. The movement vector of the bounding rigid blocks relative to the shear plane.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param xOfT A 3-dimensional real vector. The location of the point, at which we wish to compute the VGT. Because the deformation is homogeneous, all points inside the zone have the same VGT, and all points outside have VGT = 0.
#' @param tt A real number. The time t at which x(t) = xOfT.
#' @return A real 3x3 matrix.
defTransportPGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(1, 1, 1))
  else if (u[[2]] == 0)
    cbind(c(1, 0, 0), c(u[[1]], 1, 2 * u[[3]]), c(0, 0, 1))
  else {
    expU2 <- exp(u[[2]])
    rbind(
      c(1, (expU2 - 1) * u[[1]] / u[[2]], 0),
      c(0, expU2, 0),
      c(0, (expU2 - 1 / expU2) * u[[3]] / u[[2]], 1 / expU2))
  }
}



### HETEROGENEOUS TRANSPRESSION (JAEGER, 1962; ROBIN AND CRUDEN, 1994) ###

# Heterogeneous triclinic transpression (Robin and Cruden, 1994). Get xOfT elsewhere first.
defHeterogeneousVGT <- function(u, hOf1, xOfT, tt) {
  hOfT <- hOf1 * exp(u[[2]] * (tt - 1))
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(0, 0, 0))
  else {
    l22 <- 1.5 * u[[2]] * (1 - xOfT[[2]]^2 / hOfT^2)
    l32 <- u[[3]] + 3 * u[[2]] * xOfT[[2]] * xOfT[[3]] / hOfT^2
    rbind(
      c(0, u[[1]], 0),
      c(0, l22, 0),
      c(0, l32, -l22))
  }
}



### INCLINED TRANSPRESSIONS ###

defInclinedRotation <- function(strike, dip) {
  sinStr <- sin(strike)
  cosStr <- cos(strike)
  sinDip <- sin(dip)
  cosDip <- cos(dip)
  cbind(
    c(sinStr, cosStr, 0),
    c(-cosStr * sinDip, sinStr * sinDip, -cosDip),
    c(-cosStr * cosDip, sinStr * cosDip, sinDip))
}

defInclinedVector <- function(dip, aoc, mag) {
  mag * c(cos(aoc), -sin(dip) * sin(aoc), -cos(dip) * sin(aoc))
}

# Global coordinates y are related to local coordinates x by y = R x + o, so x = R^T (y - o).
# That is, o is the origin of the x-coordinate frame, rendered in y-coordinates.
defInclinedXFromY <- function(y, r, o) {
  as.numeric(t(r) %*% (y - o))
}

defInclinedYFromX <- function(x, r, o) {
  as.numeric(r %*% x + o)
}

defInclinedTensorInXFromInY <- function(tensor, r) {
  t(r) %*% tensor %*% r
}

defInclinedTensorInYFromInX <- function(tensor, r) {
  r %*% tensor %*% t(r)
}

#' Convenience function giving global VGT or PGT for inclined version of homogeneous transpressions.
#' 
#' @param strike A real number. The strike of the shear plane, in radians, but otherwise as usual in geology: measured clockwise from north.
#' @param dip A real number. The dip of the shear plane, in radians. We assume that strike and dip obey the right-hand rule.
#' @param aoc A real number. The angle of oblique convergence, in radians. Kind of like in the monoclinic case.
#' @param mag A real number. The magnitude. Kind of like in the monoclinic case.
#' @param func An R function of the same interface as defTriclinicVGT, defTriclinicPGT, defTransportVGT, etc., giving the tensor in local coordinates.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param yOfT A real 3D vector.
#' @param tt A real number. The time t at which y(t) = yOfT.
#' @param origin A real 3D vector. The shear zone origin o in y = R x + o.
#' @return A real 3x3 matrix. The output of func transformed to global coordinates.
defInclinedTensor <- function(strike, dip, aoc, mag, func, hOf1=1, yOfT=c(0, 0, 0), tt=1, origin=c(0, 0, 0)) {
  rot <- defInclinedRotation(strike, dip)
  vec <- defInclinedVector(dip, aoc, mag)
  xOfT <- defInclinedXFromY(yOfT, rot, origin)
  tensor <- func(vec, hOf1=hOf1, xOfT=xOfT, tt=tt)
  defInclinedTensorInYFromInX(tensor, rot)
}



### RIGID ELLIPSOID DYNAMICS ###

# When the ellipsoid is a spheroid, its orientation is ambiguous. However, the dynamical calculations assume that the orientation has been chosen to line up with the eigenvectors of D in a particular way. So this function returns the adjusted orientation. In practice this adjustment is important mainly when using rigid spheroids or when starting a deformable simulation from a spheroid.
defDynamicsCorrection <- function(q, a, d) {
  qNew <- q
  if (a[[1]] == a[[2]]) {
    b <- qNew %*% d %*% t(qNew)
    if (b[[2, 2]] == b[[1, 1]])
      s <- pi / 4
    else
      s <- 0.5 * atan(2 * b[[1, 2]] / (b[[2, 2]] - b[[1, 1]]))
    r <- matrix(c(cos(s), sin(s), 0, -sin(s), cos(s), 0, 0, 0, 1), 3, 3)
    qNew <- r %*% qNew
  }
  if (a[[1]] == a[[3]]) {
    b <- qNew %*% d %*% t(qNew)
    if (b[[1, 1]] == b[[3, 3]])
      s <- pi / 4
    else
      s <- 0.5 * atan(2 * b[[1, 3]] / (b[[1, 1]] - b[[3, 3]]))
    r <- matrix(c(cos(s), 0, -sin(s), 0, 1, 0, sin(s), 0, cos(s)), 3, 3)
    qNew <- r %*% qNew
  }
  if (a[[2]] == a[[3]]) {
    b <- qNew %*% d %*% t(qNew)
    if (b[[3, 3]] == b[[2, 2]])
      s <- pi / 4
    else
      s <- 0.5 * atan(2 * b[[2, 3]] / (b[[3, 3]] - b[[2, 2]]))
    r <- matrix(c(1, 0, 0, 0, cos(s), sin(s), 0, -sin(s), cos(s)), 3, 3)
    qNew <- r %*% qNew
  }
  qNew
}

# The velocity gradient tensor K = V of the Jeffery (1922) rigid ellipsoid.
defJefferyVGT <- function(q, a, d, w) {
  q <- defDynamicsCorrection(q, a, d)
  dTilde <- q %*% d %*% t(q)
  wTilde12 <- (a[[1]]^2 - a[[2]]^2) / (a[[1]]^2 + a[[2]]^2) * dTilde[[1, 2]]
  wTilde13 <- (a[[1]]^2 - a[[3]]^2) / (a[[1]]^2 + a[[3]]^2) * dTilde[[1, 3]]
  wTilde23 <- (a[[2]]^2 - a[[3]]^2) / (a[[2]]^2 + a[[3]]^2) * dTilde[[2, 3]]
  wTilde <- matrix(c(0, -wTilde12, -wTilde13, wTilde12, 0, -wTilde23, wTilde13, wTilde23, 0), 3, 3)
  w - t(q) %*% wTilde %*% q
}

#' Simulating the rotation of a rigid ellipsoid in a slow viscous flow.
#' 
#' See Jeffery (1922). Uses the classical fourth-order Runge-Kutta algorithm.
#' @param q0 A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 0.
#' @param a A real 3-dimensional vector. The semi-axis lengths a1, a2, a3, in order corresponding to the rows of Q.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 1.
defClassicalJeffery <- function(q0, a, l, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(s, q) {-defJefferyVGT(q, a, d, w)}
  rungeKutta(q0, vel, n, rotProjectedMatrix)
}

#' Simulating the rotation of a rigid ellipsoid in a slow viscous flow (Jeffery, 1922).
#' 
#' The same calculation as defClassicalJeffery, but using a left-invariant Lie group Runge-Kutta method instead of a classical Runge-Kutta method. An improved version of what's in Davis et al. (2013). This function may deliver less error (per time spent) than defClassicalJeffery, but I have not done testing to confirm that.
#' @param q0 A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 0.
#' @param a A real 3-dimensional vector. The semi-axis lengths a1, a2, a3, in order corresponding to the rows of Q.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 1.
defLeftJeffery <- function(q0, a, l, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(q) {-defJefferyVGT(q, a, d, w)}
  rungeKuttaLeft(q0, vel, n, exponential=rotExp)
}

defRightJeffery <- function(q0, a, l, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(s, qT) {defJefferyVGT(t(qT), a, d, w)}
  t(rungeKuttaRight(t(q0), vel, n, exponential=rotExp))
}

# Test that they work identically. It's really a test of rungeKuttaLeft vs. rungeKuttaRight.
#q0 <- rotUniform()
#a <- abs(rnorm(3))
#a <- a / prod(a)^(1 / 3)
#l <- replicate(3, rnorm(3))
#l <- l - diag(c(1, 1, 1)) * tr(l) / 3
#defLeftJeffery(q0, a, l, 10)
#defRightJeffery(q0, a, l, 10)



### DEFORMABLE ELLIPSOID DYNAMICS ###

# Computes the Eshelby J-tensors, assuming a1 = a2 = a3.
defSphereEshelbyJ <- function(a1) {
  jTwo <- (4 * pi) / (15 * a1^2) * rbind(c(3, 1, 1), c(1, 3, 1), c(1, 1, 3))
  jOne <- (4 * pi / 3) * c(1, 1, 1)
  list(jTwo, jOne)
}

# Computes the Eshelby J-tensors, assuming a1 = a2 > a3.
defOblateEshelbyJ <- function(a1, a3) {
  j1 <- 2 * pi * a1^2 * a3 / (a1^2 - a3^2)^(3 / 2) * (acos(a3 / a1) - a3 / a1 * (1 - a3^2 / a1^2)^(1 / 2))
  j2 <- j1
  j3 <- 4 * pi - j1 - j2
  j13 <- (j3 - j1) / (3 * (a1^2 - a3^2))
  j12 <- pi / (3 * a1^2) - j13 / 4
  j11 <- 3 * j12
  j23 <- (j3 - j2) / (3 * (a1^2 - a3^2))
  j22 <- 4 * pi / (3 * a1^2) - j12 - j23
  j33 <- 4 * pi / (3 * a3^2) - j13 - j23
  jTwo <- rbind(c(j11, j12, j13), c(j12, j22, j23), c(j13, j23, j33))
  jOne <- c(j1, j2, j3)
  list(jTwo, jOne)
}

# Computes the Eshelby J-tensors, assuming a1 > a2 = a3.
defProlateEshelbyJ <- function(a1, a2) {
  j3 <- 2 * pi * a1 * a2^2 / (a1^2 - a2^2)^(3 / 2) * (a1 / a2 * (a1^2 / a2^2 - 1)^(1 / 2) - acosh(a1 / a2))
  j2 <- j3
  j1 <- 4 * pi - j2 - j3
  j13 <- (j3 - j1) / (3 * (a1^2 - a2^2))
  j12 <- (j2 - j1) / (3 * (a1^2 - a2^2))
  j11 <- 4 * pi / (3 * a1^2) - j12 - j13
  j23 <- pi / (3 * a2^2) - j12 / 4
  j22 <- 3 * j23
  j33 <- 4 * pi / (3 * a2^2) - j13 - j23
  jTwo <- rbind(c(j11, j12, j13), c(j12, j22, j23), c(j13, j23, j33))
  jOne <- c(j1, j2, j3)
  list(jTwo, jOne)
}

# Helper function for defTypicalEshelbyJOne. Not my best work.
defTypicalEshelbyIntegral <- function(ai, a1, a2, a3) {
  f <- function(u) {
    1 / ((ai^2 + u) * ((a1^2 + u) * (a2^2 + u) * (a3^2 + u))^(1 / 2))
  }
  tryCatch(
    integrate(f, 0, Inf)$value,
    error=function(e) tryCatch(
      integrate(f, 0, 1)$value + integrate(f, 1, Inf)$value,
      error=function(e) tryCatch(
        integrate(f, 0, 1)$value + integrate(f, 1, 10) + integrate(f, 10, Inf)$value,
        error=function(e) {print(paste("error: defTypicalEshelbyIntegral failed on", ai, a1, a2, a3)); NaN})))
}

# Computes the Eshelby J-tensor [J_i], assuming a1 != a2 != a3 != a1. Does not assume a1 > a2 > a3.
defTypicalEshelbyJOne <- function(a) {
  # Set a1 < a2 < a3.
  a1 <- min(a)
  a3 <- max(a)
  a2 <- sum(a) - a1 - a3
  # The integrals for the two largest ai are most likely to be stable?
  j3 <- 2 * pi * a1 * a2 * a3 * defTypicalEshelbyIntegral(a3, a1, a2, a3)
  j2 <- 2 * pi * a1 * a2 * a3 * defTypicalEshelbyIntegral(a2, a1, a2, a3)
  j1 <- 4 * pi - j3 - j2
  # Permute back.
  c(j1, j2, j3)[order(order(a))]
}

# Computes the Eshelby J-tensor [J_ij], assuming a1 != a2 != a3 != a1. Does not assume a1 > a2 > a3.
defTypicalEshelbyJTwo <- function(a, jOne) {
  j13 <- (jOne[[3]] - jOne[[1]]) / (3 * (a[[1]]^2 - a[[3]]^2))
  j12 <- (jOne[[2]] - jOne[[1]]) / (3 * (a[[1]]^2 - a[[2]]^2))
  j11 <- 4 * pi / (3 * a[[1]]^2) - j12 - j13
  j23 <- (jOne[[3]] - jOne[[2]]) / (3 * (a[[2]]^2 - a[[3]]^2))
  j22 <- 4 * pi / (3 * a[[2]]^2) - j12 - j23
  j33 <- 4 * pi / (3 * a[[3]]^2) - j13 - j23
  rbind(c(j11, j12, j13), c(j12, j22, j23), c(j13, j23, j33))
}

# Computes the Eshelby J-tensors, using helper functions for all subcases.
defEshelbyJ <- function(a) {
  if (a[[1]] != a[[2]] && a[[2]] != a[[3]] && a[[3]] != a[[1]]) {
    jOne <- defTypicalEshelbyJOne(a)
    jTwo <- defTypicalEshelbyJTwo(a, jOne)
    list(jTwo, jOne)
  } else if (a[[1]] == a[[2]] && a[[2]] == a[[3]])
    defSphereEshelbyJ(a[[1]])
  else if (a[[1]] == a[[2]] && a[[2]] > a[[3]])
    defOblateEshelbyJ(a[[1]], a[[3]])
  else if (a[[1]] > a[[2]] && a[[2]] == a[[3]])
    defProlateEshelbyJ(a[[1]], a[[3]])
  else if (a[[2]] == a[[3]] && a[[3]] > a[[1]]) {
    js <- defOblateEshelbyJ(a[[2]], a[[1]])
    p <- rbind(c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  } else if (a[[2]] > a[[3]] && a[[3]] == a[[1]]) {
    js <- defProlateEshelbyJ(a[[2]], a[[1]])
    p <- rbind(c(0, 1, 0), c(1, 0, 0), c(0, 0, 1))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  } else if (a[[3]] == a[[1]] && a[[1]] > a[[2]]) {
    js <- defOblateEshelbyJ(a[[3]], a[[2]])
    p <- rbind(c(1, 0, 0), c(0, 0, 1), c(0, 1, 0))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  } else { # a[[3]] > a[[1]] && a[[1]] == a[[2]]
    js <- defProlateEshelbyJ(a[[3]], a[[1]])
    p <- rbind(c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  }
}

# Returns the velocity gradient tensor K of the Eshelby deformable ellipsoid.
defEshelbyVGT <- function(q, a, d, w, r) {
  js <- defEshelbyJ(a)
  jTwo <- js[[1]]
  jOne <- js[[2]]
  # Get the diagonal entries of C-tilde into cc.
  dTilde <- q %*% d %*% t(q)
  aa <- diag(c(1, 1, 1)) + (3 / (4 * pi)) * (r - 1) * jTwo %*% diag(a^2)
  bb <- diag(dTilde)
  cc <- solve(aa, bb)
  # Compute the other entries of C-tilde.
  cc12 <- dTilde[[1, 2]] / (1 + (r - 1) * (3 / (4 * pi)) * (a[[1]]^2 + a[[2]]^2) * jTwo[[1, 2]])
  cc13 <- dTilde[[1, 3]] / (1 + (r - 1) * (3 / (4 * pi)) * (a[[1]]^2 + a[[3]]^2) * jTwo[[1, 3]])
  cc23 <- dTilde[[2, 3]] / (1 + (r - 1) * (3 / (4 * pi)) * (a[[2]]^2 + a[[3]]^2) * jTwo[[2, 3]])
  cTilde <- rbind(c(cc[[1]], cc12, cc13), c(cc12, cc[[2]], cc23), c(cc13, cc23, cc[[3]]))
  # Compute K-tilde - W-tilde without computing either of those.
  kTildeMinusWTildeMinusCTilde12 <- (r - 1) / (4 * pi) * (jOne[[1]] - jOne[[2]]) * cTilde[[1, 2]]
  kTildeMinusWTildeMinusCTilde13 <- (r - 1) / (4 * pi) * (jOne[[1]] - jOne[[3]]) * cTilde[[1, 3]]
  kTildeMinusWTildeMinusCTilde23 <- (r - 1) / (4 * pi) * (jOne[[2]] - jOne[[3]]) * cTilde[[2, 3]]
  kTildeMinusWTilde <- cTilde + rbind(
    c(0, kTildeMinusWTildeMinusCTilde12, kTildeMinusWTildeMinusCTilde13),
    c(-kTildeMinusWTildeMinusCTilde12, 0, kTildeMinusWTildeMinusCTilde23),
    c(-kTildeMinusWTildeMinusCTilde13, -kTildeMinusWTildeMinusCTilde23, 0))
  # Return K = Q^T (K-tilde - W-tilde) Q + W.
  t(q) %*% kTildeMinusWTilde %*% q + w
}

#' Simulating the deformation of a deformable ellipsoid in a slow viscous flow.
#' 
#' See Eshelby, (1957); Bilby et al. (1975). Uses the classical fourth-order Runge-Kutta algorithm, based on the differential equation E-dot = -K^T E - E K (Davis et al., 2013).
#' @param e0 A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 0.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param r A real number (positive). The viscosity ratio between the clast and the host rock. For example, a competent clast might have r = 10.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 1.
defClassicalEshelby <- function(e0, l, r, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  eDot <- function(s, e) {
    rotA <- ellRotationAFromTensor(e)
    q <- defDynamicsCorrection(rotA$rotation, rotA$a, d)
    k <- defEshelbyVGT(q, rotA$a, d, w, r)
    -t(k) %*% e - e %*% k
  }
  rungeKutta(e0, eDot, n)
}

#' Simulating the deformation of a deformable ellipsoid in a slow viscous flow.
#' 
#' See Eshelby, (1957); Bilby et al. (1975). Similar to defClassicalEshelby, but uses the left-invariant Lie group, fourth-order Runge-Kutta method of Davis et al. (2013). Should be faster (more precise per time required). Requires R package 'expm'.
#' @param e0 A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 0.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param r A real number (positive). The viscosity ratio between the clast and the host rock. For example, a competent clast might have r = 10.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 1.
defLeftEshelby <- function(e0, l, r, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(fInv) {
    rotA <- ellRotationAFromTensor(t(fInv) %*% e0 %*% fInv)
    q <- defDynamicsCorrection(rotA$rotation, rotA$a, d)
    k <- defEshelbyVGT(q, rotA$a, d, w, r)
    -k
  }
  fnInv <- rungeKuttaLeft(diag(c(1, 1, 1)), vel, n)
  t(fnInv) %*% e0 %*% fnInv
}

# Tests.
# ellClassicalEshelby(diag(c(1, 1, 1)), rbind(c(1, 5, -2), c(1, 1, 0), c(0, -3, -1)), 10, 4)
# ellLeftEshelby(diag(c(1, 1, 1)), rbind(c(1, 5, -2), c(1, 1, 0), c(0, -3, -1)), 10, 4)


