


# These rotations describe the relative motions of the Farallon and Pacific plates. They come from Engebretson et al. (1984), as cited by Prentice (1987). Times are in millions of years before present.
rotationFrom2PsiThetaPhiDeg <- function(twoPsiThetaPhiDeg) {
  rho <- sin(twoPsiThetaPhiDeg[1] * degree / 2)
  phi <- twoPsiThetaPhiDeg[2] * degree
  theta <- twoPsiThetaPhiDeg[3] * degree
  v <- cartesianFromSpherical(c(rho, phi, theta))
  q0 <- cos(twoPsiThetaPhiDeg[1] * degree / 2)
  rotMatrixFromQuaternion(c(q0, v[3], v[1], v[2]))
}
farpacTs <- c(0, 5, 9, 17, 28, 37, 43, 48, 56, 61, 66, 71, 85, 119, 127, 135, 145, 165)
farpac2PsiThetaPhiDegs <- list(c(0, 0, 0), c(3.1, -19, 207), c(6.3, -69, 186), c(12.6, -79, 188), c(20, -81, 205), c(32.6, -88, 155), c(43.5, -87, 85), c(51.4, -86, 70), c(61.7, -85, 49), c(65, -86, 32), c(68.1, -86, 352), c(71.3, -86, 13), c(80.3, -85, 315), c(107.4, -78, 285), c(108.8, -75, 275), c(110.5, -71, 269), c(114.2, -68, 269), c(121.7, -63, 270))
farpacRs <- lapply(farpac2PsiThetaPhiDegs, rotationFrom2PsiThetaPhiDeg)


