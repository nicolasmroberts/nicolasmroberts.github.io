


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



/// PLOTTING FUNCTIONS ///
  
  // Returns the midpoint between the two volumetric vectors, via a trip to SO(3).
// If both vectors are at the boundary of the plot, then so is the output.
void volumetricMidpoint(double *vol1, double *vol2, double *out) {
  // double radius = pow((3.0 / (4.0 * M_PI)), 1.0 / 3.0);
  double radius = 0.6203504908994;
  double diff1, diff2;
  diff1 = fabs(scalarProduct(vol1, vol1) - radius * radius);
  diff2 = fabs(scalarProduct(vol1, vol1) - radius * radius);
  if (diff1 < 0.00000001 && diff2 < 0.00000001) {
    // The rotations are at the boundary of the equal-volume plot.
    vectorAdd(vol1, vol2, out);
    vectorScale(radius / sqrt(scalarProduct(out, out)), out, out);
  } else {
    // The rotations are not at the boundary; just average them.
    double rot1[9], rot2[9], rot[9];
    rotationFromVolumetric(vol1, rot1);
    rotationFromVolumetric(vol2, rot2);
    matrixTranspose(rot1, rot);
    matrixMatrixMultiply(rot2, rot, rot);
    double angle, axis[3];
    angleAxisFromRotation(rot, &angle, axis);
    rotationFromAngleAxis(0.5 * angle, axis, rot);
    matrixMatrixMultiply(rot, rot1, rot);
    volumetricFromRotation(rot, out);
  }
}

// In the volumetric level surface functions:
  // * A vertex is a 4-tuple of doubles: x, y, z, f.
// * A box is a list of eight vertices, so 32 doubles.
// * boxes is a list of n boxes, so 32 * n doubles.
// * Define f like "double f(double *vol, void *extra)", where vol is a 3-vector
//   and extra is extra data to be used in the furnction.
// * Pass f like "volumetricLevelBase(f, extra, boxes)".

void volumetricLevelVertexCopy(double *v, double *out) {
  out[0] = v[0];
  out[1] = v[1];
  out[2] = v[2];
  out[3] = v[3];
}

void volumetricLevelBoxSet(double *v1, double *v2, double *v3, double *v4,
                           double *v5, double *v6, double *v7, double *v8, double *buffer) {
  volumetricLevelVertexCopy(v1, &(buffer[0]));
  volumetricLevelVertexCopy(v2, &(buffer[4]));
  volumetricLevelVertexCopy(v3, &(buffer[8]));
  volumetricLevelVertexCopy(v4, &(buffer[12]));
  volumetricLevelVertexCopy(v5, &(buffer[16]));
  volumetricLevelVertexCopy(v6, &(buffer[20]));
  volumetricLevelVertexCopy(v7, &(buffer[24]));
  volumetricLevelVertexCopy(v8, &(buffer[28]));
}

// Fills boxes with the unrefined base mesh of seven boxes.
void volumetricLevelBase(double (*f)(double *, void *), void *extra, double *boxes) {
  // Build the eight mesh points on the bounding sphere.
  double pnp[4], ppp[4], nnn[4], nnp[4], npn[4], npp[4], pnn[4], ppn[4];
  vectorSet(0.3581595229189165, 0.3581595229189165, 0.3581595229189165, ppp);
  ppp[3] = (*f)(ppp, extra);
  vectorSet(-ppp[0], -ppp[1], -ppp[2], nnn);
  nnn[3] = (*f)(nnn, extra);
  vectorSet(-ppp[0], -ppp[1], ppp[2], nnp);
  nnp[3] = (*f)(nnp, extra);
  vectorSet(-ppp[0], ppp[1], -ppp[2], npn);
  npn[3] = (*f)(npn, extra);
  vectorSet(-ppp[0], ppp[1], ppp[2], npp);
  npp[3] = (*f)(npp, extra);
  vectorSet(ppp[0], -ppp[1], -ppp[2], pnn);
  pnn[3] = (*f)(pnn, extra);
  vectorSet(ppp[0], -ppp[1], ppp[2], pnp);
  pnp[3] = (*f)(pnp, extra);
  vectorSet(ppp[0], ppp[1], -ppp[2], ppn);
  ppn[3] = (*f)(ppn, extra);
  // Build their halves.
  double pnp2[4], ppp2[4], nnn2[4], nnp2[4], npn2[4], npp2[4], pnn2[4], ppn2[4];
  vectorScale(0.5, nnn, nnn2);
  nnn2[3] = (*f)(nnn2, extra);
  vectorScale(0.5, nnp, nnp2);
  nnp2[3] = (*f)(nnp2, extra);
  vectorScale(0.5, npn, npn2);
  npn2[3] = (*f)(npn2, extra);
  vectorScale(0.5, npp, npp2);
  npp2[3] = (*f)(npp2, extra);
  vectorScale(0.5, pnn, pnn2);
  pnn2[3] = (*f)(pnn2, extra);
  vectorScale(0.5, pnp, pnp2);
  pnp2[3] = (*f)(pnp2, extra);
  vectorScale(0.5, ppn, ppn2);
  ppn2[3] = (*f)(ppn2, extra);
  vectorScale(0.5, ppp, ppp2);
  ppp2[3] = (*f)(ppp2, extra);
  // Build the upward box.
  volumetricLevelBoxSet(nnp2, npp2, ppp2, pnp2, nnp, npp, ppp, pnp, &(boxes[0 * 32]));
  // Build the northward box.
  volumetricLevelBoxSet(ppn2, npn2, npp2, ppp2, ppn, npn, npp, ppp, &(boxes[1 * 32]));
  // Build the eastward box.
  volumetricLevelBoxSet(pnn2, ppn2, ppp2, pnp2, pnn, ppn, ppp, pnp, &(boxes[2 * 32]));
  // Build the southward box.
  volumetricLevelBoxSet(nnn2, pnn2, pnp2, nnp2, nnn, pnn, pnp, nnp, &(boxes[3 * 32]));
  // Build the westward box.
  volumetricLevelBoxSet(npn2, nnn2, nnp2, npp2, npn, nnn, nnp, npp, &(boxes[4 * 32]));
  // Build the downward box.
  volumetricLevelBoxSet(ppn2, pnn2, nnn2, npn2, ppn, pnn, nnn, npn, &(boxes[5 * 32]));
  // Build the central box.
  volumetricLevelBoxSet(nnn2, pnn2, ppn2, npn2, nnp2, pnp2, ppp2, npp2, &(boxes[6 * 32]));
}

// Returns midpoints of the edges, faces, center of the box.
void volumetricLevelEdgesFacesCenter(double (*f)(double *, void *), void *extra, 
    double *box, double *edges, double *faces, double *center) {
  int i;
  // Find the center of each bottom edge.
  for (i = 0; i < 4; i += 1) {
    volumetricMidpoint(
      &(box[i * 4]),
      &(box[((i + 1) % 4) * 4]),
      &(edges[i * 4]));
    edges[i * 4 + 3] = (*f)(&(edges[i * 4]), extra);
  }
  // Find the center of each top edge.
  for (i = 0; i < 4; i += 1) {
    volumetricMidpoint(
      &(box[(i + 4) * 4]),
      &(box[(((i + 1) % 4) + 4) * 4]),
      &(edges[(i + 4) * 4]));
    edges[(i + 4) * 4 + 3] = (*f)(&(edges[(i + 4) * 4]), extra);
  }
  // Find the center of each side edge.
  for (i = 0; i < 4; i += 1) {
    volumetricMidpoint(
      &(box[i * 4]),
      &(box[(i + 4) * 4]),
      &(edges[(i + 8) * 4]));
    edges[(i + 8) * 4 + 3] = (*f)(&(edges[(i + 8) * 4]), extra);
  }
  // Find the center of each face.
  volumetricMidpoint(&(edges[0 * 4]), &(edges[2 * 4]), &(faces[0 * 4]));
  faces[0 * 4 + 3] = (*f)(&(faces[0 * 4]), extra);
  volumetricMidpoint(&(edges[4 * 4]), &(edges[6 * 4]), &(faces[1 * 4]));
  faces[1 * 4 + 3] = (*f)(&(faces[1 * 4]), extra);
  for (i = 0; i < 4; i += 1) {
    volumetricMidpoint(
      &(edges[i * 4]),
      &(edges[(i + 4) * 4]),
      &(faces[(i + 2) * 4]));
    faces[(i + 2) * 4 + 3] = (*f)(&(faces[(i + 2) * 4]), extra);
  }
  // Find the center of the box.
  volumetricMidpoint(&(faces[0 * 4]), &(faces[1 * 4]), center);
  center[3] = (*f)(center, extra);
}

// Returns 1 if box crosses level, 0 if not.
int volumetricLevelBoxCrosses(double *box, double level) {
  int i, lesser = 0, greater = 0;
  double value;
  for (i = 0; i < 8; i += 1) {
    value = box[i * 4 + 3];
    if (value >= level)
      greater = 1;
    if (value <= level)
      lesser = 1;
  }
  return (lesser && greater);
}

// Partitions the given box into eight subboxes, no questions asked.
void volumetricLevelEight(double *box, double *edges, double *faces, double *center, double *boxes) {
  // Construct the southwestdown subbox.
  volumetricLevelBoxSet(&(box[0 * 4]), &(edges[0 * 4]), &(faces[0 * 4]),&(edges[3 * 4]), 
    &(edges[8 * 4]), &(faces[2 * 4]), center, &(faces[5 * 4]), &(boxes[0 * 32]));
  // Construct the southeastdown subbox.
  volumetricLevelBoxSet(&(edges[0 * 4]), &(box[1 * 4]), &(edges[1 * 4]), &(faces[0 * 4]), 
    &(faces[2 * 4]), &(edges[9 * 4]), &(faces[3 * 4]), center, &(boxes[1 * 32]));
  // Construct the northeastdown subbox.
  volumetricLevelBoxSet(&(faces[0 * 4]), &(edges[1 * 4]), &(box[2 * 4]), &(edges[2 * 4]), 
    center, &(faces[3 * 4]), &(edges[10 * 4]), &(faces[4 * 4]), &(boxes[2 * 32]));
  // Construct the northwestdown subbox.
  volumetricLevelBoxSet(&(edges[3 * 4]), &(faces[0 * 4]), &(edges[2 * 4]), &(box[3 * 4]), 
    &(faces[5 * 4]), center, &(faces[4 * 4]), &(edges[11 * 4]), &(boxes[3 * 32]));
  // Construct the southwestup subbox.
  volumetricLevelBoxSet(&(edges[8 * 4]), &(faces[2 * 4]), center, &(faces[5 * 4]), 
    &(box[4 * 4]), &(edges[4 * 4]), &(faces[1 * 4]), &(edges[7 * 4]), &(boxes[4 * 32]));
  // Construct the southeastup subbox.
  volumetricLevelBoxSet(&(faces[2 * 4]), &(edges[9 * 4]), &(faces[3 * 4]), center, 
    &(edges[4 * 4]), &(box[5 * 4]), &(edges[5 * 4]), &(faces[1 * 4]), &(boxes[5 * 32]));
  // Construct the northeastup subbox.
  volumetricLevelBoxSet(center, &(faces[3 * 4]), &(edges[10 * 4]), &(faces[4 * 4]), 
    &(faces[1 * 4]), &(edges[5 * 4]), &(box[6 * 4]), &(edges[6 * 4]), &(boxes[6 * 32]));
  // Construct the northwestup subbox.
  volumetricLevelBoxSet(&(faces[5 * 4]), center, &(faces[4 * 4]), &(edges[11 * 4]), 
    &(edges[7 * 4]), &(faces[1 * 4]), &(edges[6 * 4]), &(box[7 * 4]), &(boxes[7 * 32]));
}

// Computes the point, on the closed line segment from a to b, 
// where the given level surface hits. Returns 1 if that point exists, 0 if not.
int volumetricLevelCut(double *a, double *b, double level, double *out) {
  if (a[3] == b[3])
    if (level == a[3]) {
      vectorCopy(a, out);
      return 1;
    } else
      return 0;
  else {
    double t = (level - a[3]) / (b[3] - a[3]);
    if (t < 0.0 || t > 1.0)
      return 0;
    else {
      // out = a + t (b - a).
      vectorSubtract(b, a, out);
      vectorScale(t, out, out);
      vectorAdd(a, out, out);
      return 1;
    }
  }
}

// Given a box and a level surface, and a buffer of 12 output points (36 
// doubles). Intersects the level surface with the box, computing an 
// approximate polygon (not necessarily planar). Fills the output with the 
// vertices of the polygon, in no particular order, and returns the 
// number of those points.
int volumetricLevelCuts(double *box, double level, double *points) {
  // Find every point where the level surface cuts the 12 edges.
  double cuts[36];
  int i, c = 0;
  int edges[] = {0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7};
  for (i = 0; i < 12; i += 1)
    if (volumetricLevelCut(&(box[edges[2 * i] * 4]), &(box[edges[2 * i + 1] * 4]), level, &(cuts[c * 3])))
      c += 1;
  // Remove duplicates.
  double diff[3];
  int cNew = 0, found, j;
  for (i = 0; i < c; i += 1) {
    found = 0;
    j = 0;
    while (!found && j < cNew) {
      vectorSubtract(&(points[j * 3]), &(cuts[i * 3]), diff);
      found = (scalarProduct(diff, diff) < 0.00000001);
      j += 1;
    }
    if (!found) {
      vectorCopy(&(cuts[i * 3]), &(points[cNew * 3]));
      cNew += 1;
    }
  }
  if (cNew < 3)
    return 0;
  else
    return cNew;
}

// Helper function for volumetricLevelPolygon. Returns whether point b is 
// counterclockwise from point a on the unit circle.
int volumetricLevelWins(double cosa, double sina, double cosb, double sinb) {
  return (sina > 0.0 && sinb <= 0.0) || (sina == 0.0 && sinb < 0.0) ||
    (sina > 0.0 && cosa > cosb) || (sina < 0.0 && sinb < 0.0 && cosa < cosb);
}

// Given unordered points from volumetricLevelCuts, puts them into either
// clockwise or counterclockwise order.
// In theory, could fail if first three cuts are colinear (two ends of an
// edge and somewhere in the middle of that edge). Consider improving.
void volumetricLevelPolygon(int n, double *cuts) {
  // Find a vector normal to the plane containing the first three cuts.
  double u[3], v[3], pole[3];
  vectorSubtract(&(cuts[0 * 3]), &(cuts[1 * 3]), u);
  vectorSubtract(&(cuts[2 * 3]), &(cuts[1 * 3]), v);
  vectorCross(u, v, pole);
  vectorNormalize(pole, pole);
  // Compute the direction from the mean to the first cut.
  double mean[3], sines[12], cosines[12], d0[3], d[3];
  vectorMean(n, cuts, mean);
  vectorSubtract(&(cuts[0 * 3]), mean, d0);
  vectorNormalize(d0, d0);
  //sines[0] = 0.0;
  //cosines[0] = 1.0;
  // Compute the relative direction from the mean to each other cut.
  int i;
  for (i = 1; i < n; i += 1) {
    vectorSubtract(&(cuts[i * 3]), mean, d);
    vectorNormalize(d, d);
    vectorCross(d0, d, u);
    sines[i] = scalarProduct(u, pole);
    cosines[i] = scalarProduct(d0, d);
  }
  // Selection-sort the cuts CCW about the pole, starting from the first cut.
  double swap[3];
  int p, iBest;
  for (p = 1; p < n; p += 1) {
    iBest = p;
    for (i = p + 1; i < n; i += 1)
      if (volumetricLevelWins(cosines[i], sines[i], cosines[iBest], sines[iBest]))
        iBest = i;
    // Swap that unsorted cut into the next sorted spot.
    vectorCopy(&(cuts[p * 3]), swap);
    vectorCopy(&(cuts[iBest * 3]), &(cuts[p * 3]));
    vectorCopy(swap, &(cuts[iBest * 3]));
    // Let the sine and cosine of the unsorted swappee follow.
    sines[iBest] = sines[p];
    cosines[iBest] = cosines[p];
  }
}

// Allocates memory for use in volumetricLevelSurface, equivalent to:
// double aux[7 * 8^depth * 32], out[7 * 8^depth * 32], poly[7 * 8^depth * 36];
// int ns[7 * 8^depth];
// where depth == nonAdapt + adapt.
// Upon failure, *aux is NULL. Upon success, *aux, but none of the other three
// pointers, must later be deallocated using free().
void volumetricLevelMalloc(int nonAdapt, int adapt, double **aux, double **out, int **ns, double **polygons) {
  // Compute 8^depth.
  int i, power = 1;
  for (i = 0; i < (nonAdapt + adapt); i += 1)
    power *= 8;
  // Try to allocate the memory.
  int numDoubles = 7 * power * (32 + 32 + 36);
  int numInts = 7 * power;
  *aux = malloc(numDoubles * sizeof(double) + numInts * sizeof(int));
  if (*aux == NULL)
    fprintf(stderr, "error: volumetricLevelMalloc: could not allocate memory; try decreasing nonAdapt or adapt\n");
  else {
    *out = &((*aux)[32 * 7 * power]);
    *polygons = &((*aux)[64 * 7 * power]);
    *ns = (int *)(&((*aux)[100 * 7 * power]));
  }
}

// Given a function f to contour, a level of f, a non-adaptive refinement 
// depth >= 0, and an adaptive refinement depth >= 0. 
// Given buffers as produced by volumetricLevelMalloc. Subdivides volumetric 
// space into boxes and forms 0 or 1 polygon per box. On output, the integer 
// buffer holds the point counts of the polygons (including those with < 3 
// vertices), and the polygons buffer holds the points of the polygons, 
// packed. Returns the number of polygons.
int volumetricLevelSurface(double (*f)(double *, void *), void *extra, double level, 
    int nonAdapt, int adapt, double *boxesAux, double *boxesOut, int *ns, double *polygons) {
  // Prepare the ping-ponging buffers.
  double *aux, *out, *swap;
  if ((nonAdapt + adapt) % 2 == 0) {
    aux = boxesAux;
    out = boxesOut;
  } else {
    aux = boxesOut;
    out = boxesAux;
  }
  // Set up the base mesh of seven boxes.
  volumetricLevelBase(f, extra, out);
  // Non-adaptively refine all of the boxes.
  int d, b, numBoxes = 7;
  double edges[12 * 4], faces[6 * 4], center[4];
  for (d = 0; d < nonAdapt; d += 1) {
    for (b = 0; b < numBoxes; b += 1) {
      volumetricLevelEdgesFacesCenter(f, extra, &(out[b * 32]), edges, faces, center);
      volumetricLevelEight(&(out[b * 32]), edges, faces, center, &(aux[b * 8 * 32]));
    }
    numBoxes *= 8;
    swap = aux;
    aux = out;
    out = swap;
  }
  // Adaptively refine only the boxes that the level surface crosses.
  int newNumBoxes;
  for (d = 0; d < adapt; d += 1) {
    newNumBoxes = 0;
    for (b = 0; b < numBoxes; b += 1)
      if (volumetricLevelBoxCrosses(&(out[b * 32]), level)) {
        volumetricLevelEdgesFacesCenter(f, extra, &(out[b * 32]), edges, faces, center);
        volumetricLevelEight(&(out[b * 32]), edges, faces, center, &(aux[newNumBoxes * 32]));
        newNumBoxes += 8;
      }
    numBoxes = newNumBoxes;
    swap = aux;
    aux = out;
    out = swap;
  }
  // Now out == boxesOut contains the boxes; convert each to a polygon.
  int n = 0;
  for (b = 0; b < numBoxes; b += 1) {
    ns[b] = volumetricLevelCuts(&(boxesOut[b * 32]), level, &(polygons[n * 3]));
    if (ns[b] >= 3)//!!2015/04/12
      volumetricLevelPolygon(ns[b], &(polygons[n * 3]));
    n += ns[b];
  }
  return numBoxes;
}



/// HIGHER-LEVEL PLOTTING ///

typedef struct {
  double center[9];
  double covarInv[9];
} rotationEllipsoidExtra;

double rotationEllipsoidFunction(double *vol, void *extra) {
  double *center = ((rotationEllipsoidExtra *)extra)->center;
  double *covarInv = ((rotationEllipsoidExtra *)extra)->covarInv;
  double rot[9], vec[3], covarInvVec[3];
  rotationFromVolumetric(vol, rot);
  tangentVectorFromRotation(rot, center, vec);
  vectorMatrixMultiply(covarInv, vec, covarInvVec);
  return scalarProduct(vec, covarInvVec);
}

int rotationEllipsoidSurface(double *center, double *covarInv, double level, int nonAdapt, 
    int adapt, double *boxesAux, double *boxesOut, int *ns, double *polygons) {
  rotationEllipsoidExtra extra;
  matrixCopy(center, extra.center);
  matrixCopy(covarInv, extra.covarInv);
  return volumetricLevelSurface(rotationEllipsoidFunction, (void *)&extra, level, 
    nonAdapt, adapt, boxesAux, boxesOut, ns, polygons);
}

typedef struct {
  double c[4];
  double r;
  int n;
  double *rots;
} rotationKambExtra;

double rotationKambFunction(double *vol, void *extra) {
  rotationKambExtra *e = (rotationKambExtra *)extra;
  double rot[9], a, density = 0.0;
  int i;
  rotationFromVolumetric(vol, rot);
  for (i = 0; i < e->n; i += 1) {
    a = rotationDistance(rot, &(e->rots[i * 9]));
    if (a < e->r)
      density += e->c[0] + e->c[1] * a + e->c[2] * a * a + e->c[3] * a * a * a;
  }
  return density;
}

// rots is n * 9 doubles, representing n rotation matrices in column-major order.
// But these n rotations are actually n / groupSize orientations, if groupSize > 1.
int rotationKambSurface(int n, double *rots, double k, double mult, int degree, int nonAdapt, 
    int adapt, int groupSize, double *boxesAux, double *boxesOut, int *ns, double *polygons) {
  // Compute the Kamb parameters.
  double r, p, sigma, level;
  r = xMinusSinXSolution(M_PI * k * k / (n + groupSize * k * k), 0.00000001);
  p = (r - sin(r)) / M_PI * groupSize;
  sigma = sqrt((n / groupSize) * p * (1.0 - p));
  level = sigma * mult;
  // Prepare the weighting coefficients.
  rotationKambExtra extra;
  if (degree == 3) {
    extra.c[3] = (r - sin(r)) / (3.0 * r * sin(r) + 6.0 * cos(r) + 0.25 * r * r * r * r - 6.0);
    extra.c[2] = -1.5 * r * extra.c[3];
    extra.c[1] = 0.0;
    extra.c[0] = r * r * r * extra.c[3] * 0.5;
  } else if (degree == 1) {
    extra.c[3] = 0.0;
    extra.c[2] = 0.0;
    extra.c[1] = (sin(r) - r) / (cos(r) + 0.5 * r * r - 1.0);
    extra.c[0] = -r * extra.c[1];
  } else {
    extra.c[3] = 0.0;
    extra.c[2] = 0.0;
    extra.c[1] = 0.0;
    extra.c[0] = 1.0;
  }
  // Generate the polygons.
  extra.n = n;
  extra.rots = rots;
  extra.r = r;
  return volumetricLevelSurface(rotationKambFunction, (void *)&extra, level, 
    nonAdapt, adapt, boxesAux, boxesOut, ns, polygons);
}

// Prints the polygons produced by rotationKambSurface to the given file, in
// Mathematica syntax of a list of polygons, where each polygon is a list of
// volumetric points: {{{x, y, z}, ...}, ..., {{x, y, z}, ...}}.
void rotationKambPrintMathematica(FILE *write, int numBoxes, int *ns, double *poly) {
  // Find the first nonempty polygon.
  int b = 0, i, c = 0;
  while (b < numBoxes && ns[b] == 0)
    b += 1;
  if (b < numBoxes) {
    // Print the first polygon.
    fprintf(write, "{{{%lf, %lf, %lf}",
    poly[(c + 0) * 3], poly[(c + 0) * 3 + 1], poly[(c + 0) * 3 + 2]);
    for (i = 1; i < ns[b]; i += 1)
      fprintf(write, ", {%lf, %lf, %lf}",
    poly[(c + i) * 3], poly[(c + i) * 3 + 1], poly[(c + i) * 3 + 2]);
    c += ns[b];
    // Print the rest of the polygons.
    b += 1;
    while (b < numBoxes) {
      if (ns[b] > 0) {
        // Print the end of the previous polygon and the start of this one.
        fprintf(write, "},\n{{%lf, %lf, %lf}",
        poly[(c + 0) * 3], poly[(c + 0) * 3 + 1], poly[(c + 0) * 3 + 2]);
        for (i = 1; i < ns[b]; i += 1) {
          fprintf(write, ", {%lf, %lf, %lf}", poly[(c + i) * 3], poly[(c + i) * 3 + 1], poly[(c + i) * 3 + 2]);
        }
        c += ns[b];
      }
      b += 1;
    }
    // End the last polygon and the list of polygons.
    fprintf(write, "}}");
  }
}


