#ifndef DISTANCE_H
#define DISTANCE_H
#include <assert.h>

// calculate dot product of two three-dimensional vectors
double dotProduct3D(Point a, Point b);

// Check if point is in triangle tria (NOT USED CURRENTLY)
//void checkPointInTriangle(pMesh mesh, pPoint P, pTria tria, double* re, double* dotProduct);

// (NOT USED CURRENTLY)
//Point calculeProjection(pMesh mesh, pTria tria, pPoint point);

/* calculates the distance between a point and a triangle by using barycentric coordinates */
double distPointToTriangle(pMesh mesh, pTria tria, pPoint P0);

/* Calculating the average distance of a point to the triangle */
double averageDistancePTT(pMesh mesh, pTria tria, pPoint P0);

#endif
