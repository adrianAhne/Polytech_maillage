//
// Created by adrian on 03/03/16.
//

#ifndef DISTANCEMESHFUNCTIONS_H
#define DISTANCEMESHFUNCTIONS_H
#include <assert.h>

// calculate dot product of two three-dimensional vectors
double dotProduct3D(Point a, Point b);

// Check if point is in triangle tria
void checkPointInTriangle(pMesh mesh, pPoint P, pTria tria, double* re, double* dotProduct);

double distancePointToTriangle(pMesh mesh, pTria tria, pPoint point);

Point calculeProjection(pMesh mesh, pTria tria, pPoint point);

double distPointToTriangle(pMesh mesh, pTria tria, pPoint P0);

#endif
