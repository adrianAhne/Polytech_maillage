//
// Created by adrian on 03/03/16.
//

#ifndef DISTANCEMESHFUNCTIONS_H
#define DISTANCEMESHFUNCTIONS_H


// calculate dot product of two three-dimensional vectors
double dotProduct3D(Point a, Point b);

// Check if point is in triangle tria
int checkPointInTriangle(pMesh mesh, pPoint P, pTria tria);

#endif //DISTANCEMESHFUNCTIONS_H
