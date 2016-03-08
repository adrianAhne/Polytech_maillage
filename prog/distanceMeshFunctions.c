//
// Created by adrian on 03/03/16.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "groupfunctions.h"
#include "distanceMeshFunctions.h"



// calculate dot product of two three-dimensional vectors
double dotProduct3D(Point a, Point b)
{
    int sum = 0;
    int i;

    for (i = 0; i < 3; i++) {
        sum += a.c[i] * b.c[i];
    }

    return sum;
}


// Check if point is in triangle tria
int checkPointInTriangle(pMesh mesh, pPoint P, pTria tria)
{
    // points of the triangle
    Point P1, P2, P3;
    P1 = mesh->point[(tria->v[0])];
    P2 = mesh->point[(tria->v[1])];
    P3 = mesh->point[(tria->v[2])];

// Compute vectors
//    v0 = P3 - P1
//    v1 = P2 - P1
//    v2 = P - P1
    Point v0, v1, v2;

    v0.c[0] = P3.c[0] - P1.c[0];
    v0.c[1] = P3.c[1] - P1.c[1];
    v0.c[2] = P3.c[2] - P1.c[2];

    v1.c[0] = P2.c[0] - P1.c[0];
    v1.c[1] = P2.c[1] - P1.c[1];
    v1.c[2] = P2.c[2] - P1.c[2];

    v2.c[0] = P->c[0] - P1.c[0];
    v2.c[1] = P->c[1] - P1.c[1];
    v2.c[2] = P->c[2] - P1.c[2];



    // Compute dot products
    double dot00, dot01, dot02, dot11, dot12;
    dot00 = dotProduct3D(v0, v0);
    dot01 = dotProduct3D(v0, v1);
    dot02 = dotProduct3D(v0, v2);
    dot11 = dotProduct3D(v1, v1);
    dot12 = dotProduct3D(v1, v2);

    // Compute barycentric coordinates
    double invDenom, u, v;
    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v < 1);

}


