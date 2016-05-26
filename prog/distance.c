#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "rigidtransformation.h"
#include "distance.h"
#include "bucket.h"
#include "hash.h"


// calculate dot product of two three-dimensional vectors
double dotProduct3D(Point a, Point b)
{
    double sum = 0;
    int i;
    for (i = 0; i < 3; i++) 
		sum += a.c[i] * b.c[i];

    return sum;
}

/*
// Check if point is in triangle tria (NOT USED CURRENTLY)
void checkPointInTriangle(pMesh mesh, pPoint P, pTria tria, double* re, double* dotProduct)
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
    double dot[5];
    dot00 = dotProduct3D(v0, v0);
    dot[0] = dot00;
    dot01 = dotProduct3D(v0, v1);

    dot02 = dotProduct3D(v0, v2);
    dot11 = dotProduct3D(v1, v1);
    dot12 = dotProduct3D(v1, v2);

    // Compute barycentric coordinates
    double invDenom, u, v;
    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    printf("u = %f, v = %f, 1-u-v=%f", u, v, 1-u-v);
    // Check if point is in triangle
    
    re[0] = u;
    re[1] = v;
    re[2] = 1-u-v;


}*/


/*
Point calculeProjection(pMesh mesh, pTria tria, pPoint PO)
{

    // check if normales of triangles are already calculated
    if (!mesh->triaNorm)
    {
        normalesOfTriangles(mesh);
    }

    Point P1P0, Np, POproj, proj;
    double lenPointToProj;

    // Vector P1P0
    P1P0.c[0] = mesh->point[tria->v[0]].c[0] - PO->c[0];
    P1P0.c[1] = mesh->point[tria->v[0]].c[1] - PO->c[1];
    P1P0.c[2] = mesh->point[tria->v[0]].c[2] - PO->c[2];

    // The normal of triangle
    Np.c[0] = mesh->triaNorm[tria->ref].n[0];
    Np.c[1] = mesh->triaNorm[tria->ref].n[1];
    Np.c[2] = mesh->triaNorm[tria->ref].n[2];

    // cos(alpha) between the normal and the vector P1P0
    double cosalpha;
    cosalpha = dotProduct3D(P1P0, Np)/((sqrt(pow(P1P0.c[0], 2) + pow(P1P0.c[1], 2) + pow(P1P0.c[2], 2))) *
    sqrt(pow(Np.c[0], 2) + pow(Np.c[1], 2) + pow(Np.c[2], 2)));

    // lenght of the vector P0 to its projection
    lenPointToProj = sqrt(pow(P1P0.c[0], 2) + pow(P1P0.c[1], 2) + pow(P1P0.c[2], 2)) * cosalpha;

    // Vector P0 to projection
    POproj.c[0] = - lenPointToProj * ( Np.c[0] / sqrt(pow(Np.c[0], 2) + pow(Np.c[1], 2) + pow(Np.c[2], 2) ));
    POproj.c[1] = - lenPointToProj * ( Np.c[1] / sqrt(pow(Np.c[0], 2) + pow(Np.c[1], 2) + pow(Np.c[2], 2) ));
    POproj.c[2] = - lenPointToProj * ( Np.c[2] / sqrt(pow(Np.c[0], 2) + pow(Np.c[1], 2) + pow(Np.c[2], 2) ));

    // Calcul projection
    proj.c[0] = PO->c[0] + POproj.c[0];
    proj.c[1] = PO->c[1] + POproj.c[1];
    proj.c[2] = PO->c[2] + POproj.c[2];

    return proj;

}
*/



/* calculates the distance between a point and a triangle by using barycentric coordinates */
double distPointToTriangle(pMesh mesh, pTria tria, pPoint P0)
{

    Point P1, P2, P3;
    P1 = mesh->point[(tria->v[0])];
    P2 = mesh->point[(tria->v[1])];
    P3 = mesh->point[(tria->v[2])];

	// Compute vectors
	//    v0 = P3 - P1
	//    v1 = P2 - P1
	//    v2 = P - P1
    Point v0, v1, v2, D;
    
    v0.c[0] = P1.c[0];
    v0.c[1] = P1.c[1];
    v0.c[2] = P1.c[2];

    v1.c[0] = P2.c[0] - P1.c[0];
    v1.c[1] = P2.c[1] - P1.c[1];
    v1.c[2] = P2.c[2] - P1.c[2];

    v2.c[0] = P3.c[0] - P1.c[0];
    v2.c[1] = P3.c[1] - P1.c[1];
    v2.c[2] = P3.c[2] - P1.c[2];
    
    D.c[0] = P1.c[0] - P0->c[0];
    D.c[1] = P1.c[1] - P0->c[1];
    D.c[2] = P1.c[2] - P0->c[2];


    // Calculate relations with the help of barycentric coordinates
    double a, b, c, d, e, f, sqrDistance, invDet, det, s, t, tmp0, tmp1, numer, denom, dist;
    a = dotProduct3D(v0, v0);
    b = dotProduct3D(v0, v1);
    c = dotProduct3D(v1, v1);
    d = dotProduct3D(v0, D);
    e = dotProduct3D(v1, D);
    f = dotProduct3D(D, D);
    
    det = a*c - b*b; 
    s = b*e - c*d;
    t = b*d - a*e;
    //printf("%f %f %f %f %f %f %f %f %f\n", a,b,c,d,e,f,det,s,t);


	// For the following tree this small picture might be useful to visualise what the tree is doing. We calculate some relations/barycentric coordinates concerning the triangle P1, P2, P3 and now depending on the signe of this calculated values we can determine the region/direction in which the point is lying concering triangle and so calculate the distance

    //  \     |
    //   \reg2|
    //    \   |
    //     \  |
    //      \ |
    //       \|
    //        *P3
    //        |\
    //        | \
    //  reg3  |  \ reg1
    //        |   \
    //        |reg0\ 
    //        |     \ 
    //        |      \ P2
    // -------*-------*------->s
    //        |P1      \ 
    //  reg4  | reg5    \ reg6

    if ((s+t) <= det) {
        if (s < 0) {
            if (t < 0) {
                // region4
                //printf("Region4\n");
                if (d < 0) {
                    t = 0;
                    if (-d >= a) {
                        s = 1;
                        sqrDistance = a + 2*d + f;
                    } else {
                        s = -d/a;
                        sqrDistance = d*s + f;
                    }
                } else {
                    s = 0;
                    if (e >= 0) {
                        t = 0;
                        sqrDistance = f;
                    } else {
                        if (-e >= c) {
                            t = 1;
                            sqrDistance = c + 2*e + f;
                        } else {
                            t = -e/c;
                            sqrDistance = e*t + f;
                        }
                    }
                }
                // end of region 4
            } else {
                // region 3
                //printf("Region3\n");
                s = 0;
                if (e >= 0) {
                    t = 0;
                    sqrDistance = f;
                } else {
                    if (-e >= c) {
                        t = 1;
                        sqrDistance = c + 2*e +f;
                    } else {
                        t = -e/c;
                        sqrDistance = e*t + f;
                    }
                }
            }
        }
        // end of region 3
        else {
            if (t < 0) {
                // region 5
                //printf("Region5\n");
                t = 0;
                if (d >= 0) {
                    s = 0;
                    sqrDistance = f;
                } else {
                    if (-d >= a) {
                        s = 1;
                        sqrDistance = a + 2*d + f;
                    } else {
                        s = -d/a;
                        sqrDistance = d*s + f;
                    }
                }
            } else {
                // region 0
                invDet = 1/det;
                s = s*invDet;
                t = t*invDet;
                sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
            }
        }
    } else {
        if (s < 0) {
            // region 2
            //printf("Region2\n");
            tmp0 = b + d;
            tmp1 = c + e;
            if (tmp1 > tmp0) { // minimum on edge s+t=1
                numer = tmp1 - tmp0;
                denom = a - 2*b + c;
                if (numer >= denom) {
                    s = 1;
                    t = 0;
                    sqrDistance = a + 2*d + f;
                } else {
                    s = numer/denom;
                    t = 1-s;
                    sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                }
            } else { // minimum on edge s=0
                s = 0;
                if (tmp1 <= 0) {
                    t = 1;
                    sqrDistance = c + 2*e + f;
                } else {
                    if (e >= 0) {
                        t = 0;
                        sqrDistance = f;
                    } else {
                        t = -e/c;
                        sqrDistance = e*t + f;
                    }
                }
            }
        }
        // end of region 2
        else {
            if (t < 0) {
                // region 6
                //printf("Region6\n");
                tmp0 = b + e;
                tmp1 = a + d;
                if (tmp1 > tmp0) {
                    numer = tmp1 - tmp0;
                    denom = a-2*b+c;
                    if (numer >= denom) {
                        t = 1;
                        s = 0;
                        sqrDistance = c + 2*e + f;
                    } else {
                        t = numer/denom;
                        s = 1 - t;
                        sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    }
                } else {
                    t = 0;
                    if (tmp1 <= 0) {
                        s = 1;
                        sqrDistance = a + 2*d + f;
                    } else {
                        if (d >= 0) {
                            s = 0;
                            sqrDistance = f;
                        } else {
                            s = -d/a;
                            sqrDistance = d*s + f;
                        }
                    }
                }
            }
            // end of region 6
            else {
                // region 1
                //printf("Region1\n");
                numer = c + e - b - d;
                if (numer <= 0) {
                    s = 0;
                    t = 1;
                    sqrDistance = c + 2*e + f;
                } else {
                    denom = a - 2*b + c;
                    if (numer >= denom) {
                        s = 1;
                        t = 0;
                        sqrDistance = a + 2*d + f;
                    } else {
                        s = numer/denom;
                        t = 1-s;
                        sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    }
                }
            }
        }
    }

    if (sqrDistance < 0)
        sqrDistance = 0;
    
    dist = sqrt(sqrDistance);
    //printf("Distance from point to triangle : %f\n", dist);

    return dist;
}


/* Calculating the average distance of a point to the triangle */
double averageDistancePTT(pMesh mesh, pTria tria, pPoint P0)
{
    int i;
    Tria temp;
    double dAverage = 0;
    for (i = 0; i < 3; i++)
    {
        temp.v[0] = tria->v[i%3];
        temp.v[1] = tria->v[(i+1)%3];
        temp.v[2] = tria->v[(i+2)%3];
        
        dAverage += distPointToTriangle(mesh, &temp, P0);
    }

    dAverage /= 3;
    return dAverage;
}

