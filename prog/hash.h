#ifndef HASH_H
#define HASH_H

/*
Hedge struct doc:
na and nb are the endpoints of the edge

ia = min(na,nb)
ib = max(na,nb)

This objects Hedge are ordered in a global array

key = (KA*min(na,nb)+KB*max(na,nb))%hsize
KA = 7, KB = 11, hsize = mesh->np
*/

typedef struct {
	int ia, ib, adj1, adj2, nxt;
} Hedge;
typedef Hedge * pHedge;

/* Construct the table hachage */
int hashHedge(pMesh mesh, Hedge * tab);

/* Calculate neighbour/adjacent triangles of a given triangle (only at edges -> max 3 triangles)*/
int setAdj(pMesh mesh, Hedge * tab);

/* Brute force method to find the index of a triangle for a given point */
int localiseTriangleBruteForce(pMesh mesh, pPoint point);

/* calculates coordinates barycentric u, v, 1-u-v for a given point p and saves them in the pointer cb */
void baryCoord(pMesh mesh, int triangle, pPoint p, double cb[3]);

/* find triangle who includes the point p by starting at a starttriangle and approaches the searched triangles via adjacents */
int locelt(pMesh mesh, int startTriangle, pPoint p, double cb[3]);

/*	 Following the algorithm given in github/doc/main.pdf
	 Compute the distance between a point and a surface using buckets and adjacences
	 VertToTria gives an array which tells us to which triangle a given point belongs to */
double distanceUsingBucket(pMesh mesh, pPoint p, int *VertToTria , pBucket bucket);

#endif
