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


/* FUNCTION distbuck 
		This function will calculate the minimum distance between a point of the mesh A and the mesh B 
		PARAMETERS: a pointer to the point p of the mesh A, a pointer to the bucket of the mesh B, a pointer to the mesh B 
		RETURN : the distance d 
*/ 
double distbuck( pPoint p , pBucket bucket_meshB , pMesh meshB );

#endif
