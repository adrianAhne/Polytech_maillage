#ifndef GROUPFUNCTIONS_H
#define GROUPFUNCTIONS_H

#define PI         3.14159265
#define FLT_MAX    1e7
#define D_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define D_MIN(a,b)   ( ((a) > (b)) ? (b) : (a) )


// Tupac: 2D mesh rotation
void rotation2D(pMesh mesh, double angle);

void rotation3D(pMesh mesh, double angleX, double angleY, double angleZ);

// Tupac: defines the center of a 2D mesh
// I don't know if it works
int center(pMesh mesh, double *c);

/*** FUNCTION SUPERPOSITION ***/
/*	Parameters : 3 meshs: the 2 meshs to combine + the combinaison of the 2
	return : 1 if it works, 0 if not
	This function will join a mesh to another 
	Before this function The function loadMesh has to be done on the two meshs to combine, 
	in order to fill the number of Vertices, Triangles and edges of each one.
	After that we will go through the file of Mesh1 as a text file and find where is the last point
*/

int Superposition(pMesh Mesh1, pMesh Mesh2, pMesh Mesh_final) ;

void translation2D(pMesh mesh, double lengthX, double lengthY);

// calculate a new mesh translated by length length
void translation3D(Mesh *mesh, double lengthX, double lengthY, double lengthZ);

#endif
