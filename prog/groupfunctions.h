#ifndef GROUPFUNCTIONS_H
#define GROUPFUNCTIONS_H

#define PI 3.14159265

// Tupac: 2D mesh rotation
void rotation2D(Mesh *mesh, float angle);

// Tupac: defines the center of a 2D mesh
// I don't know if it works
void center2D(Mesh *mesh, float *xc, float *yc);

/*** FUNCTION SUPERPOSITION ***/
/*	Parameters : 3 meshs: the 2 meshs to combine + the combinaison of the 2
	return : 1 if it works, 0 if not
	This function will join a mesh to another 
	Before this function The function loadMesh has to be done on the two meshs to combine, 
	in order to fill the number of Vertices, Triangles and edges of each one.
	After that we will go through the file of Mesh1 as a text file and find where is the last point
*/

int Superposition(pMesh Mesh1, pMesh Mesh2, pMesh Mesh_final ) ;

void translation2D(Mesh *mesh, float lengthX, float lengthY);

#endif
