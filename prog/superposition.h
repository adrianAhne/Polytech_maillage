
#ifndef SUPERPOSITION_H
#define SUPERPOSITION_H


/*** Change 2D to 3D ***/
/* Parameter : a mesh with dim == 2 
	This function will add a square to the tab of points for the third dimension
*/
void Change2Dto3D( pMesh Mesh );


/*** FUNCTION SUPERPOSITION ***/
/*	Parameters : 3 meshs: the 2 meshs to combine + the combinaison of the 2
	return : 1 if it works, 0 if not
	This function will join a mesh to another 
	Before this function The function loadMesh has to be done on the two meshs to combine, 
	in order to fill the number of Vertices, Triangles and edges of each one.
	After that we will go through the file of Mesh1 as a text file and find where is the last point
*/

int Superposition(pMesh Mesh1, pMesh Mesh2, pMesh Mesh_final ) ;







#endif


