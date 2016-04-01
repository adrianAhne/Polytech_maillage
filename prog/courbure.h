#ifndef COURBURE_H
#define COURBURE_H

#define PI         3.14159265

/*** FUNCTION COURBURE 3D***/
/* 
	This functions will calculate the curve for each point 
	For that we will add up all the angle of the triangles around the point and make the difference to 2PI 
	Parameters : a mesh 
	return : it create the sol file in the current repository
						it return if it works ,  if not 
*/
int courbure3D(pMesh mesh ) ;

/*** FUNCTION COURBURE 2D***/ 
/* In this function we will calculate the curvature for each point in a 2D Mesh.
		In order to do that we will use the gaussian curvature: this is the angular between two edge.
		PARAMETERS : the 2D mesh ( pMesh )
		RETURN : 1 if everythings happen alright 
*/
int courbure2D( pMesh mesh );


#endif
