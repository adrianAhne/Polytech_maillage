#ifndef CURVATURE_H
#define CURVATURE_H

#define PI         3.14159265

/*
	Calculates the curvature for each point by adding up
	all the angles of the triangles around the point
	and substract this sum from 2*PI
	It creates the sol file in the current working directory	
	3 dimensional
*/
int curvature3D(pMesh mesh);

/*
	Calculates the curvature for each point by adding up
	all the angles of the triangles around the point
	and substract this sum from 2*PI
	It creates the sol file in the current working directory	
	2 dimensional
*/
int curvature2D(pMesh mesh);

#endif
