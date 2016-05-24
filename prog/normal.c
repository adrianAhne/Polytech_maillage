#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "normal.h"


/* Calculates the normals of the triangles with the help of the cross product*/
void normalesOfTriangles(pMesh mesh)
{
	// allocate memory
	mesh->triaNorm = (pTriaNorm)calloc(mesh->nt+1, sizeof(TriaNorm));
	assert(mesh->triaNorm);
	
	mesh->Normal = (pNormal)calloc(mesh->np+1, sizeof(Normal));
	assert(mesh->Normal);

	int i, j, k;
	pTria currentTria;
	Point P1,P2,P3,N;
	double weight, norm;

	// Loop over all triangles calculating the normales, the area/surface and store data into mesh->triaNorm
	for(i=1; i <= mesh->nt; i++)
	{
		currentTria = &mesh->tria[i];
		
		P1 = mesh->point[(currentTria->v[0])];
		P2 = mesh->point[(currentTria->v[1])];
		P3 = mesh->point[(currentTria->v[2])];
		
		// Normal of the triangle
		N.c[0] = (P2.c[1]-P1.c[1])*(P3.c[2]-P1.c[2]) - (P2.c[2]-P1.c[2])*(P3.c[1]-P1.c[1]);
		N.c[1] = (P2.c[2]-P1.c[2])*(P3.c[0]-P1.c[0]) - (P2.c[0]-P1.c[0])*(P3.c[2]-P1.c[2]);
		N.c[2] = (P2.c[0]-P1.c[0])*(P3.c[1]-P1.c[1]) - (P2.c[1]-P1.c[1])*(P3.c[0]-P1.c[0]);

		(mesh->triaNorm[i]).n[0] = N.c[0];
		(mesh->triaNorm[i]).n[1] = N.c[1];
		(mesh->triaNorm[i]).n[2] = N.c[2];
		//printf("%f %f %f\n", N.c[0],N.c[1],N.c[2]);

		// Weight (= area) of the triangle 
		weight = 0.5 * sqrt( pow(N.c[0],2)+pow(N.c[1],2)+pow(N.c[2],2) );
		(mesh->triaNorm[i]).weight = weight;
	}


	// Loop over all points and calculate weighted normales for each point by using the touching/neighbour triangles of the point
	for (i = 1; i < mesh->nt ; i++) {
		for (j = 0; j < 3 ; j++) {
			mesh->Normal[mesh->tria[i].v[j]].n[0] += mesh->triaNorm[i].n[0];
			mesh->Normal[mesh->tria[i].v[j]].n[1] += mesh->triaNorm[i].n[1];
			mesh->Normal[mesh->tria[i].v[j]].n[2] += mesh->triaNorm[i].n[2];
		}
	}
	
	for (i = 1; i < mesh->np; ++i)
	{
		norm = sqrt( pow(mesh->Normal[i].n[0],2) + pow(mesh->Normal[i].n[1],2) + pow(mesh->Normal[i].n[2],2) );
		for (j=0;j<3;j++) {
			mesh->Normal[i].n[j] *= 1 / norm;
		}
	}

	mesh->nn = mesh->np;
}



