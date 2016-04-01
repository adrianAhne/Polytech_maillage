#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "normale.h"

void normalesOfTriangles(pMesh mesh)
{

	mesh->triaNorm = (pTriaNorm)calloc(mesh->nt+1, sizeof(TriaNorm));
	assert(mesh->triaNorm);
	
	mesh->Normal = (pNormal)calloc(mesh->np+1, sizeof(Normal));
	assert(mesh->Normal);

	int i, j, k;
	pTria currentTria;
	Point P1,P2,P3,N;
	double weight, norm;

	// Loop over all triangles calculating the normales, the area and store data into mesh->triaNorm
	for(i=1; i <= mesh->nt; i++)
	{
		currentTria = &mesh->tria[i];
		
		/*
		In this part, I'm calculating the normal of the face of the triangle
		After, we'll average each normal
		For this technique, we'll have to calculate the weighting of each face (= surface ratio)

		Abstract : normal of 1 triangle (composed of 3 points P1, P2 and P3) :

		The cross product of two sides of the triangle equals the surface normal. 
		So, if V = P2 - P1 and W = P3 - P1, and N is the surface normal, then:

		Nx=(Vy∗Wz)−(Vz∗Wy)=((P2y-P1y)*(P3z-P1z)-(P2z-P1z)*(P3y-P1y))
		Ny=(Vz∗Wx)−(Vx∗Wz)=((P2z-P1z)*(P3x-P1x)-(P2x-P1x)*(P3z-P1z))
		Nz=(Vx∗Wy)−(Vy∗Wx)=((P2x-P1x)*(P3y-P1y)-(P2y-P1y)*(P3x-P1x))

		Weight=AreaOfTheTriangle=0.5*sqrt((Vy*Wz-Vz*Wy)^2+(Vz*Wx-Vx*Wz)^2+(Vx*Wy-Vx*Wx)^2)

		Translation in code :
		*/

		P1 = mesh->point[(currentTria->v[0])];
		P2 = mesh->point[(currentTria->v[1])];
		P3 = mesh->point[(currentTria->v[2])];
		
		// Normal of the triangle
		N.c[0] = (P2.c[1]-P1.c[1])*(P3.c[2]-P1.c[2]) - (P2.c[2]-P1.c[2])*(P3.c[1]-P1.c[1]);
		N.c[1] = (P2.c[2]-P1.c[2])*(P3.c[0]-P1.c[0]) - (P2.c[0]-P1.c[0])*(P3.c[2]-P1.c[2]);
		N.c[2] = (P2.c[0]-P1.c[0])*(P3.c[1]-P1.c[1]) - (P2.c[1]-P1.c[1])*(P3.c[0]-P1.c[0]);

		// Weight of the triangle (let's say it the area)
		weight = 0.5 * sqrt( pow(N.c[0],2)+pow(N.c[1],2)+pow(N.c[2],2) );
		
		// Normale is calculated at point P1 for the triangle; save both points defining the normale
		(mesh->triaNorm[i]).n[0] = N.c[0];
		(mesh->triaNorm[i]).n[1] = N.c[1];
		(mesh->triaNorm[i]).n[2] = N.c[2];
		printf("%f %f %f\n", N.c[0],N.c[1],N.c[2]);
		(mesh->triaNorm[i]).weight = weight;
		

	}


	// REMARQUES
	// Faire plutot une boucle sur les triangles
	// Calculer l'aire à chaque boucle après
	//
	// A la fin du calcul de normales, normaliser chaque normale pour chaque point par sa norme

	
	// Loop over all points and calculate weighted normales for each point by using neighbour triangles
	// WARNING : Not optimized and not normalized normales ! 
	
	
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



