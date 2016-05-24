#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // you need to add -lm for link edition
//#include <omp.h> // parallel computing
#include "mesh.h"
#include "hash.h"
#include "distanceMeshFunctions.h"
#include "bucket.h"


// calculate hausdorff distance between two meshes
double Hausdorff(pMesh meshA, pMesh meshB)
{

	int i,j;
	int* VertToTriaA = (int*)calloc(meshA->np+1, sizeof(int));
	int* VertToTriaB = (int*)calloc(meshB->np+1, sizeof(int));
	hashTria(&meshA, VertToTriaA);
	hashTria(&meshB, VertToTriaB);
	printf("dd\n");
	
	// calculate shortest distance from the mesh A to mesh B
	double distAB = distanceUsingBucket(meshB, &meshA->point[1], VertToTriaB);
	double distCurrent;
	for(i=2; i<meshA->np; i++)
	{
		printf("%d\n", i);
		distCurrent = distanceUsingBucket(meshB, &meshA->point[i], VertToTriaB);
		if (distCurrent > distAB)
	  		distAB = distCurrent;
	}
	printf("d444\n");

	// calculate shortest distance from meshB to meshA
	double distBA = distanceUsingBucket(meshA, &meshB->point[1], VertToTriaA);
	for(i=2; i<meshB->np; i++)
	{
		distCurrent = distanceUsingBucket(meshA, &meshB->point[i], VertToTriaA);
		if (distCurrent > distBA)
	  		distBA = distCurrent;
	}
	printf("g\n");
	// hausdorff distance is the maximum of both calculated distances
	double distHausdorff = distAB > distBA ? distAB : distBA;

	free(VertToTriaA);
	free(VertToTriaB);

	return distHausdorff;
}
