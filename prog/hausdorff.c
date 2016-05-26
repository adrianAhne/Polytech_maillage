#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // you need to add -lm for link edition
//#include <omp.h> // parallel computing
#include "mesh.h"
#include "bucket.h"
#include "hash.h"
#include "distance.h"
#include "ball.h"


// calculate hausdorff distance between two meshes
double Hausdorff(pMesh meshA, pMesh meshB )
{
  printf("dd\n");

  
	int i,j,N;
	int* VertToTriaA = (int*)calloc(meshA->np+1, sizeof(int));
	int* VertToTriaB = (int*)calloc(meshB->np+1, sizeof(int));
	hashTria(meshA, VertToTriaA);
	hashTria(meshB, VertToTriaB);
	 		printf("dd\n"); 
	Bucket bucketA, bucketB;
	Point point;
	point.c[0] = 0 ;
	point.c[1] = 0 ;
	point.c[2] = 0 ;

	positive_boundingbox( meshA ,  &point);
	positive_boundingbox(meshB ,  &point);

	fprintf(stdout,"\n  -- Creation bucket MESH \n\nPlease type the number of subdivision : \n");
	fflush(stdin);	// just to be able to do a scanf without conflicts
	fscanf(stdin,"%d",&N);
	bucketB.size = N ;
	init_bucket( &bucketB , meshB ); 
	fill_bucket( &bucketB , meshB );
	



	fprintf(stdout,"\n  -- Creation bucket MESH \n\nPlease type the number of subdivision : \n");
	fflush(stdin);	// just to be able to do a scanf without conflicts
	fscanf(stdin,"%d",&N);
	bucketA.size = N ;
	init_bucket( &bucketA , meshA ); 
	fill_bucket( &bucketA , meshA );

	Hedge *tabA = (Hedge*)calloc(3*meshA->nt+1,sizeof(Hedge));	
	hashHedge(meshA, tabA);
	setAdj(meshA, tabA);
  Hedge *tabB = (Hedge*)calloc(3*meshB->nt+1,sizeof(Hedge));	
	hashHedge(meshB, tabB);
	setAdj(meshB, tabB);

	// calculate shortest distance from the mesh A to mesh B
	double distAB = distanceUsingBucket(meshB, &meshA->point[1], VertToTriaB, &bucketB );
	double distCurrent;
	for(i=2; i<=meshA->np; i++)
	{
		printf(" INDICE DE LA BOUCLE = %d\n", i);
		distCurrent = distanceUsingBucket(meshB, &meshA->point[i], VertToTriaB , &bucketB);
		printf("DISTANCE ENTRE %d et MESH B = %f \n", i, distCurrent);
		if (distCurrent > distAB)
	  		distAB = distCurrent;
	}
	printf("d444\n");

	// calculate shortest distance from meshB to meshA
	double distBA = distanceUsingBucket(meshA, &meshB->point[1], VertToTriaA, &bucketA);
	for(i=2; i<=meshB->np; i++)
	{
		distCurrent = distanceUsingBucket(meshA, &meshB->point[i], VertToTriaA , &bucketA);
		if (distCurrent > distBA)
	  		distBA = distCurrent;
	}
	printf("g\n");
	// hausdorff distance is the maximum of both calculated distances
	printf("distAB = %f et dist BA = %f \n", distAB, distBA);
	double distHausdorff = (distAB > distBA) ? distAB : distBA;

	free(VertToTriaA);
	free(VertToTriaB);

	return distHausdorff;
}
