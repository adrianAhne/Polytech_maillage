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

  
	int i,j,N;
	/*int* VertToTriaA = (int*)calloc(meshA->np+1, sizeof(int));
	int* VertToTriaB = (int*)calloc(meshB->np+1, sizeof(int));
	hashTria(meshA, VertToTriaA);
	hashTria(meshB, VertToTriaB);
*/
	// create bucket
	Bucket bucketA, bucketB;
	Point point;
	point.c[0] = 0 ;
	point.c[1] = 0 ;
	point.c[2] = 0 ;

	// bound the points between the interval [0,1]
	positive_boundingbox( meshA ,  &point);
	positive_boundingbox(meshB ,  &point);

	fprintf(stdout,"\n  -- Creation bucket MESH_A \n\nPlease type the number of subdivision : \n");
	fflush(stdin);	// just to be able to do a scanf without conflicts
	fscanf(stdin,"%d",&N);
	bucketA.size = N ;
	init_bucket( &bucketA , meshA ); 
	fill_bucket( &bucketA , meshA );



	fprintf(stdout,"\n  -- Creation bucket MESH_B \n\nPlease type the number of subdivision : \n");
	fflush(stdin);	// just to be able to do a scanf without conflicts
	fscanf(stdin,"%d",&N);
	bucketB.size = N ;
	init_bucket( &bucketB , meshB ); 
	fill_bucket( &bucketB , meshB );
	



	// create table of hachage
	Hedge *tabA = (Hedge*)calloc(3*meshA->nt+1,sizeof(Hedge));	
	hashHedge(meshA, tabA);
	setAdj(meshA, tabA);
 	Hedge *tabB = (Hedge*)calloc(3*meshB->nt+1,sizeof(Hedge));	

	hashHedge(meshB, tabB);
	setAdj(meshB, tabB);

	// calculate shortest distance from the mesh A to mesh B
	double distAB = distbuck( &meshA->point[1] , &bucketB , meshB );
	double distCurrent;
	for(i=2; i<=meshA->np; i++)
	{
		distCurrent = distbuck( &meshA->point[i] , &bucketB , meshB );
		//fprintf(stdout," INDICE DE LA BOUCLE = %d \nDISTANCE ENTRE LE PT %d DE MAILLAGE A ET LE MAILLAGE B = %f\n",i,i,distCurrent);
		if (distCurrent > distAB)
	  		distAB = distCurrent;
	}
	//printf( "distAB = %f \n " , distAB ) ;

	// calculate shortest distance from meshB to meshA
	double distBA = distbuck( &meshB->point[1] , &bucketA , meshA );
	for(i=2; i<=meshB->np; i++)
	{
		distCurrent = distbuck( &meshB->point[i] , &bucketA , meshA );
		if (distCurrent > distBA)
	  		distBA = distCurrent;
	}
	//printf( "distBA = %f \n" , distBA ) ;
	// hausdorff distance is the maximum of both calculated distances
	double distHausdorff ;
	if ( distAB > distBA ) 
		return distAB ;
	else 
		return distBA ;

	/*free(VertToTriaA);
	free(VertToTriaB);
*/
	return distHausdorff;
}
