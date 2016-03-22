/*	
		CONSTRUCTION BUCKET
		
		Jean-Tupac Quiroga
		Adrian Ahne
		Alexandre POULAIN
		
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // you need to add -lm for link edition
#include "mesh.h"
#include "bucket.h"


void init_bucket ( pBucket bucket , pMesh mesh)
{
	
	bucket->head = calloc(pow(bucket->size,3),sizeof(int));
	bucket->link = calloc(mesh->np , sizeof(int));
	
}

void fill_bucket( pBucket bucket , pMesh mesh , int N )
{
	int indice,i,j,k, key  ;
	

	/* go through all the points of the mesh */
	for ( indice = 1 ; indice <= mesh->np ; indice ++ )
	{ 	
	
		
		/* calcul i,j,k */
		i = max(0,(int)(N*mesh->point[indice].c[0])-1) ; 
		j = max(0,(int)(N*mesh->point[indice].c[1])-1) ;
		k = max(0,(int)(N*mesh->point[indice].c[2])-1) ;
		
		/* calcul of the key */ 
		key = (j*N+k)*N+i ;
		fprintf(stdout," key = %d \n " , key ) ;
		
		/* we test the key in the bucket */
		
		if ( bucket->head[key] == 0 )
		{
			fprintf(stdout,"  ICI \n");
			bucket->head[key] = indice ;
		}
		else
		{
			bucket->link[indice] = bucket->head[key] ;
			bucket->head[key] = indice ;
		}
		
	}
	/* end for */
	
}

void free_bucket (pBucket bucket)
{
	free( bucket->head );
	free( bucket->link );
}







