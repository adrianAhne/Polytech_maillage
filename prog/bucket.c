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
	
	bucket->head = calloc(pow(bucket->size,3)+1,sizeof(int));
	bucket->link = calloc(mesh->np+1 , sizeof(int));
	
}

void fill_bucket( pBucket bucket , pMesh mesh)
{
	int indice,i,j,k, key;
	

	/* go through all the points of the mesh */
	for ( indice = 1 ; indice <= mesh->np ; indice ++ )
	{
	
    /* calcul i,j,k */
		i = max( 0, (int)(bucket->size*mesh->point[5207].c[0]) -1) ;
		j = max( 0, (int)(bucket->size*mesh->point[5207].c[1]) -1) ;
		k = max( 0, (int)(bucket->size*mesh->point[5207].c[2]) -1) ;
		
		
		/* calcul of the key */ 
		key = (j*bucket->size+k)*bucket->size+i ;
			
		
		/* we test the key in the bucket */
		
		if ( !bucket->head[key]  )
		{
			
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

/* FUNCTION point_min 
In order to determine the minimum point  
PARAMETERS : a pMesh , a char for the axis
*/
double point_min (pMesh mesh , char axis)
{
	int i ;
	double min ;
	 
	
	/* for the x axis */
	
	if ( axis == 'x' )
	{
		min = mesh->point[1].c[0];
		for ( i = 1 ; i <= mesh->np ; i++ )
		{
			if ( mesh->point[i].c[0] < min )
				min = mesh->point[i].c[0] ;
		}
	}
	if ( axis == 'y' )
	{
		min = mesh->point[1].c[1];
		for ( i = 1 ; i <= mesh->np ; i++ )
		{
			if ( mesh->point[i].c[1] < min )
				min = mesh->point[i].c[1] ;
		}
	}
	if ( axis == 'z' )
	{
		min = mesh->point[1].c[2];
		for ( i = 1 ; i <= mesh->np ; i++ )
		{
			if ( mesh->point[i].c[2] < min )
				min = mesh->point[i].c[2] ;
		}
	}
	return min ;
}

/* FUNCTION positive_boundingbox 
This function will translate the mesh using the 3D translation function in order to have a bounding box in the positive part of the axis (x,y,z)
PARAMETERS : a pointer to a mesh 
*/
void positive_boundingbox( pMesh mesh , pPoint point)
{
	double xmin , ymin, zmin ;
	
	/* First of all we have to determine the length of the translation */
	/* so we determine the bounding box */
	xmin = point_min ( mesh , 'x');
	ymin = point_min ( mesh , 'y');
	zmin = point_min ( mesh , 'z');
	//fprintf(stdout , " xmin = %f \n ymin = %f \n zmin = %f \n", xmin , ymin , zmin ) ;
	
	if ( xmin >= 0 )
		xmin = 0 ;
	if ( ymin >= 0 )
		ymin = 0 ;
	if ( zmin >= 0 )
		zmin = 0 ;
		
	/* now we translate */
	translation3D(mesh, abs(xmin), abs(ymin), abs(zmin));	
	
	/* and the point */
	point->c[0] += xmin ;
	point->c[1] += ymin ;
	point->c[2] += zmin ;
	
	/* from now on we have a positive bounding box: GOOD */
	
	
}

/* FUNCTION use_bucket 
		This function will use the coordinates of a point and associate the key.
		With this key, we will define the neighbourhood  of the point.
		Parameters : the bucket the mesh  and a point
		Return : the point in the bucket ( indice )
*/

int use_bucket( pBucket bucket , pMesh mesh ,  pPoint point , double increment )
{
		
		
		int i,j,k,key,N,cherche,indice = 0 , ret ;
    pPoint  pp1,pp2;
  
			/* calcul i,j,k */
		i = max(0,(int)(bucket->size*point->c[0])-1) + increment ; 
		j = max(0,(int)(bucket->size*point->c[1])-1) + increment;
		k = max(0,(int)(bucket->size*point->c[2])-1) + increment;
		
		key = (j*bucket->size+k)*bucket->size+i ;
		//fprintf(stdout,"  key = %d \n",key);
		
		
		
		/* we check if the tab head has a point in it */
		if ( bucket->head[key] )
		{
      fprintf(stdout,"  Point in Bucket head \n ");
      ret     = bucket->head[key];
      cherche = ret ;
      pp1     = &mesh->point[cherche];
      
      
			/* on affiche les points du même subdomain */
			while( bucket->link[ cherche ])
			{
				cherche = bucket->link[cherche] ;
        pp2 = &mesh->point[cherche];
        fprintf(stdout,"  Point in Bucket link  x = %f  y = %f  z =  %f \n", pp2->c[0],pp2->c[1],pp2->c[2]);
        indice ++ ;
				
			}
			
		} 
		/* else there is no neighbourgh in this subdomain so we explore the subdomain around */
		else
    {
      fprintf(stdout,"  The point is not in the current cell \n");
      increment ++ ;
      ret = use_bucket( bucket , mesh , point , increment ) ; 
    }
     pp1     = &mesh->point[ret];
    fprintf(stdout," indice = %d \n" , indice );
    fprintf(stdout,"  Point in Bucket head  x = %f  y = %f  z =  %f \n", pp1->c[0],pp1->c[1],pp1->c[2]);
    return ret ;
    
		
}



void free_bucket (pBucket bucket)
{
	free( bucket->head );
	free( bucket->link );
}







