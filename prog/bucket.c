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
#include "rigidtransformation.h"

/* initialise the bucket*/
void init_bucket ( pBucket bucket , pMesh mesh)
{
	bucket->head = calloc(pow(bucket->size,3)+1,sizeof(int));
	bucket->link = calloc(mesh->np+1 , sizeof(int));
}


/* Fill the 3 objects of a bucket
	First we will go through all the points of the mesh and calculate for each point:
						i = max(0,(int)(N*p->c[0])-1) 
						j = max(0,(int)(N*p->c[1])-1) 
						k = max(0,(int)(N*p->c[2])-1)
	Each of the integer is in the interval [0, N-1]
	Via this number we can associate each point with a key = (k*N+j)*N+i				
	With this key we go to bucket->head[key] and check : 
		- if  bucket->head[key] == 0, we put p in it 
		- if  bucket->head[key] > 0 , there is already a point (= Old_point) 
			So we put the old point at bucket->link[p] = bucket->head[key] 
			and we put the new point in head bucket->head[key] = p	   */
void fill_bucket( pBucket bucket , pMesh mesh)
{
	int indice,i,j,k, key;
	

	/* go through all the points of the mesh */
	for ( indice = 1 ; indice <= mesh->np ; indice ++ )
	{
	
    	/* calcul i,j,k */
		i = max( 0, (int)(bucket->size*mesh->point[indice].c[0]) -1) ;
		j = max( 0, (int)(bucket->size*mesh->point[indice].c[1]) -1) ;
		k = max( 0, (int)(bucket->size*mesh->point[indice].c[2]) -1) ;

		/* calcul of the key */ 
		key = (j*bucket->size+k)*bucket->size+i ;
		
		/* we test the key in the bucket */
		if ( !bucket->head[key]  )
			bucket->head[key] = indice ;
		else
		{
			bucket->link[indice] = bucket->head[key] ;
			bucket->head[key] = indice ;
		}
	 }
}


/* 
	Determine the minimum point
*/
double point_min (pMesh mesh , char axis)
{
	int i ;
	double min ;
	 
	
	// for the x axis
	if ( axis == 'x' )
	{
		min = mesh->point[1].c[0];
		for ( i = 1 ; i <= mesh->np ; i++ )
		{
			if ( mesh->point[i].c[0] < min )
				min = mesh->point[i].c[0] ;
		}
	}
	
	// for the y axis
	if ( axis == 'y' )
	{
		min = mesh->point[1].c[1];
		for ( i = 1 ; i <= mesh->np ; i++ )
		{
			if ( mesh->point[i].c[1] < min )
				min = mesh->point[i].c[1] ;
		}
	}
	
	// for the z axis
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
*/
void positive_boundingbox( pMesh mesh , pPoint point)
{
	double xmin , ymin, zmin ;
	int i ;
	/* First of all we have to determine the length of the translation */
	/* so we determine the bounding box */
	xmin = point_min ( mesh , 'x');
	ymin = point_min ( mesh , 'y');
	zmin = point_min ( mesh , 'z');
	
	if ( xmin < 0 || ymin < 0 || zmin < 0 ) 
	{
		
		/* now we translate */
		translation3D(mesh, fabs(xmin), fabs(ymin), fabs(zmin));
		
	
		/* and the point */
		point->c[0] = (point->c[0] + fabs(xmin))/10 ;
		point->c[1] = (point->c[1] + fabs(ymin))/10 ;
		point->c[2] = (point->c[2] + fabs(zmin))/10 ;
		
		
		for ( i = 1; i <=  mesh->np ; i++ )
		{
	
			mesh->point[i].c[0] = mesh->point[i].c[0] /10 ;
			mesh->point[i].c[1] = mesh->point[i].c[1] /10 ;
			mesh->point[i].c[2] = mesh->point[i].c[2] /10 ;
		} 
	}
}


/* 
	This function will use the coordinates of a point and associate the key.
	With this key, we will define the neighbourhood  of the point.
*/
int use_bucket( pBucket bucket , pMesh mesh ,  pPoint point , double increment )
{
		
		
	int i,j,k,key,N,cherche,indice = 0 , ret ;
	pPoint  pp1,pp2;
  

	/* calcul i,j,k */
	i = max(0,(int)((bucket->size)*(point->c[0]))-1 + increment) ; 
	j = max(0,(int)((bucket->size)*(point->c[1]))-1 + increment) ;
	k = max(0,(int)((bucket->size)*(point->c[2]))-1 + increment) ;
			
	key = (j*bucket->size+k)*bucket->size+i ;
		
	/* we check if the tab head has a point in it */
	if ( bucket->head[key] )
	{
      	ret     = bucket->head[key];
      	cherche = ret ;
      	pp1     = &mesh->point[cherche];
      
     	/* we take the points of the same subdomaine */
		while( bucket->link[ cherche ])
		{
			cherche = bucket->link[cherche] ;
    	    pp2 = &mesh->point[cherche];
    	}			
	} 
	/* else there is no neighbourgh in this subdomain so we explore the subdomain around */
	else
	{
		increment ++ ;
		ret = use_bucket( bucket , mesh , point , increment ) ; 
	}
	
	pp1     = &mesh->point[ret];

	return ret ;	
}


/* Function uses bucket to return the value of the associated key to the case where the point is */
int bucket_retour_key( pBucket bucket , pMesh mesh ,  pPoint point , double increment )
{
				
	int i,j,k,key,N,cherche,indice = 0 , ret ;
	pPoint  pp1,pp2;
  
	
	/* calcul i,j,k */
	i = max(0,(int)((bucket->size)*(point->c[0]))-1 + increment) ; 
	j = max(0,(int)((bucket->size)*(point->c[1]))-1 + increment) ;
	k = max(0,(int)((bucket->size)*(point->c[2]))-1 + increment) ;
	

		
	key = (j*bucket->size+k)*bucket->size+i ;

	/* we check if the tab head has a point in it */
	if ( bucket->head[key] )
	{
      	ret     = key;
      	cherche = bucket->head[key] ;
      
     
		/* we take the points of the same subdomaine */
		while( bucket->link[ cherche ])
		{
			cherche = bucket->link[cherche] ;			
			pp2 = &mesh->point[cherche];
		}		
	} 
	/* else there is no neighbourgh in this subdomain so we explore the subdomain around */
	else
	{
		increment ++ ;
		ret = bucket_retour_key( bucket , mesh , point , increment ) ; 
	}
	
	return ret ;	
}

/* free bucket */
void free_bucket (pBucket bucket)
{
	free( bucket->head );
	free( bucket->link );
}







